# ================================
# EMC Pipeline for Zhu et al. (2010) — Toy Example
# ================================

#--- 0. Dependencies ----------------------------------------------------------
if (!requireNamespace("MASS", quietly=TRUE)) install.packages("MASS")
if (!requireNamespace("truncnorm", quietly=TRUE)) install.packages("truncnorm")
if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly=TRUE)) install.packages("reshape2")


library(MASS)        # for mvrnorm
library(truncnorm)   # for rtruncnorm
library(ggplot2)
library(reshape2)

#--- 1. Simulate Zhu-style Toy Functional Data --------------------------------
simulate_sim1 <- function(n = 200,
                          T = 4,
                          df_basis = 10,
                          noise_sd = 0.05,
                          batch_sd = 0.5,
                          seed = 123,
                          plot = TRUE,
                          sample_ids = 1:4) {
  set.seed(seed)
  
  # 1) time‐grid and cosine basis on [0,1]
  tgrid <- seq(0, 1, length = T)
  phi <- matrix(NA, nrow = T, ncol = df_basis)
  phi[,1] <- 1
  for (k in 1:(df_basis-1)) {
    phi[,k+1] <- sqrt(2)*cos(k*pi*tgrid)
  }
  
  J <- 4  # four functional predictors
  
  # 2) basis‐scores (amps)
  amps <- array(rnorm(n*J*df_basis), dim = c(n, J, df_basis))
  
  # 3) batch effects
  batch    <- sample(1:2, n, replace = TRUE)
  beta_bch <- rnorm(2, 0, batch_sd)
  
  # 4) build full curves
  curves <- array(0, dim = c(n, J, T))
  for (i in 1:n) {
    for (j in 1:J) {
      curves[i,j,] <- amps[i,j,] %*% t(phi) +
        beta_bch[batch[i]] +
        rnorm(T, 0, noise_sd)
    }
  }
  
  # 5) flatten into X: n × (J×T)
  X <- matrix(NA, nrow = n, ncol = J * T)
  for (i in 1:n) {
    X[i, ] <- as.vector(curves[i, , ])
  }
  
  # 6) fixed‐effects design matrix S
  S1 <- runif(n)
  S2 <- rbinom(n,1,0.5)
  S  <- cbind(1, S1, S2)  # intercept + two covariates
  
  # 7) generate binary y: only functions 2 & 4 predictive
  gamma2 <- rep(1.5, df_basis)
  gamma4 <- rep(-1.0, df_basis)
  eta_fun <- sapply(1:n, function(i) {
    sum(amps[i,2,] * gamma2) + sum(amps[i,4,] * gamma4)
  })
  alpha  <- c(-0.5, 1.0, 0.8)
  eta    <- alpha[1] + alpha[2]*S1 + alpha[3]*S2 + eta_fun + beta_bch[batch]
  p      <- 1 / (1 + exp(-eta))
  y      <- rbinom(n,1,p)
  
  # 8) group index for EM–CMC (each block of T columns is one func)
  group_index <- rep(1:J, each = T)
  
  # 9) Plotting (optional)
  if (plot) {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    
    for (j in 1:J) {
      curve_matrix <- matrix(NA, nrow = length(sample_ids), ncol = T)
      for (i in 1:length(sample_ids)) {
        curve_matrix[i, ] <- curves[sample_ids[i], j, ]
      }
      matplot(tgrid, t(curve_matrix), type = 'l', lty = 1, col = 1:4,
              xlab = "Time", ylab = paste("Function", j),
              main = paste("Sample curves for Func", j))
    }
    
    # Plot mean curves by batch for each functional predictor
    for (j in 1:J) {
      mean_batch1 <- colMeans(curves[batch == 1, j, ])
      mean_batch2 <- colMeans(curves[batch == 2, j, ])
      plot(tgrid, mean_batch1, type = "l", col = "blue", ylim = range(c(mean_batch1, mean_batch2)),
           ylab = paste("Mean of Func", j), xlab = "Time",
           main = paste("Batch effect: Func", j))
      lines(tgrid, mean_batch2, col = "red")
      legend("topright", legend = c("Batch 1", "Batch 2"), col = c("blue", "red"), lty = 1)
    }
    
   
  }
  
  # 10) Return the result
  return(list(
    y           = y,
    S           = S,
    X           = X,
    curves      = curves,
    phi         = phi,
    tgrid       = tgrid,
    batch       = batch,
    group_index = group_index
  ))
}

#--- sanity‐check plot for first 5 curves of predictor 1 -----------------------
# dat <- simulate_sim1(n=50)
# matplot(dat$tgrid, t(dat$curves[1:5,1,]), type="l",
#         xlab="t", ylab="f1(t)", main="First 5 subjects, predictor 1")


#--- 2. Helper Functions ------------------------------------------------------

# 2.1 log marginal posterior of tau (approx)
log_posterior_tau <- function(Y, X, S, tau, sigma2_b, group_index,
                              nu1=10, nu0=0.01,
                              sigma2_0=100, sigma2_1=100) {
  n <- length(Y)
  p <- ncol(X); q <- ncol(S); J <- length(tau)
  d <- p/J
  
  D_tau <- diag(rep(sapply(1:J,
                           function(j) nu1*tau[j] + nu0*(1-tau[j])),
                    each=d))
  W      <- diag(p)
  Sigma0 <- D_tau %*% W %*% D_tau
  
  V_alpha_inv <- crossprod(S) + diag(1/sigma2_1, q)
  V_alpha     <- solve(V_alpha_inv)
  
  V0_inv <- diag(p)/sigma2_b + Sigma0/sigma2_0
  V0     <- solve(V0_inv)
  
  V <- diag(n) + S %*% V_alpha %*% t(S) + X %*% V0 %*% t(X)
  
  # proxy “Z” for marginal like
  Zm <- matrix((2*Y-1)*0.1, n, 1)
  lp  <- -0.5*( t(Zm) %*% solve(V,Zm) + determinant(V,TRUE)$modulus )
  
  log_prior <- sum(dbinom(tau,1,0.5,log=TRUE))
  as.numeric(lp + log_prior)
}

# 2.2 sample latent Z
sample_latent_Z <- function(Y, mu_z) {
  n <- length(Y); Z <- numeric(n)
  idx1 <- which(Y==1); idx0 <- which(Y==0)
  if(length(idx1)>0)
    Z[idx1] <- rtruncnorm(length(idx1),a=-Inf,b=0,mean=mu_z[idx1],sd=1)
  if(length(idx0)>0)
    Z[idx0] <- rtruncnorm(length(idx0),a=0,  b=Inf,mean=mu_z[idx0],sd=1)
  Z
}

# 2.3 MH update sigma2_b
update_sigma2b <- function(sigma2b, tau, Y, X, S, group_index,
                           delta=0.1, d1=2, d2=1) {
  prop <- exp(rnorm(1,log(sigma2b),delta))
  lp_new <- log_posterior_tau(Y,X,S,tau,prop,group_index)
  lp_old <- log_posterior_tau(Y,X,S,tau,sigma2b,group_index)
  log_pr <- (-(d1+1))*(log(prop)-log(sigma2b)) + d2*(1/prop - 1/sigma2b)
  if(log(runif(1)) < (lp_new-lp_old+log_pr)) prop else sigma2b
}

  # 2.4 Gibbs update α
update_alpha <- function(Z, X, S, bl, sigma2_1=100) {
  resid <- Z - X %*% bl
  V_inv <- crossprod(S) + diag(1/sigma2_1,ncol(S))
  V     <- solve(V_inv)
  m     <- V %*% crossprod(S,resid)
  as.vector(mvrnorm(1,m,V))
}

# 2.5 Gibbs update b_l
update_bl <- function(Z, X, S, alpha, b0, sigma2_b) {
  resid <- Z - S %*% alpha
  V     <- solve(crossprod(X) + diag(1/sigma2_b,ncol(X)))
  m     <- V %*% (crossprod(X,resid) + b0/sigma2_b)
  as.vector(mvrnorm(1,m,V))
}

# 2.6 Gibbs update b0
update_b0 <- function(bl, tau, sigma2_b, sigma2_0=100,
                      nu1=10, nu0=0.01, group_index) {
  p <- length(bl); J <- length(tau); d <- p/J
  D_tau <- diag(rep(sapply(1:J,
                           function(j) nu1*tau[j]+nu0*(1-tau[j])),
                    each=d))
  V0_inv<- diag(p)/sigma2_b + D_tau^2/sigma2_0
  V0    <- solve(V0_inv)
  m0    <- V0 %*% (bl/sigma2_b)
  as.vector(mvrnorm(1,m0,V0))
}


#--- 3. EMC Moves -----------------------------------------------------------

# 3.1 Mutation (switch/swap)
mutation_tau <- function(tau, Y, X, S, sigma2_b, temp, group_index, xi=0.5) {
  J <- length(tau); tnew <- tau
  if(runif(1)<xi) {
    pos <- sample(1:J,1); tnew[pos] <- 1-tnew[pos]
  } else {
    ones <- which(tau==1); zeros <- which(tau==0)
    if(length(ones)>0 && length(zeros)>0) {
      i1 <- sample(ones,1); i0 <- sample(zeros,1)
      tnew[i1] <- 0; tnew[i0] <- 1
    }
  }
  lp_new <- log_posterior_tau(Y,X,S,tnew,sigma2_b,group_index)
  lp_old <- log_posterior_tau(Y,X,S,tau, sigma2_b,group_index)
  if(log(runif(1)) < (lp_new-lp_old)/temp) tnew else tau
}

# 3.2 Crossover
crossover_tau <- function(ti,tj,Y,X,S,sigma2_b,ti_temp,tj_temp,group_index) {
  J  <- length(ti); cp <- sample(1:(J-1),1)
  new_i <- c(ti[1:cp], tj[(cp+1):J])
  new_j <- c(tj[1:cp], ti[(cp+1):J])
  lpi_new <- log_posterior_tau(Y,X,S,new_i,sigma2_b,group_index)
  lpi_old <- log_posterior_tau(Y,X,S,ti,   sigma2_b,group_index)
  lpj_new <- log_posterior_tau(Y,X,S,new_j,sigma2_b,group_index)
  lpj_old <- log_posterior_tau(Y,X,S,tj,   sigma2_b,group_index)
  if(log(runif(1)) < ( (lpi_new-lpi_old)/ti_temp + (lpj_new-lpj_old)/tj_temp ))
    list(tau_i=new_i, tau_j=new_j)
  else
    list(tau_i=ti,    tau_j=tj)
}

# 3.3 Exchange
exchange_tau <- function(ti,tj,Y,X,S,sigma2_b,ti_temp,tj_temp,group_index) {
  lpi <- log_posterior_tau(Y,X,S,ti,sigma2_b,group_index)
  lpj <- log_posterior_tau(Y,X,S,tj,sigma2_b,group_index)
  if(log(runif(1)) < ( (lpj-lpi)/ti_temp + (lpi-lpj)/tj_temp ))
    list(tau_i=tj, tau_j=ti)
  else
    list(tau_i=ti, tau_j=tj)
}


#--- 4. Main EMC -------------------------------------------------------------
run_emc <- function(Y, X, S, group_index,
                    max_iter=200, burn_in=1000,
                    n_chains=4, max_temp=10, init_ratio=1.2,
                    zeta=0.6, xi=0.5) {
  n <- length(Y); p <- ncol(X); q <- ncol(S); J <- length(unique(group_index))
  temps <- max_temp / init_ratio^(0:(n_chains-1))
  temps <- temps / min(temps)
  
  # init chains
  tau_ch  <- matrix(rbinom(J*n_chains,1,0.5), n_chains, J)
  sigma2b <- rep(1, n_chains)
  alpha   <- matrix(0, n_chains, q)
  b0      <- matrix(0, n_chains, p)
  bl      <- matrix(0, n_chains, p)
  
  # storage (cold chain)
  nsamp     <- max_iter - burn_in
  tau_samps <- matrix(0, nsamp, J)
  sig_samps <- numeric(nsamp)
  alp_samps <- matrix(0, nsamp, q)
  bl_samps  <- matrix(0, nsamp, p)
  
  for(iter in 1:max_iter) {
    cat(sprintf("current iteration %d\n", iter))
    # 1) latent Z + sigma2b
    for(c in 1:n_chains) {
      mu_z <- S %*% alpha[c,] + X %*% bl[c,]
      Z_lat <- sample_latent_Z(Y, mu_z)
      sigma2b[c] <- update_sigma2b(sigma2b[c], tau_ch[c,], Y, X, S, group_index)
    }
    # 2) mutation / crossover
    if(runif(1)<zeta) {
      for(c in 1:n_chains)
        tau_ch[c,] <- mutation_tau(tau_ch[c,], Y, X, S, sigma2b[c], temps[c], group_index, xi)
    } else {
      prs <- split(sample(1:n_chains), rep(1:(n_chains/2), each=2))
      for(pr in prs) {
        res <- crossover_tau(tau_ch[pr[1],], tau_ch[pr[2],],
                             Y, X, S, sigma2b[pr[1]],
                             temps[pr[1]], temps[pr[2]],
                             group_index)
        tau_ch[pr[1],] <- res$tau_i
        tau_ch[pr[2],] <- res$tau_j
      }
    }
    # 3) exchange
    for(i in 1:(n_chains-1)) {
      res <- exchange_tau(tau_ch[i,], tau_ch[i+1,],
                          Y, X, S, sigma2b[i],
                          temps[i], temps[i+1],
                          group_index)
      tau_ch[i,]   <- res$tau_i
      tau_ch[i+1,] <- res$tau_j
    }
    # 4) update α, b0, bl
    for(c in 1:n_chains) {
      mu_z <- S %*% alpha[c,] + X %*% bl[c,]
      Z_lat <- sample_latent_Z(Y, mu_z)
      alpha[c,] <- update_alpha(Z_lat, X, S, bl[c,])
      b0[c,]    <- update_b0(bl[c,], tau_ch[c,], sigma2b[c], group_index=group_index)
      bl[c,]    <- update_bl(Z_lat, X, S, alpha[c,], b0[c,], sigma2b[c])
    }
    # Store cold‐chain after burn‐in
    if(iter > burn_in) {
      idx <- iter - burn_in
      tau_samps[idx,] <- tau_ch[n_chains,]
      sig_samps[idx]  <- sigma2b[n_chains]
      alp_samps[idx,] <- -alpha[n_chains,]
      bl_samps[idx,]  <- -bl[n_chains,]
    }
  }
  
  list(
    tau_samples    = tau_samps,
    sigma2b_samples= sig_samps,
    alpha_samples  = alp_samps,
    bl_samples     = bl_samps,
    tau_post       = colMeans(tau_samps),
    tau_map        = as.numeric(colMeans(tau_samps)>0.5)
  )
}

#--- 5. Run & Summarize -------------------------------------------------------
set.seed(42)
dat <- simulate_sim1(n=1000, T=4, df_basis=10, noise_sd=0.1, batch_sd=0.5,plot = TRUE,sample_ids = c(1, 2, 3, 4,5))
res <- run_emc(
  Y           = dat$y,
  X           = dat$X,
  S           = dat$S,
  group_index = dat$group_index,
  max_iter    = 1000,
  burn_in     = 100,
  n_chains    =   4,
  max_temp    =  10,
  init_ratio  = 1.2
)

par(mfrow=c(1,1))
# Posterior inclusion probabilities
barplot(res$tau_post,
        main = "Posterior Inclusion Probabilities",
        xlab = "Function Index", ylab = "P(τ_j = 1)",
        ylim = c(0,1), col = "steelblue")
abline(h=0.5, lty=2, col="red")

# MAP‐selected functions
cat("MAP‐selected functions:", which(res$tau_map==1), "\n")

# ================================
# Evaluation for Sim1 (n = 1000, no train-test split)
# ================================

# Get posterior means from cold chain
bl_mean     <- colMeans(res$bl_samples)
alpha_mean  <- colMeans(res$alpha_samples)

# Compute predicted probabilities
eta <- as.vector(dat$S %*% alpha_mean + dat$X %*% bl_mean)
prob <- 1 / (1 + exp(-eta))
yhat <- ifelse(prob > 0.5, 1, 0)

# Ground truth
ytrue <- dat$y

# Evaluation metrics
sensitivity <- sum(yhat == 1 & ytrue == 1) / sum(ytrue == 1)
specificity <- sum(yhat == 0 & ytrue == 0) / sum(ytrue == 0)
misclass_rate <- mean(yhat != ytrue)
accuracy <- mean(yhat == ytrue)

# Report
cat(sprintf("Evaluation Metrics on Full Dataset (n = 1000):\n"))
cat(sprintf("  Sensitivity:           %.2f%%\n", 100 * sensitivity))
cat(sprintf("  Specificity:           %.2f%%\n", 100 * specificity))
cat(sprintf("  Accuracy:              %.2f%%\n", 100 * accuracy))
cat(sprintf("  Misclassification Rate: %.2f%%\n", 100 * misclass_rate))




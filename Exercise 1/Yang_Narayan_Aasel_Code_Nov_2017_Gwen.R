
# This is a replication of Yang, Narayan and Assael (2006, Marketing Science)
# Author: Gwen Ahn (under Puneet Manchanda's guidance)

# load libraries ----------------------------------------------------------

library(bayesm)
library(coda)
library(MCMCpack)
library(matrixcalc)
library(mvtnorm)   
library(truncnorm)
library(clusterGeneration)


# set parameters ----------------------------------------------------------

set.seed(100)

I = 100 # number of households i = 1, ...., I
J = 50 # number of programs (observations/household)
ncovar <- 2 # 2 covariates 

Sigma = matrix(c(1.0, 0.5, 0.5, 1.0), 2, 2)

# beta
mean_beta = c(1,1,1,1) # set beta_bar = 1 for beta_0_h, beta_1_h, beta_0_w, beta_1_w
Psi = matrix(rep(0.05, length(mean_beta) * length(mean_beta)), length(mean_beta)) + 
  diag(0.05, length(mean_beta))
beta_household = rmvnorm(100, mean_beta, Psi) # individual household specific beta

# structural parameters
w = c(w_hw = 0.5, w_wh = 0.8) # interdependence of preferences

# reduced form coefficients
Pi <- vector()
Pi[1] = 1/(1-w[1] * w[2])
Pi[2] = w[1]/(1-w[1] * w[2])
Pi[3] = w[2]/(1-w[1] * w[2])
Pi[4] = 1/(1-w[1] * w[2])

# reduced form variance
eye_w = diag(1, 2) - matrix(c(0, w, 0), 2, 2, byrow = T)
Omega = solve(eye_w) %*% Sigma %*% t(solve(eye_w))


# simulate the data -------------------------------------------------------

# generate y based on X values and parameters in each 
sim_data <- NULL # list to contain dataframes of each households
prop_neg <- NULL
for(i in 1:I){
  # generate X (intercept and x1) for husband/wife
  # changing mean and variance of normal affects proportion of 0s that we observe
  x_husb = cbind(rep(1, J), rnorm(J, -1, 1))
  x_wife = cbind(rep(1, J), rnorm(J, -1, 1))
  x = cbind(x_husb, x_wife)
  
  # x_h'beta_h, x_w'beta_w
  husb_xb <- x[ , 1:2] %*% (beta_household[i, 1:2]) 
  wife_xb <- x[ , 3:4] %*% (beta_household[i, 3:4])
  
  # deterministic part of latent utility
  husb_y <- Pi[1] * husb_xb + Pi[2] * wife_xb 
  wife_y <- Pi[3] * husb_xb + Pi[4] * wife_xb
  
  # generate latent utility 
  y_latent <- NULL
  for(j in 1:J){
    y_latent <- rbind(y_latent, mvrnorm(1, cbind(husb_y[j], wife_y[j]), Omega))
  }
  
  # observed y_husb and y_wife
  y_husb <- ifelse(y_latent[ ,1] > 0, y_latent[ ,1], 0) 
  y_wife <- ifelse(y_latent[ ,2] > 0, y_latent[ ,2], 0)
  sim_data[[i]] <- list(y_husb = y_husb, y_wife = y_wife, 
                        x_husb = x_husb, x_wife = x_wife)
  
  prop_neg <- rbind(prop_neg, c(sum(y_husb <= 0)/50, sum(y_wife <= 0)/50))
  # print(i)
}

apply(prop_neg, 2, mean)
# [1] 0.2154 0.2046 when x_husb/x_wife = cbind(rep(1, J), rnorm(J, 0, 1))
# [1] 0.5198 0.5198 when x_husb/x_wife = cbind(rep(1, J), rnorm(J, -1, 1))

# Set initial values and create dataframes to hold the draws --------------

niter <- 2000 
set.seed(100)

# beta_bar
beta_bar_draws <- matrix(vector(), nrow = niter, ncol = ncovar * 2)
beta_current <- mean_beta # beta_bar 


# individual beta_i's for each household (husband0, husband1, ..., husband_ncovar, wife0, ...)
beta_hh_draws <- array(NA, dim = c(niter, 100, ncovar*2))
beta_i_current <- matrix(rep(1, I * ncovar * 2), ncol = ncovar * 2, byrow = T)

# Psi (variance of beta_i's)
Psi_draws <- matrix(vector(), nrow = niter, ncol = (2*ncovar) * (2*ncovar))
Psi_current <- riwish(10, rcorrmatrix(4))


# structural parameters
w_draws <- matrix(vector(), nrow = niter, ncol = 2)
colnames(w_draws) <- c("w_hw", "w_wh")
w_current <- runif(2) 

# transformed coefficients (in the reduced form)
Pi_draws <- matrix(vector(), nrow = niter, ncol = 4)
colnames(Pi_draws) <- c("Pi11", "Pi12", "Pi21", "Pi22")

Pi_current <- vector()
Pi_current[1] <- 1/(1-w_current[1]*w_current[2])
Pi_current[2] <- w_current[1]/(1-w_current[1]*w_current[2])
Pi_current[3] <- w_current[2]/(1-w_current[1]*w_current[2])
Pi_current[4] <- 1/(1-w_current[1]*w_current[2])


# variance Sigma (from structural equations)
Sigma_draws <- matrix(vector(), nrow = niter, ncol = 4)
colnames(Sigma_draws) <- c("sigma_11", "sigma_12", "sigma_21", "sigma_22")
Sigma_current <- riwish(5, matrix(c(1,.3, .3, 1), 2, 2))

# transformed variance Omega (in the reduced form)
eye_w_current <- diag(1, 2) - matrix(c(0, w_current, 0), nrow = 2, byrow = T)
Omega_current <- solve(eye_w_current) %*% Sigma_current %*% t(solve(eye_w_current))
Omega_draws <- matrix(vector(), nrow = niter, ncol = ncovar * ncovar)
colnames(Omega_draws) <- c("omega_11", "omega_12", "omega_21", "omega_22")


# Make the draws ----------------------------------------------------------

for(r in 1:niter){
  
  quadratic.d <- 0
  quadratic.dy <- 0
  Omega_1st_term <- 0
  
  for (i in 1:I){
    y_husb <- sim_data[[i]]$y_husb
    y_wife <- sim_data[[i]]$y_wife
    x_husb <- sim_data[[i]]$x_husb
    x_wife <- sim_data[[i]]$x_wife
    nobs <- length(y_husb)
  
    # data augmentation step ----------------------------------------------
    # first generate the latent quantities y_h_latent, y_w_latent  
    # draw from conditional normal distribution
    husb_xb = x_husb %*% (beta_i_current[i, 1:ncovar]) 
    wife_xb = x_wife %*% (beta_i_current[i, (ncovar +1):(2*ncovar)])
    
    # latent mean of (y_h*, y_w*)
    y_latent_h_hat <- Pi_current[1] * husb_xb + Pi_current[2] * wife_xb
    y_latent_w_hat <- Pi_current[3] * husb_xb + Pi_current[4] * wife_xb
    
    # how do we condition if y_husb < 0 & y_wife < 0
    
    y_h_latent_bar <- y_latent_h_hat +  Omega_current[1,2]/Omega_current[2,2] * (y_wife - y_latent_w_hat)
    y_w_latent_bar <- y_latent_w_hat +  Omega_current[2,1]/Omega_current[1,1] * (y_husb - y_latent_h_hat)
    
    cond_omega_husb <- Omega_current[1, 1] - (Omega_current[2,1] * Omega_current[1,2]/Omega_current[2,2])
    cond_omega_wife <- Omega_current[2, 2] - (Omega_current[1,2] * Omega_current[2,1]/Omega_current[1,1])
    
    # generate draws from conditional truncated normal + latent variables
    # generate from (-10, 0) of the normal|spouse's latent utility
    neg_husb <- rtruncnorm(1, rep(-10, nobs), rep(0, nobs), 
                           y_h_latent_bar, rep(sqrt(cond_omega_husb), nobs))
    neg_wife <-  rtruncnorm(1, rep(-10, nobs), rep(0, nobs), 
                            y_w_latent_bar, rep(sqrt(cond_omega_wife), nobs))
    
    # latent husband/wife utility
    y_h_latent <- ifelse(y_husb > 0, y_husb, neg_husb)
    y_w_latent <- ifelse(y_wife > 0, y_wife, neg_wife)
    
    # drawing for parameters ----------------------------------------------
    # drawing for structural parameters w_hw, w_wh (1)
    inv.omega.current <- solve(Omega_current)
    
    for(j in 1:J){
      d <- matrix(c(husb_xb[j], wife_xb[j], 0, 0,  
                    0, 0, husb_xb[j], wife_xb[j]), nrow = 2, byrow = T)
      quadratic.d <- quadratic.d + t(d) %*% inv.omega.current %*% d
      quadratic.dy <- quadratic.dy + t(d) %*% inv.omega.current %*% c(y_h_latent[j], y_w_latent[j])
    }
    
    # drawing for Sigma (1)
    Omega_1st_term <- Omega_1st_term + 
      t(cbind(y_h_latent, y_w_latent) - cbind(y_latent_h_hat, y_latent_w_hat)) %*% 
      (cbind(y_h_latent, y_w_latent) - cbind(y_latent_h_hat, y_latent_w_hat))
    
    # draw for beta_i 
    inv.Sigma.current <- solve(Sigma_current)
    
    S_beta_i_1st <- 0
    M_beta_i_1st <- 0
    
    for(j in 1:J){
      x <- matrix(c(x_husb[j, ], 0, 0, 0, 0, x_wife[j, ]), nrow = 2, byrow = T)
      S_beta_i_1st <- S_beta_i_1st + t(x) %*% inv.Sigma.current %*% x
      M_beta_i_1st <- M_beta_i_1st + t(x) %*% inv.Sigma.current %*% (eye_w_current) %*% as.matrix(c(y_h_latent[j], y_w_latent[j]), nrow = 2)
    }
    
      S_beta_i <- solve((S_beta_i_1st + solve(Psi_current)))
      M_beta_i <- S_beta_i %*% (M_beta_i_1st + solve(Psi_current) %*% beta_current)
      
      beta_i_current[i, ] <- rmvnorm(1,M_beta_i,S_beta_i)
    }
  
  ## drawing for structural parameters w_hw, w_wh (2)
  S <- solve(quadratic.d + diag(1,4))
  U <- S %*% (quadratic.dy)
  
  # drawing from Pi
  temp <- rmvnorm(1,U,S); # update Pi
  w_hw =  temp[2]/temp[1]
  w_wh = temp[3]/temp[1]
  w_current <- c(w_hw, w_wh) 
  Pi_current[1] <- 1/(1-w_hw * w_wh)
  Pi_current[2] <- w_hw/(1-w_hw * w_wh)
  Pi_current[3] <- w_wh/(1-w_hw * w_wh)
  Pi_current[4] <- 1/(1-w_hw * w_wh)
  
  ## draw Sigma (2)
  eye_w_current <-  diag(1,2) - matrix(c(0, w_hw, w_wh, 0), 2, 2, byrow = TRUE)
  Omega_current <- riwish(I*J+10, Omega_1st_term + diag(1, 2)) # instead of 10 suggested
  Sigma_current <- eye_w_current %*% Omega_current %*% t(eye_w_current)
  
  ## draw beta_bar ##
  beta_bar_S <- solve(solve(Psi_current/I) +  diag(0.01, 2*ncovar))
  sum_beta_i <- apply(beta_i_current, 2, sum)
  beta_bar_mean_U <- beta_bar_S %*% (solve(Psi_current) %*% sum_beta_i)
  beta_current <- matrix(rmvnorm(1, beta_bar_mean_U, beta_bar_S), 4, 1)
  
  ## draw V_beta ##
  diff_beta_i_mean <- beta_i_current - matrix(rep(beta_current, I), nrow = I, byrow = T)
  V_a <- t(diff_beta_i_mean) %*% diff_beta_i_mean + diag(1, nrow(Psi_current)) #instead of 10 suggested
  V_b <- I + 10
  Psi_current <- riwish(V_b,V_a)
  
  
  beta_bar_draws[r, ] <- beta_current
  Omega_draws[r, ] <- as.vector(Omega_current)
  Sigma_draws[r, ] <- as.vector(Sigma_current)
  Pi_draws[r , ] <- Pi_current
  w_draws[r, ] <- w_current
  Psi_draws[r, ] <- as.vector(Psi_current)
  beta_hh_draws[r, ,] <- beta_i_current

  print(r)
}


# Summarize the draws -----------------------------------------------------


beta_bar_draws <- as.mcmc(beta_bar_draws)
Psi_draws <- as.mcmc(Psi_draws)
w_draws <- as.mcmc(w_draws)
Pi_draws <- as.mcmc(Pi_draws)
Sigma_draws <- as.mcmc(Sigma_draws)
Oemga_draws <- as.mcmc(Omega_draws)

summary(beta_bar_draws)
quartz()
plot(beta_bar_draws)
# effectiveSize(beta_bar_draws)
acf(beta_bar_draws)

# summary(Psi_draws)
# plot(Psi_draws)
# acf(Psi_draws)

summary(w_draws)
par(mar=c(2,2,2,2))
plot(w_draws)
effectiveSize(w_draws)
acf(w_draws)

# Pi
# summary(Pi_draws)
# plot(Pi_draws)

summary(Sigma_draws)
plot(Sigma_draws)


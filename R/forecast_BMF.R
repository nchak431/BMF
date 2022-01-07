source("mat_gen_BMF.R")

library(matrixcalc)
library(clusterGeneration)
library(fBasics)
library(scoringRules)
library("CVglasso")
#==========input=============
k1 <- 3
k2 <- 30
nobs <- 100
theta_true <- 0.5
iteration <- 3000
burn <- 1000
repl <- 1
lag <- 1
spec_rad<- 0.7
edge_density <- 0.04
horizon <- 10
N <- 1000
#============================================================================
vec <- permutation(k1,k2)[[1]]  
vec1 <- permutation(k1,k2)[[2]]
vec2 <- permutation(k1,k2)[[3]]

p <- (3*k1) + k2
p1 <- (k1+k2)
#===========Generate VAR transition matrix================
A <- Gen_A(k1,k2,theta_true,lag, spec_rad, edge_density)[[1]]
A_block_one <- Gen_A(k1,k2,theta_true,lag, spec_rad, edge_density)[[2]]
A11_true <- (A_block_one[1:k1,1:k1]/(theta_true^2))
A12_true <- (A_block_one[1:k1,(k1+1):(k1+k2)]/(theta_true^2))
A21_true <- (A_block_one[(k1+1):(k1+k2),1:k1])
A22_true <- (A_block_one[(k1+1):(k1+k2),(k1+1):(k1+k2)])
A1 <- cbind(A11_true, A12_true)
A2 <- cbind(A21_true, A22_true)
A_true <- rbind(A1,A2) # true A
sigma_true <- Gen_Sigma (k1,k2)[[1]]
desired_SNR <- 2
k <- norm(A, type="F")/norm(desired_SNR*sigma_true, type="F")
sigma_true <- k * sigma_true

forecast_med <- array(0, c(p, horizon, repl))
ferr_med <- array(0, c(horizon,1, repl))
ferr_med_low <- array(0, c(horizon,1, repl))

  r=1
  W_sim  <- array(0, c(p, p, (iteration-burn)))
  theta_sim  <- array(0, c(1, 1, (iteration-burn)))
  sigma_sim <- array(0, c(p, p, (iteration-burn)))
  b1_sim <- array(0, c((k1*k1), 1, (iteration-burn)))
  b2_sim <- array(0, c((k1*k2), 1, (iteration-burn)))
  b3_sim <- array(0, c((k1*k2), 1, (iteration-burn)))
  b4_sim <- array(0, c((k2*k2), 1, (iteration-burn)))
  
  w_final <- matrix(0,p,p)
  theta_final <- 0
  sigma_final <- matrix(0,p,p)
  b1_final <- rep(0,(k1*k1))
  b2_final <- rep(0,(k1*k2))
  b3_final <- rep(0,(k1*k2))
  b4_final <- rep(0,(k2*k2))

  count <- 0
  #===============Generate data================
  y_true <- list()
  y <- list()
  y_mod <- list()
  err <- list()
  y_true[[1]] <- mvrnorm(1,rep(0,p),sigma_true)
  y[[1]]<-y_true[[1]][vec]
  y_mod[[1]]<- y[[1]][vec1]
  for ( t in 2:(nobs + horizon + 1))
  {
    err[[t]] <- mvrnorm(1,rep(0,p),sigma_true)
    y_true[[t]] <- as.vector((A %*% y_true[[t-1]]) + err[[t]])
    y[[t]]<-y_true[[t]][vec]
    y_mod[[t]]<-y[[t]][vec1]
    
  }
  
  y_full <- as.matrix(do.call(cbind,y))
  y_mod_full <- do.call(cbind,y_mod)
  data <-y_full[,1:(nobs + 1)]
  data_mod <-as.matrix(y_mod_full)[,1:(nobs + 1)]
  data_t <- t(data)
  forecast_true <- y_full[,(nobs + 1 + 1):(nobs + 1 + horizon)]
    
  X <- data_t[1:nobs,]
  Y <- data_t[2:(nobs + 1),]
  LSE <- lm(Y ~ X - 1)$coef
  #======Initial values===========
  A11_initial <- t(LSE)[1:k1,1:k1]
  b1<- as.vector(vec(A11_initial))
  b2 <- numeric((k1*k2))
  b3 <- numeric((k1*k2))
  
  fit_cv <- cv.glmnet(X, Y, alpha=1, family = "mgaussian", nfolds = 5, type.measure = "mse")
  lambda_min <- fit_cv$lambda.min
  fit <- glmnet(X, Y, lambda=lambda_min, family = "mgaussian")
  sparse_coef <- do.call(cbind, fit$beta)
  Phi_initial <- as.matrix(sparse_coef)
  A22_initial <- t(Phi_initial)[(k1+1):(k1+k2),(k1+1):(k1+k2)]
  b4 <- as.vector(vec(A22_initial))
  #=====================================
  source("functions.R")
  
  
  theta_func <- theta_matrix (theta_true,b1,b2,b3,b4,k1,k2,vec1)
  theta_mat <- theta_func[[1]]
  u <- theta_func[[2]]
  v <- theta_func[[3]]
  w11 <- theta_func[[4]]
  w12 <- theta_func[[5]]
  w21 <- theta_func[[6]]
  w22 <- theta_func[[7]]
  w <- theta_func[[8]]
  w_mod <- theta_func[[9]]
  
  S_mat <- S_matrix(data,nobs,k1,k2)
  s1 <-S_mat[[1]]
  s2 <- S_mat[[2]]
  s3 <- S_mat[[3]]
  s11<- S_mat[[4]]
  s12 <- S_mat[[5]]
  s21 <- S_mat[[6]]
  s22 <- S_mat[[7]]
  ss11<- S_mat[[8]]
  ss12<- S_mat[[9]]
  ss21 <- S_mat[[10]]
  ss22 <- S_mat[[11]]
  s <- S_mat[[12]]
  #============Initial values of Sigma===================
  omega <-list()
  sigma_adj <- numeric()
  
  for(i in 1:k1)
  {
    sigma1 <- 0.3*diag(3)
    omega[[i]] <- solve(sigma1)
    sigma_adj[i] <- 1/(t(u)%*% solve(sigma1)%*% u)  
  }
  
  sigma2 <- rep(1,k2)
  #==============hyperparameters=================
  alpha <- 1
  beta <- 1
  Q <- 0.3*diag(3)
  df <- 1
  
  sigma_inv_func <- sigma_inverse(omega,sigma2,k1,k2)
  sigma_11 <- sigma_inv_func[[1]]
  sigma_22 <- sigma_inv_func[[2]]
  sigma_inv<- sigma_inv_func[[3]]
  Final_Sigma_Inv <- sigma_inv_func[[4]]
  sigma_mat <- sigma_inv_func[[5]]
  
  A11_func <- para_A11(theta_mat,s,sigma_inv,sigma_11,s11,ss11,s12,w12)
  gamma_A11 <- A11_func[[1]]
  d_A11 <- A11_func[[2]]
  q11 <- 0
  tao11 <- 0.5
  
  A12_func <- para_A12(k1,k2,u,sigma_inv,sigma_11,s12,ss12,s22,w11)
  gamma_A12 <- A12_func[[1]]
  d_A12 <- A12_func[[2]]
  q12 <- 0.96
  tao12 <- 0.5
  
  A21_func <- para_A21(k1,k2,v,sigma_22,s,s11,ss21,s12,w22)
  gamma_A21 <- A21_func[[1]]
  d_A21 <- A21_func[[2]]
  q21 <- 0.96
  tao21 <- 0.5
  
  A22_func <- para_A22(v,sigma_22,s22,ss22,s12,w21)
  gamma_A22 <- A22_func[[1]]
  d_A22 <- A22_func[[2]]
  q22 <- 0.7
  tao22 <- 0.5
  
  #============Gibbs sampler============
  for (i in 1:iteration)
  {
    for(j1 in 1:k1)
    {
      for(j in (((j1-1)*k1) + 1):(((j1-1)*k1) + k1))
      {
        a1 <- sum(b1*gamma_A11[,j])-(b1[j]*gamma_A11[j,j])-d_A11[j]
        a2 <- gamma_A11[j,j] +(1/((tao11^2)*sigma_adj[j-((j1-1)*k1)]))
        pr <- q11/(q11+(((1-q11)*exp(a1^2/(2*a2)))/(tao11*(sqrt(a2*(sigma_adj[j-((j1-1)*k1)]))))))
        x1 <- -(a1/a2)
        x2 <- sqrt(1/a2)
        b1[j] <- draw_A11(x1,x2,pr,j)
        
      }
    }
    
    for(l1 in 1:k2)
    {
      for(l in (((l1-1)*k1) + 1):(((l1-1)*k1) + k1))
      {
        a1 <- sum(b2*gamma_A12[,l])-(b2[l]*gamma_A12[l,l])-d_A12[l]
        a2 <- gamma_A12[l,l] +(1/((tao12^2)*sigma_adj[l-((l1-1)*k1)]))
        pr <- q12/(q12+(((1-q12)*exp(a1^2/(2*a2)))/(tao12*(sqrt(a2*(sigma_adj[l-((l1-1)*k1)]))))))
        x1 <- -(a1/a2)
        x2 <- sqrt(1/a2)
        b2[l] <- draw_A12(x1,x2,pr,l)
        
      }
    }
    
    for(m1 in 1:k1)
    {
      for(m in (((m1-1)*k2) + 1):(((m1-1)*k2) + k2))
      {
        a1 <- sum(b3*gamma_A21[,m])-(b3[m]*gamma_A21[m,m])-d_A21[m]
        a2 <- gamma_A21[m,m] +(1/((tao21^2)*sigma2[m-((m1-1)*k2)]))
        pr <- q21/(q21+(((1-q21)*exp(a1^2/(2*a2)))/(tao21*(sqrt(a2*(sigma2[m-((m1-1)*k2)]))))))
        x1 <- -(a1/a2)
        x2 <- sqrt(1/a2)
        b3[m] <- draw_A21(x1,x2,pr,m)
        
      }
    }
    
    
    for(n1 in 1:k2)
    {
      for(n in (((n1-1)*k2) + 1):(((n1-1)*k2) + k2))
      {
        a1 <- sum(b4*gamma_A22[,n])-(b4[n]*gamma_A22[n,n])-d_A22[n]
        a2 <- gamma_A22[n,n] + (1/((tao22^2)*sigma2[n-((n1-1)*k2)]))
        pr <- q22/(q22+(((1-q22)*exp(a1^2/(2*a2)))/(tao22*(sqrt(a2*(sigma2[n-((n1-1)*k2)]))))))
        x1 <- -(a1/a2)
        x2 <- sqrt(1/a2)
        b4[n] <- draw_A22(x1,x2,pr,n)
        
      }
    }
    
    d_theta <- theta_dist_coeff(s1,s2,s3,Final_Sigma_Inv,b1,b2,b3,b4,k1,k2,vec1)
    theta <- draw_theta(N,d_theta)
    theta_func <- theta_matrix (theta,b1,b2,b3,b4,k1,k2,vec1)
    theta_mat <- theta_func[[1]]
    u <- theta_func[[2]]
    v <- theta_func[[3]]
    w11 <- theta_func[[4]]
    w12 <- theta_func[[5]]
    w21 <- theta_func[[6]]
    w22 <- theta_func[[7]]
    w <- theta_func[[8]]
    w_mod <- theta_func[[9]]
    A11_mat <- matrix(b1,nrow=k1,ncol=k1,byrow=F)
    A12_mat <- matrix(b2,nrow=k1,ncol=k2,byrow=F)
    A_upper <- cbind(A11_mat,A12_mat)
    A21_mat <- matrix(b3,nrow=k2,ncol=k1,byrow=F)
    A22_mat <- matrix(b4,nrow=k2,ncol=k2,byrow=F)
    A_lower <- cbind(A21_mat,A22_mat)
    sigma_func <- parameter_sigma (k1,k2,nobs,u,Q,df,alpha,beta,data_mod,w_mod,A21_mat,A22_mat,tao11,tao21,tao22,A_upper,A_lower)
    omega <- sigma_func[[1]]
    sigma2 <- sigma_func[[2]]
    sigma_adj <- sigma_func[[3]]
    sigma_inv_func <- sigma_inverse(omega,sigma2,k1,k2)
    sigma_11 <- sigma_inv_func[[1]]
    sigma_22 <- sigma_inv_func[[2]]
    sigma_inv<- sigma_inv_func[[3]]
    Final_Sigma_Inv <- sigma_inv_func[[4]]
    sigma_mat <- sigma_inv_func[[5]]
    
    A11_func <- para_A11(theta_mat,s,sigma_inv,sigma_11,s11,ss11,s12,w12)
    gamma_A11 <- A11_func[[1]]
    d_A11 <- A11_func[[2]]
    
    A12_func <- para_A12(k1,k2,u,sigma_inv,sigma_11,s12,ss12,s22,w11)
    gamma_A12 <- A12_func[[1]]
    d_A12 <- A12_func[[2]]
    
    A21_func <- para_A21(k1,k2,v,sigma_22,s,s11,ss21,s12,w22)
    gamma_A21 <- A21_func[[1]]
    d_A21 <- A21_func[[2]]
    
    A22_func <- para_A22(v,sigma_22,s22,ss22,s12,w21)
    gamma_A22 <- A22_func[[1]]
    d_A22 <- A22_func[[2]]
    
  if (i > burn) 
    {
      count <- count + 1
      W_sim [, , count] <- w
      theta_sim [ , , count] <- theta
      sigma_sim [ , , count] <- as.matrix(sigma_mat)
      b1_sim [ , , count] <- b1
      b2_sim [ , , count] <- b2
      b3_sim [ , , count] <- b3
      b4_sim [ , , count] <- b4
      
      w_final <- w_final + w
      theta_final <- theta_final + theta
      sigma_final <- sigma_final + sigma_mat
      b1_final <- b1_final + b1
      b2_final <- b2_final + b2
      b3_final <- b3_final + b3
      b4_final <- b4_final + b4
    }
    
  }
  
  w_final = w_final/count
  theta_est = theta_final/count
  sigma_final = sigma_final/count
  b1_final = b1_final/count
  b2_final = b2_final/count
  b3_final = b3_final/count
  b4_final = b4_final/count    
w_est <- apply(W_sim, c(1,2), mean)
sigma_est=apply(sigma_sim, c(1,2), mean)

#=========================================
true_b1 <- vec(A_block_one[1:k1,1:k1]/(theta_true^2))
true_b2 <- vec(A_block_one[1:k1,(k1+1):(k1+k2)]/(theta_true^2))
true_b3 <- vec(A_block_one[(k1+1):(k1+k2),1:k1])
true_b4 <- vec(A_block_one[(k1+1):(k1+k2),(k1+1):(k1+k2)])

zero <- function(x)
  {
    z <- length(which(x==0))  
    if (z > (length(x)/2))
    {
      p <- 0
    }else{
      p <- sum(x)/(length(x) - length(which(x==0)))  
    }
    p
  }
  
  corrected_b1 <- apply(b1_sim[ , , ], 1,  FUN=function(x) zero(x) )
  corrected_b2 <- apply(b2_sim[ , , ], 1,  FUN=function(x) zero(x) )
  corrected_b3 <- apply(b3_sim[ , , ], 1,  FUN=function(x) zero(x) )
  corrected_b4 <- apply(b4_sim[ , , ], 1,  FUN=function(x) zero(x) )
A11_est <- matrix(corrected_b1,nrow=k1,ncol=k1,byrow=F)
A12_est <- matrix(corrected_b2,nrow=k1,ncol=k2,byrow=F)
A21_est <- matrix(corrected_b3,nrow=k2,ncol=k1,byrow=F)
A22_est <- matrix(corrected_b4,nrow=k2,ncol=k2,byrow=F)
A_upper <- cbind(A11_est, A12_est)
A_lower <- cbind(A21_est, A22_est)
A_est <- rbind(A_upper, A_lower)
#====================model selection performance=====================
sensitivity_b1 <- metric(true_b1,corrected_b1)[[1]]
specificity_b1 <- metric(true_b1,corrected_b1)[[2]]
sensitivity_b2 <- metric(true_b2,corrected_b2)[[1]]
specificity_b2 <- metric(true_b2,corrected_b2)[[2]]
sensitivity_b3 <- metric(true_b3,corrected_b3)[[1]]
specificity_b3 <- metric(true_b3,corrected_b3)[[2]]
sensitivity_b4 <- metric(true_b4,corrected_b4)[[1]]
specificity_b4 <- metric(true_b4,corrected_b4)[[2]]
sensitivity_A_mat <- metric_A(A_block_one,A_est)[[1]]
specificity_A_mat <- metric_A(A_block_one,A_est)[[2]]
#==================Relative error================
diff_A <- A_est - A_true
rel_error_A <- norm(diff_A,'f')/norm(A_true,'f')

diff_A11 <- A11_est - A11_true
rel_error_A11 <- norm(diff_A11,'f')/norm(A11_true,'f')

diff_A12 <- A12_est - A12_true
rel_error_A12 <- norm(diff_A12,'f')/norm(A12_true,'f')

diff_A21 <- A21_est - A21_true
rel_error_A21 <- norm(diff_A21,'f')/norm(A21_true,'f')

diff_A22 <- A22_est - A22_true
rel_error_A22 <- norm(diff_A22,'f')/norm(A22_true,'f')
#Rel error of W
w1 <- A[vec,]
w_true <-w1[,vec] 
diff_w <- w_est-w_true
rel_error_w <- norm(diff_w,'f')/norm(w_true,'f')
#Rel error of sigma
temp1 <- sigma_true[vec,]
temp2 <- temp1[,vec]
temp3 <- temp2[vec1,]
true_sigma1 <- as.matrix(temp3[,vec1])

diff_sigma <- sigma_est - true_sigma1
rel_error_sigma <- norm(diff_sigma,'f')/norm(true_sigma1,'f')
#========================================
pred <- matrix(0,nrow=((3*k1)+k2),ncol=nobs+1)
err <- matrix(0,nrow=((3*k1)+k2),ncol=nobs+1)
for ( t in 2:(nobs+1))
{
  pred[,t] <- as.vector((w_est %*% y_full[,t-1]))
  err[,t] <- y_full[,t] - pred[,t] 
  
}
#Estimating error-covariance using graphical lasso
  lasso <- CVglasso(X = t(err), S = NULL, nlam = 10, lam.min.ratio = 0.01,
                    lam = NULL, diagonal = FALSE, path = FALSE, tol = 1e-04,
                    maxit = 10000, adjmaxit = NULL, K = 5, crit.cv = c("loglik", "AIC",
                                                                       "BIC"), start = c("warm", "cold"), cores = 1)
  
  cov_mat <- lasso$Sigma
  pred_sim <- array(0, c(p, horizon, (iteration-burn)))
  for(i in 1:(iteration-burn)){
  pred_y <- list()
  pred_y[[1]] <- mvrnorm(1,(W_sim[,,i] %*% y_full[,nobs+1]),cov_mat)
  for ( t in 2:horizon)
  {
    pred_y[[t]] <- mvrnorm(1,(W_sim[,,i] %*% pred_y[[t-1]]),cov_mat)
    
  }
  pred_y <- as.matrix(do.call(cbind,pred_y))
  
  pred_sim [ , , i] <- pred_y
}
mean <- apply(pred_sim,c(1,2),mean)
sd <- apply(pred_sim,c(1,2),sd)
# Obtaining CRPS and LPS
crps_mat <- matrix(0,nrow = ((3*k1)+k2),ncol=horizon)
logs_mat <- matrix(0,nrow = ((3*k1)+k2),ncol=horizon)
for ( i in 1:horizon)
{
  crps_mat[,i] <- crps(y = forecast_true[,i], family = "normal", mean = mean[,i], sd = sd[,i])
  logs_mat[,i] <- logs(y = forecast_true[,i], family = "normal", mean = mean[,i], sd = sd[,i])
}

#=======================forecast errors==========================
  forecast <- forecast_cred(pred_sim,forecast_true,horizon)
  forecast_med[,,r] <- forecast$forecast_med
  ferr_med[,,r] <- forecast$ferr_med
  ferr_med_low[,,r] <- forecast$ferr_med_low # forecast error for quarterly variables

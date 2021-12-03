library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(mvnfast)
library(miscTools)
source("data_gen.R")
sourceCpp("Rcpp_functions.cpp")
library(matrixcalc)
library(clusterGeneration)
library(fBasics)
#==============================input======================================
k1 <- 30
k2 <- 10
nobs <- 100
theta_true <- 0.5
iteration <- 3000
burn <- 1000
lag <- 1
spec_rad<- 0.7
edge_density <- 0.04
horizon <- 20
N <- 1000
#============================================================================
vec <- permutation(k1,k2)[[1]]
vec1 <- permutation(k1,k2)[[2]]
vec2 <- permutation(k1,k2)[[3]]

p <- (3*k1) + k2
p1 <- (k1+k2)

A <- Gen_A(k1,k2,theta_true,lag, spec_rad, edge_density)[[1]]
sigma_true <- Gen_Sigma (k1,k2)[[1]]
desired_SNR <- 2
k <- norm(A, type="F")/norm(desired_SNR*sigma_true, type="F")
sigma_true <- k * sigma_true

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
  forecast_true_rep[,,r] <-forecast_true
    
  X <- data_t[1:nobs,]
  Y <- data_t[2:(nobs + 1),]
  LSE <- lm(Y ~ X - 1)$coef
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
  
  source("functions.R")
  
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
  
 
   omega <-list()
  sigma_adj <- numeric()
  
  for(i in 1:k1)
  {
    #sigma1 <- genPositiveDefMat(3, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,1), lambdaLow=1, ratioLambda=2)$Sigma
    sigma1 <- 0.3*diag(3)
    omega[[i]] <- solve(sigma1)
    sigma_adj[i] <- 1/(t(u)%*% solve(sigma1)%*% u)  
  }
  
  sigma2 <- rep(1,k2)
  alpha <- 1
  beta <- 1
  #Q <- genPositiveDefMat(3, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,1), lambdaLow=1, ratioLambda=2)$Sigma
  Q <- 0.3*diag(3)
  #df <- 5
  df <- 1
  
  
  #####################################TIME===============================
  
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
 
   
  sigma_inv_func <- sigma_inverse(omega,sigma2,k1,k2)
  sigma_11 <-as.matrix( sigma_inv_func[[1]])
  sigma_22 <- sigma_inv_func[[2]]
  sigma_inv<- sigma_inv_func[[3]]
  Final_Sigma_Inv <- sigma_inv_func[[4]]
  sigma_mat <- sigma_inv_func[[5]]
  
  A11_func <- p_A11(theta_mat,s,sigma_inv,k1,sigma_11,ss11,s12,w12)
  gamma_A11 <- A11_func$gamma_A11
  d_A11 <- A11_func$d
  q11 <- 0.7
  tao11 <- 0.5
  
  A12_func <- p_A12(k1,k2,u,sigma_inv,sigma_11,s12,ss12,s22,w11)
  gamma_A12 <- A12_func$gamma_A12
  d_A12 <- A12_func$d
  q12 <- 0.96
  tao12 <- 0.5
  
  A21_func <- p_A21( k1, k2, v, sigma_22, s, s11, ss21, s12,w22)
  gamma_A21 <- A21_func$gamma_A21
  d_A21 <- A21_func$d
  q21 <- 0.96
  tao21 <- 0.5
  
  A22_func <-  p_A22( v, sigma_22, s22,ss22, s12,w21)
  gamma_A22 <- A22_func$gamma_A22
  d_A22 <- A22_func$d
  q22 <- 0.7
  tao22 <- 0.5
  
  for (ind in 1:iteration)
  {
    b1 <- as.vector(gen_b1(k1,gamma_A11,b1,d_A11,sigma_adj,tao11,q11))
    b2 <- as.vector(gen_b2(k2,k1,gamma_A12,b2,d_A12,sigma_adj,tao12,q12))
    b3 <- as.vector(gen_b3(k1,k2,gamma_A21,b3,d_A21,sigma2,tao21,q21))
    b4 <- as.vector(gen_b4(k2,gamma_A22,b4,d_A22,sigma2,tao22,q22))
  
    theta_para <- theta_dist( b1, b2, b3, b4, k1, k2, Final_Sigma_Inv, s1, s2, s3)
    x_mat <- theta_para$x_mat
    theta_y <- theta_para$theta_y
    theta <- draw_theta_v1(N,d_theta,x_mat,theta_y)
    theta_func <- A_mat( theta, b1, b2, b3, b4, k1, k2)
    A11 <- theta_func$A11
    A12 <- theta_func$A12
    A21 <- theta_func$A21
    A22 <- theta_func$A22
    theta_mat <- theta_func$theta_mat
    u <- theta_func$u
    v <- theta_func$v
    w_func <- w_mat( A11,A12,A21,A22, theta_mat,u,v, k1,k2)
    w11 <- w_func$w11
    w12 <- w_func$w12
    w21 <- w_func$w21
    w22 <- w_func$w22
    w <- w_func$w
    temp <- w[vec1,]
    w_mod <- temp[,vec1]
    A_upper <- cbind(A11,A12)
    A_lower <- cbind(A21,A22)
    sigma_func <- parameter_sigma (k1,k2,nobs,u,Q,df,alpha,beta,data_mod,w_mod,A21,A22,tao11,tao21,tao22,A_upper,A_lower)
    omega <- sigma_func[[1]]
    sigma2 <- sigma_func[[2]]
    sigma_adj <- sigma_func[[3]]
    
    sigma_inv_func <- sigma_inverse(omega,sigma2,k1,k2)
    sigma_11 <- sigma_inv_func[[1]]
    sigma_22 <- sigma_inv_func[[2]]
    sigma_inv<- sigma_inv_func[[3]]
    Final_Sigma_Inv <- sigma_inv_func[[4]]
    sigma_mat <- sigma_inv_func[[5]]
    
    A11_func <- p_A11(theta_mat,s,sigma_inv,k1,sigma_11,ss11,s12,w12)
    gamma_A11 <- A11_func$gamma_A11
    d_A11 <- A11_func$d
    
    A12_func <- p_A12(k1,k2,u,sigma_inv,sigma_11,s12,ss12,s22,w11)
    gamma_A12 <- A12_func$gamma_A12
    d_A12 <- A12_func$d
    
    A21_func <- p_A21( k1, k2, v, sigma_22, s, s11, ss21, s12,w22)
    gamma_A21 <- A21_func$gamma_A21
    d_A21 <- A21_func$d
    
    A22_func <-  p_A22( v, sigma_22, s22,ss22, s12,w21)
    gamma_A22 <- A22_func$gamma_A22
    d_A22 <- A22_func$d
    
    if (ind > burn) 
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
  
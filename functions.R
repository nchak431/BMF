library(Matrix)
library(MCMCpack)
library(dplyr)

theta_matrix <- function(theta,b1,b2,b3,b4,k1,k2,vec1)
{
  theta_mat <- matrix(c(theta^2,theta^4,theta^3,1,(theta^2),theta,theta,theta^3,(theta^2)), nrow=3,ncol =3, byrow = T)
  u <- c(theta^2,1,theta)
  v <- c(1,theta^2,theta)
  A11 <- matrix(b1,nrow=k1,ncol=k1,byrow=F)
  A12 <- matrix(b2,nrow=k1,ncol=k2,byrow=F)
  A21 <- matrix(b3,nrow=k2,ncol=k1,byrow=F)
  A22 <- matrix(b4,nrow=k2,ncol=k2,byrow=F)
  w11 <- kronecker(theta_mat,A11)
  w12 <- kronecker(u,A12)
  w21 <- kronecker(t(v),A21)
  w22 <- A22
  w1 <- cbind(w11, w12)
  w2 <- cbind(w21, w22)
  w<- rbind(w1,w2)
  temp <- w[vec1,]
  w_mod <- temp[,vec1]
  data <- list(theta_mat,u,v,w11,w12,w21,w22,w,w_mod)
  data
}


sigma_inverse <- function(omega,sigma2,k1,k2)
{
  omega1 <- bdiag(omega)
  omega2 <- diag((1/sigma2))
  sigma_Inv <- as.matrix(bdiag(omega1,omega2))
  sigma_mat <- as.matrix(solve(sigma_Inv))
  temp <- sigma_Inv[,vec2]
  Final_Sigma_Inv <- temp[vec2,]
  sigma_11 <- Final_Sigma_Inv[1:(3*k1),1:(3*k1)]
  sigma_22 <- Final_Sigma_Inv[((3*k1)+1):((3*k1)+k2), ((3*k1)+1):((3*k1)+k2)]
  
  sigma_inv <- NULL
  for ( i in 1:3)
  {
    sigma_inv[[i]] <- list()
    for ( j in 1: 3)
    {
      sigma_inv[[i]][[j]] <- sigma_11[(((i-1)*k1)+1):(i*k1),(((j-1)*k1)+1):(j*k1)]
    }
  }
  data <- list(sigma_11,sigma_22,sigma_inv,Final_Sigma_Inv,sigma_mat)
  data
}

S_matrix <- function(data,n,k1,k2)
{
  s1 <- data[,1:n] %*% t(data[,1:n])
  s2 <- data[,2:(n+1)] %*% t(data[,1:n])
  s3 <- data[,2:(n+1)] %*% t(data[,2:(n+1)])
  
  s11 <- s1[1:(3*k1), 1:(3*k1)]
  s12 <- s1[1:(3*k1), ((3*k1)+1):((3*k1)+k2)]
  s21 <- s1[((3*k1)+1):((3*k1)+k2), 1:(3*k1)]
  s22 <- s1[((3*k1)+1):((3*k1)+k2), ((3*k1)+1):((3*k1)+k2)]
  
  ss11 <- s2[1:(3*k1), 1:(3*k1)]
  ss12 <- s2[1:(3*k1), ((3*k1)+1):((3*k1)+k2)]
  ss21 <- s2[((3*k1)+1):((3*k1)+k2), 1:(3*k1)]
  ss22 <- s2[((3*k1)+1):((3*k1)+k2), ((3*k1)+1):((3*k1)+k2)]
  
  s <- NULL
  for ( i in 1: 3)
  {
    s[[i]] <- list()
    
    for (j in 1: 3)
    {
      s[[i]][[j]] <- s11[(((i-1)*k1)+1):(i*k1),(((j-1)*k1)+1):(j*k1)]
    }
  }
  data <- list(s1,s2,s3,s11,s12,s21,s22,ss11,ss12,ss21,ss22,s)
  data
}


#---------------------Parameters for A11-------------------------

para_A11 <- function(theta_mat,s,sigma_inv,sigma_11,s11,ss11,s12,w12)
{
  gamma_A11 <- matrix(0,nrow=(k1*k1),ncol=(k1*k1))  
  for (i in 1:3)
  {
    for(j in 1:3)
    {
      mat1 <- matrix(0,nrow=k1,ncol=k1)
      mat2 <- matrix(0,nrow=k1,ncol=k1)
      for(k in 1:3)
      {
        mat1<- mat1+ s[[i]][[k]] * (t(theta_mat))[k,j]
        mat2 <- mat2 + sigma_inv[[j]][[k]] * theta_mat[k,i]
        
      }
      gamma1 <- kronecker(t(mat1),mat2)
      gamma_A11 <- gamma_A11 + gamma1  
    }
  }
  
  P_mat <- sigma_11 %*% (ss11  - (w12 %*% t(s12)))
  
  P <- NULL
  for ( i in 1:3)
  {
    P[[i]] <- list()
    for ( j in 1: 3)
    {
      P[[i]][[j]] <- P_mat[(((i-1)*k1)+1):(i*k1),(((j-1)*k1)+1):(j*k1)]
    }
  }
  
  D <- matrix(0,nrow=k1,ncol=k1)  
  for (i in 1:3)
  {
    for(j in 1:3)
    {
      D <- D + (t(theta_mat))[i,j] * P[[j]][[i]]
      
    }
  }
  d_A11 <- vec(as.matrix(D))
  data <- list(gamma_A11,d_A11)
  data
}

#---------------------Parameters for A12-------------------------

para_A12 <- function(k1,k2,u,sigma_inv,sigma_11,s12,ss12,s22,w11)
{
  gamma <- matrix(0,nrow=k1,ncol=k1)  
  
  for(i in 1:3)
  {
    
    for(j in 1:3)
    {
      gamma <- gamma + ((u[i] * u[j]) * sigma_inv[[i]][[j]])
      
    }
  }
  gamma_A12 <- kronecker(t(s22),gamma)
  
  P_mat <- sigma_11 %*%(ss12 -(w11 %*% s12 ))
  P <- list()
  for(i in 1:3)
  {
    P[[i]] <- P_mat[(((i-1)*k1)+1):(i*k1),1:k2]
  }
  
  sum <- 0
  
  for(i in 1:3)
  {
    sum <- sum + u[i]*P[[i]] 
    
  }
  d_A12 <- vec(as.matrix(sum))
  data <- list(gamma_A12,d_A12)
  data
}

#---------------------Parameters for A21-------------------------

para_A21 <- function(k1,k2,v,sigma_22,s,s11,ss21,s12,w22)
{
  gamma_A21 <- matrix(0,nrow=(k1*k2),ncol=(k1*k2))  
  
  
  for(i in 1:3)
  {
    for(j in 1:3)
    {
      gamma_A21 <- gamma_A21 + kronecker((v[i] * t(s[[i]][[j]])), (sigma_22 * v[j]))
      
    }
  }
  
  P_mat <- sigma_22 %*% (ss21 -(w22 %*% t(s12) ))
  P <- list()
  for(i in 1:3)
  {
    P[[i]] <- P_mat[1:k2,(((i-1)*k1)+1):(i*k1)]
  }
  
  sum <- 0
  
  for(i in 1:3)
  {
    sum <- sum + (v[i]*P[[i]] )
    
  }
  d_A21 <- vec(as.matrix(sum))
  data <- list(gamma_A21,d_A21)
  data
}

#---------------------Parameters for A22-------------------------

para_A22 <- function(v,sigma_22,s22,ss22,s12,w21)
{
  
  mean <- (ss22 %*% solve(s22))-(w21 %*% s12 %*% solve(s22))
  m <- vec(mean)
  gamma_A22 <- solve(kronecker(solve(s22),solve(sigma_22)))
  d_A22 <- gamma_A22 %*% m
  data <- list(gamma_A22,d_A22)
  data
}


#------------------Drawing entries of A22 from mixture dist-------------------

draw <- function(x,y,p,n)
{
  a <- rbinom(1,1,p)
  if(a==1)
  {
    draw <- 0 
  }else{
    draw <- rnorm(1,x,y)  
  }
  draw
  
}



#-------------estimating coeff of pos dist of theta-----------------
theta_dist_coeff <- function(s1,s2,s3,Final_Sigma_Inv,b1,b2,b3,b4,k1,k2,vec1)
{
  theta_vec <- seq(0.1,1,length.out=20)
  state <- c(-4,-3,-2,-1,1,2,3,4)
  W <- list()
  x <- matrix(0,nrow=length(theta_vec),ncol=length(state))
  q1 <- numeric()
  d <- numeric()
  
  for(i in 1:length(theta_vec))
  {
    W[[i]]<- theta_matrix(theta_vec[i],b1,b2,b3,b4,k1,k2,vec1)[[8]] 
    q1[i] <- tr(as.matrix((Final_Sigma_Inv %*% s3) - (2*(t(W[[i]]) %*% Final_Sigma_Inv %*% s2)) + (t(W[[i]]) %*% Final_Sigma_Inv %*% W[[i]] %*% s1)  ))
    
    for (j in 1: length(state))
    {
      x[i,j]=(theta_vec[i])^state[j]
      
    }
  }
  
  d_theta <- as.vector(lm(q1 ~ x)$coef )
  d_theta <- d_theta[-1]
  d_theta
  
}

#---------------drawing theta--------------------

draw_theta <- function(N,d_theta)
{
  seq <- c(-4,-3,-2,-1,1,2,3,4)
  lnp <- numeric()
  x1 <- seq((1/(2*N)),((2*(N-1))/(2*N)),by=(1/N))
  for (i in 1:length(x1))
  {
    q1 <- (x1[i] ^ seq)
    lnp[i] <- sum(d_theta*q1)*(-1/2)
    
  }
  if (max(lnp) < -750 | max(lnp) > 700) lnp <- lnp + (700 - max(lnp))
  theta <- sample(x1,1,prob=exp(lnp))
  theta
  
}

draw_theta_v1 <- function(N,d_theta,x_mat,theta_y)
{
  d_theta <- as.vector(lm(theta_y ~ x_mat)$coef )
  d_theta <- d_theta[-1]
  seq <- c(-4,-3,-2,-1,1,2,3,4)
  lnp <- numeric()
  x1 <- seq((1/(2*N)),((2*(N-1))/(2*N)),by=(1/N))
  for (i in 1:length(x1))
  {
    q1 <- (x1[i] ^ seq)
    lnp[i] <- sum(d_theta*q1)*(-1/2)
    
  }
  if (max(lnp) < -750 | max(lnp) > 700) lnp <- lnp + (700 - max(lnp))
  theta <- sample(x1,1,prob=exp(lnp))
  theta
  
}

#-----------------parameters for Sigma----------------

parameter_sigma <- function(k1,k2,n,u,Q,df,alpha,beta,data_mod,w_mod,A21_mat,A22_mat,tao11,tao21,tao22,A_upper,A_lower)
{
  s1_mod <- data_mod[,1:n] %*% t(data_mod[,1:n])
  s2_mod <- data_mod[,2:(n + 1)] %*% t(data_mod[,1:n])
  s3_mod <- data_mod[,2:(n + 1)] %*% t(data_mod[,2:(n + 1)])
  
  mat <- s3_mod - (s2_mod%*%(t(w_mod))) - t(s2_mod%*%(t(w_mod))) +(w_mod %*% s1_mod %*% t(w_mod))
  mat_1st_block <- mat[1:(3*k1), 1:(3*k1)]
  e1 <- c(1,0,0)
  e2 <- c(0,1,0)
  c <- matrix(t(c(e1,e2,u)),nrow=3,ncol=3,byrow = F)
  
  v <- list()
  omega <- list()
  sigma_adj <- numeric()
  sigma2 <- numeric()
  
  for(i in 1:k1)
  {
    t <- mat_1st_block[(((i-1)*3)+1):(i*3),(((i-1)*3)+1):(i*3)]
    t1 <-  solve(c) %*% t %*% solve(t(c))
    d <- t1 + Q
    A <- rwish((n + df - 1),solve(d[1:2,1:2]))
    #m <- (-((as.matrix(d[1:2,3]) + as.matrix(d[3,1:2]))/2))
    m <- -(as.matrix(d[1:2,3]))
    sum <- sum(A_upper[i,]^2) 
    count <- (k1+k2)-length(which(A_upper[i,]==0))
    alpha1 <- ((n+df)/2) + (count/2)
    beta1 <- tr((sum/(tao11*tao11)) + d[3,3] - (t(m) %*% solve(d[1:2,1:2]) %*% m))/2
    b <- rgamma(1, shape=alpha1, rate = beta1)
    a <-  mvrnorm(1,((solve(as.matrix(d[1:2,1:2]))) %*% m),solve((b*d[1:2,1:2])))
    v1 <- matrix(c((A+ (a%*%t(a))*b),(a*b)),nrow=2,byrow=F)
    v2 <- matrix(c(t(a*b),b),nrow=1,byrow=F)
    v[[i]] <- rbind(v1,v2)
    omega[[i]] <- solve(t(c)) %*% v[[i]] %*% solve(c)
    #sigma1[[i]] <- solve(omega)
    sigma_adj[i] <- 1/(t(u)%*% omega[[i]] %*% u)
    
  }
  
  for(i in 1:k2)
  {
    p <- mat[((3*k1)+i),((3*k1)+i)]
    sum1 <- sum(A21_mat[i,]^2) 
    sum2 <- sum(A22_mat[i,]^2)
    count <- (k1+k2)-length(which(A_lower[i,]==0))
    q1 <- (alpha + (n/2)+(count/2))
    q2 <- (beta + (p/2) + (sum1/(2*tao21*tao21)) + (sum2/(2*tao22*tao22)))
    sigma2[i] <- rinvgamma(1,q1,q2)
    
  }
  
  data <- list(omega,sigma2,sigma_adj)  
  data  
  
}

forecast_cred <- function(sim,true,horizon)
{
  sim_q11 <- apply(sim,c(1,2),quantile,probs=0.025)
  sim_q1 <- apply(sim,c(1,2),quantile,probs=0.25)
  sim_med <- apply(sim,c(1,2),quantile,probs=0.50)
  sim_q3 <- apply(sim,c(1,2),quantile,probs=0.75)
  sim_q33 <- apply(sim,c(1,2),quantile,probs=0.975)
  q11 <- vec(sim_q11);q1 <- vec(sim_q1);med <- vec(sim_med);q3 <- vec(sim_q3);q33 <- vec(sim_q33)
  err <- numeric()
  for ( i in 1:horizon)
  {
    err[i] <- (sum((true[,i]-sim_med[,i])^2))/length(true[,i])
  }
  err_low <- numeric()
  for ( j in 1:horizon)
  {
    err_low[j] <- sum((true[,j][((3*k1)+1):((3*k1)+k2)]-sim_med[,j][((3*k1)+1):((3*k1)+k2)])^2)/k2
  }  
  return(list(forecast_med = sim_med,ferr_med=err,ferr_med_low=err_low ))
}

metric<- function(true,est)
{
  a <- true*est
  TP <- length(a) - (length(which(a==0)))
  true_nonzero <- (length(true) - (length(which(true==0))))
  FN <- true_nonzero - TP
  b <- true + est
  FP <- (length(b) - (length(which(b==0))))- true_nonzero
  TN <- (length(which(true==0))) - FP
  d <- (TP* TN) - (FP * FN)
  c1 <- (TP + FP)
  c2 <- (TP + FN)
  c3 <- (TN + FP)
  c4 <- (TN + FN)
  SN <- TP/c2
  SP <- TN/c3
  data <- list(SN,SP)
  data
}


metric_A <- function(A_block_one,A_est)
{
  m <- hadamard.prod(A_block_one,A_est)
  m1 <- colSums(m!=0)
  TP_A_mat <- sum(m1)
  m2 <- colSums(A_block_one!=0)
  true_nonzero <- sum(m2)
  FN_A_mat <- true_nonzero - TP_A_mat
  n <- A_block_one + A_est
  n1 <- colSums(n!=0)
  FP_A_mat <- sum(n1)- true_nonzero
  true_zero <- (p1^2) - true_nonzero
  TN_A_mat <- true_zero - FP_A_mat
  d <- (TP_A_mat* TN_A_mat) - (FP_A_mat * FN_A_mat)
  c1 <- (TP_A_mat + FP_A_mat)
  c2 <- (TP_A_mat + FN_A_mat)
  c3 <- (TN_A_mat + FP_A_mat)
  c4 <- (TN_A_mat + FN_A_mat)
  SN <- TP_A_mat/c2
  SP <- TN_A_mat/c3
  data <- list(SN,SP)
  data
}

# For implemenetation of `MFBVAR'
mfb_rearr <- function(k1,k2,horizon)
{
  arr <- NULL   
  for(i in 1:horizon)
  {
    m1 <- NULL;m2 <- NULL;m3 <- NULL;m4 <- NULL
    for(j in 1:k1)
    {
      m <- ((j-1)*(3*horizon))+(3*i)
      m1 <- c(m1,m)  
      m <- ((j-1)*(3*horizon))+((3*i)-2)
      m2 <- c(m2,m)
      m <- ((j-1)*(3*horizon))+((3*i)-1)
      m3 <- c(m3,m)
    }  
    for(j in 1:k2)
    {
      m <- (3*horizon*k1) + ((j-1)*(3*horizon))+(3*i)
      m4 <- c(m4,m)  
      
    }  
     
    arr <- c(arr,m1,m2,m3,m4) 
    arr
  }
  return(arr=arr)
}



mfb <- function(data_mod,k1,k2,horizon,arr,true)
{
  mf_list <- list() 
  for(i in 1:(k1+k2))  
  {
    if(i <= k1) {
      index <- c((3*i)-2,(3*i)-1,(3*i))  
      sub <- data_mod[index,]
      mf_list[[i]] <- vec(sub)  
    }else {
      sub <- matrix(data_mod[((2*k1)+i),],ncol=1)
      mf_list[[i]] <- vec(sub)      
    } 
  }
  alfred_to_ts <- function(x,freq){
    ts(x[,1],frequency=freq)
  }
  mf <- mapply(alfred_to_ts,x=mf_list,freq=c(rep(12,k1),rep(4,k2)))
  quantile1 <- numeric()
  quantile2 <- numeric()
  for(i in 1:(k1+k2))
  {
    quantile1[i] <- quantile(mf[[i]], probs = 0.025)
    quantile2[i] <- quantile(mf[[i]], probs = 0.975)
  }
  int <- cbind(quantile1,quantile2)
  moments <- interval_to_moments(int)
  prior <- set_prior(Y = mf, n_lags = 4, n_reps = 1000)
  prior 
  prior1 <-  update_prior(prior, n_fcst = (3*horizon))
  est1 <- estimate_mfbvar(prior1, prior = "minn", variance = "iw")
  
  prior2 <- update_prior(prior,d = "intercept",prior_psi_mean = moments$prior_psi_mean,
                         prior_psi_Omega = moments$prior_psi_Omega)
  prior2 <-  update_prior(prior2, n_fcst = (3*horizon))
  est2 <- estimate_mfbvar(prior2, prior = "ss", variance = "iw")
  est3 <- estimate_mfbvar(prior2, prior = "ssng", variance = "iw")
  pred_1 <- predict(est1,aggregate_fcst = FALSE, pred_bands = 0.9)
  pred_2 <- predict(est2,aggregate_fcst = FALSE, pred_bands = 0.9)
  pred_3 <- predict(est3,aggregate_fcst = FALSE, pred_bands = 0.9)
  
  #===========only quarterly for nowcast=======================
  
  arr1 <- NULL   
  for(j in 1:(3*horizon))
  {
    m1 <- NULL
    for(i in 1:k2)
    {
      m <- j+(3*horizon*(i-1))
      m1 <- c(m1,m)  
      
    }
    arr1 <- c(arr1,m1) 
    arr1
  }
  ############# Minn ###############
pred1=pred_1[-(1:(3*horizon*k1)),]
pred1 <- as.matrix(pred1[,4:6])
pred1 <- pred1[arr1,]
# mean = median(under normality) = pred1[,2]
sd_minn <- (pred1[,3]-pred1[,2])/2
mfb_mean_minn <- matrix(pred1[,2],nrow=k2,ncol=(3*horizon),byrow=F)
mfb_sd_minn <- matrix(sd_minn,nrow=k2,ncol=(3*horizon),byrow=F)

 err_low_1 <- numeric()
  for ( j in 1:horizon)
  {
    set <- c((3*j-2),(3*j-1),(3*j))
    for(i in set){
      err_low_1[i] <- sum((true[,j][((3*k1)+1):((3*k1)+k2)]-mfb_mean_minn[,i])^2)/k2
    }
  }
 

crps_mat_minn <- matrix(0,nrow = k2,ncol=(3*horizon))
logs_mat_minn <- matrix(0,nrow = k2,ncol=(3*horizon))

for ( j in 1:horizon)
{
  set <- c((3*j-2),(3*j-1),(3*j))
  for(i in set){
    crps_mat_minn[,i] <- crps(y = true[,j][((3*k1)+1):((3*k1)+k2)], family = "normal", mean = mfb_mean_minn[,i], sd = mfb_sd_minn[,i])
    logs_mat_minn[,i] <- logs(y = true[,j][((3*k1)+1):((3*k1)+k2)], family = "normal", mean = mfb_mean_minn[,i], sd = mfb_sd_minn[,i])
    
  }
}
  
crps_mat_minn
logs_mat_minn

###################### SS #################################
pred2 <- pred_2[-(1:(3*horizon*k1)),]
pred2 <- as.matrix(pred2[,4:6])
pred2 <- pred2[arr1,]
# mean = median(under normality) = pred2[,2]
sd_ss <- (pred2[,3]-pred2[,2])/2
mfb_mean_ss <- matrix(pred2[,2],nrow=k2,ncol=(3*horizon),byrow=F)
mfb_sd_ss <- matrix(sd_ss,nrow=k2,ncol=(3*horizon),byrow=F)

 err_low_2 <- numeric()
  for ( j in 1:horizon)
  {
    set <- c((3*j-2),(3*j-1),(3*j))
    for(i in set){
      err_low_2[i] <- sum((true[,j][((3*k1)+1):((3*k1)+k2)]-mfb_mean_ss[,i])^2)/k2
    }
  }
 


crps_mat_ss <- matrix(0,nrow = k2,ncol=(3*horizon))
logs_mat_ss <- matrix(0,nrow = k2,ncol=(3*horizon))

for ( j in 1:horizon)
{
  set <- c((3*j-2),(3*j-1),(3*j))
  for(i in set){
    crps_mat_ss[,i] <- crps(y = true[,j][((3*k1)+1):((3*k1)+k2)], family = "normal", mean = mfb_mean_ss[,i], sd = mfb_sd_ss[,i])
    logs_mat_ss[,i] <- logs(y = true[,j][((3*k1)+1):((3*k1)+k2)], family = "normal", mean = mfb_mean_ss[,i], sd = mfb_sd_ss[,i])
    
  }
}

############################# SSNG ################################

pred3 <- pred_3[-(1:(3*horizon*k1)),]
pred3 <- as.matrix(pred3[,4:6])
pred3 <- pred3[arr1,]
# mean = median(under normality) = pred2[,2]
sd_ssng <- (pred3[,3]-pred3[,2])/2
mfb_mean_ssng <- matrix(pred3[,2],nrow=k2,ncol=(3*horizon),byrow=F)
mfb_sd_ssng <- matrix(sd_ssng,nrow=k2,ncol=(3*horizon),byrow=F)

 err_low_3 <- numeric()
  for ( j in 1:horizon)
  {
    set <- c((3*j-2),(3*j-1),(3*j))
    for(i in set){
      err_low_3[i] <- sum((true[,j][((3*k1)+1):((3*k1)+k2)]-mfb_mean_ssng[,i])^2)/k2
    }
  }
 

crps_mat_ssng <- matrix(0,nrow = k2,ncol=(3*horizon))
logs_mat_ssng <- matrix(0,nrow = k2,ncol=(3*horizon))

for ( j in 1:horizon)
{
  set <- c((3*j-2),(3*j-1),(3*j))
  for(i in set){
    crps_mat_ssng[,i] <- crps(y = true[,j][((3*k1)+1):((3*k1)+k2)], family = "normal", mean = mfb_mean_ssng[,i], sd = mfb_sd_ssng[,i])
    logs_mat_ssng[,i] <- logs(y = true[,j][((3*k1)+1):((3*k1)+k2)], family = "normal", mean = mfb_mean_ssng[,i], sd = mfb_sd_ssng[,i])
    
  }
}
  return(list(pred_1=pred_1[,4:6],pred_2=pred_2[,4:6],pred_3=pred_3[,4:6],crps_minn=crps_mat_minn,logs_minn=logs_mat_minn,
              crps_ss=crps_mat_ss,logs_ss=logs_mat_ss,crps_ssng=crps_mat_ssng,logs_ssng=logs_mat_ssng,
              ferr_mfb_low1=err_low_1,ferr_mfb_low2=err_low_2,ferr_mfb_low3=err_low_3))
}

#####################
### GENERATES VAR TRANSITION MATRIX for BMF
#####################
library(MASS)
library(Matrix)

gensparse_A = function(  # returns a list of d matrices A_1, ..., A_d
  p    # generate p x p VAR
  ,k1 #high freq variables
  ,k2 #low freq variables
  ,d=1   # generate VAR(d)
  ,max_eig # = 0.6   spectral norm of A
  ,edge_density # = 0.1  if diferent for different lags, provide a vector of length d
  ,nonzero_diag # = 1  ensure all diagonal entries are non-zero; if different for different lags, provide a d-vector of 0 and 1
  ,stationary = 1  # ensure stationarity
  ,network.family = "random"
  ,theta
  # A matrix filled randomly or with specific structure
){
  #set.seed = 1729

    A = list()
    for (lag in 1:d){
      e_d = ceiling(p^2 * ifelse(length(edge_density) == 1, edge_density, edge_density[lag]))
      temp = sample(p^2, e_d)
      temp_v = rep(0, p^2)
      temp_v[temp] = runif(e_d, max_eig, 10*max_eig)*ifelse(rbinom(e_d,1,0.7),-1,1)
      A[[lag]] = array(temp_v/lag, c(p,p))
      if (ifelse(length(nonzero_diag) == 1, nonzero_diag, nonzero_diag[lag]))
        diag(A[[lag]]) = runif(p, max_eig, 10)
    }
    
    if (stationary){
      A.big = array(0, c(p*d, p*d))
      
      for (i in 1:d)
        A.big[1:p, ((i-1)*p+1):(i*p)] = max_eig*A[[i]]
      
      if (d > 1)
        diag(A.big[(p+1):(p*d), 1:((d-1)*p)]) = 1
      hf <- list()
      for (i in 1:d)
      {
        hf[[i]] <- A.big[(((i-1)*k1)+1):(i*k1), (((i-1)*k1)+1):(i*k1)]
        hf[[i]] <- matrix(runif((k1^2),(abs(max(hf[[i]]))-2),abs(max(hf[[i]])))*ifelse(rbinom(k1^2,1,0.7),-1,1),nrow=k1,ncol=k1) #no sparsity in high-freq block
        A.big[(((i-1)*k1)+1):(i*k1), (((i-1)*k1)+1):(i*k1)] <- hf[[i]]
      }             
      for (i in 1:d)
      {
        temp = max(abs(eigen(hf[[i]])$values))  
        count = 0
        while (temp > 0.31){
          count = count+1
          hf[[i]] = hf[[i]]*0.95
          temp = max(abs(eigen(hf[[i]])$values)) 
          A.big[(((i-1)*k1)+1):(i*k1), (((i-1)*k1)+1):(i*k1)] <- hf[[i]]
        }  
      }
      eig = max(abs(eigen(A.big)$values))
      target= 0.9 - (2*(theta^2)*0.3)
      count = 0
      while (eig>target){
        count = count+1
        A.big[(k1+1):p,] = A.big[(k1+1):p,]*0.95
        A.big[1:k1,(k1+1):p] = A.big[1:k1,(k1+1):p]*0.95
        eig = max(abs(eigen(A.big)$values))
      }
    }
    for (i in 1:d)
      A[[i]] = A.big[1:p, ((i-1)*p+1):(i*p)]
     return(list(A=A,Eigen=temp,
              Signal=round(max(abs(A.big[1:p,])), 2))
  )
}


VAR_Mat <- function (p,k1,k2, lag, spec_rad, edge_density,theta)
{ 
  a <- gensparse_A(p, k1,k2,lag, spec_rad, edge_density, 0, 1, "random",theta)$A
  Mat <- do.call(cbind, a)
  return (Mat)
}

Gen_A = function(k1,k2,theta,lag, spec_rad, edge_density)
{
  p <- k1+k2
  mu <- c(1/(theta* theta), 1/theta) 
  gamma <- c((theta* theta), theta)
  A_block_one <- VAR_Mat (p,k1,k2,lag, spec_rad, edge_density,theta)
  temp <- vec(A_block_one[1:k1,(k1+1):(k1+k2)])
  temp[which(abs(temp) > 0 & abs(temp) < 0.5)]=0
  A_block_one[1:k1,(k1+1):(k1+k2)] = matrix(temp,nrow=k1,byrow=F)
  A_block_one[1:k1,1:(k1+k2)] <- A_block_one[1:k1,1:(k1+k2)]*(theta ^ 2)
  
  temp <- vec(A_block_one[(k1+1):(k1+k2),1:k1])
  temp[which(abs(temp) > 0 & abs(temp) < 0.5)]=0
  A_block_one[(k1+1):(k1+k2),1:k1] = matrix(temp,nrow=k2,byrow=F)
  
  A_block_two <-kronecker(A_block_one[1:k1,], mu)
  A_block_three <-kronecker(A_block_one[,1:k1], t(gamma))
  v <- kronecker(t(gamma), mu)
  A_block_four <- kronecker(A_block_one[1:k1,1:k1], v)
  A1 <- rbind(A_block_one, A_block_two)
  A2 <- rbind(A_block_three, A_block_four)
  A <- cbind(A1,A2)
  data <- list(A,A_block_one)
  data
  }


#=======================================Generate Sigma Matrix=================================
Gen_Sigma = function(k1,k2)
{
  set.seed(1233)
  c <- list()
  d <- list()
  p <- k1 + k2
  sigma_sq <- runif(p, 0.1, 1)
  ro <- rep(0.1,k1)
  sig_block1 <- diag(sigma_sq)
  for(i in 1:k1)
  {
    c[[i]] <- matrix(c(((ro[i]^2)*sigma_sq[i]),(ro[i]*sigma_sq[i])),nrow=2,ncol=1)
  }
  sig_block2 <- bdiag(c)
  null <- matrix(0,(2*k1),k2)
  sig_block2 <- cbind(sig_block2, null)
  sig_block3 <- t(sig_block2)
  for(i in 1:k1)
  {
    d[[i]] <- sigma_sq[i] * matrix(c(1,ro[i],ro[i],1), nrow=2)
  }
  sig_block4 <- bdiag(d)
  Sig1 <- rbind(sig_block1, sig_block2)
  Sig2 <- rbind(sig_block3, sig_block4)
  sigma_true <- cbind(Sig1,Sig2)
  #sigma_true
  data <- list(sigma_true,sigma_sq)
  data
}


#------permutation vectors----------

permutation <- function(k1,k2)
{
  y0<- NULL
  for(i1 in 1:k1)
  {
    x <- i1
    y0 <- c(y0,x)
  }
  y0
  
  y1<- NULL
  for(i1 in 1:k1)
  {
    x <- k1+k2+(2*i1-1)
    y1 <- c(y1,x)
  }
  y1
  
  y2<- NULL
  for(i1 in 1:k1)
  {
    x <- k1+k2+(2*i1)
    y2 <- c(y2,x)
  }
  y2
  
  z<- NULL
  for(i1 in 1:k2)
  {
    x <- k1+i1
    z <- c(z,x)
  }
  z
  vec <- c(y0,y1,y2,z)
  
  m1<- NULL
  for(i1 in 1:k1)
  {
    m <- c(k1+i1,2*k1+i1,i1)
    m1 <- c(m1,m)
  }
  m1
  
  m2<- NULL
  for(i1 in 1:k2)
  {
    m <- 3*k1+i1
    m2 <- c(m2,m)
  }
  
  vec1 <- c(m1,m2)
  
  a<- NULL
  for(i1 in 1:k1)
  {
    m <- 3*i1
    a <- c(a,m)
  }
  a
  
  b<- NULL
  for(i1 in 1:k1)
  {
    m <- 3*i1 - 2
    b <- c(b,m)
  }
  b
  
  d<- NULL
  for(i1 in 1:k1)
  {
    m <- 3*i1 - 1
    d <- c(d,m)
  }
  d
  
  m2<- NULL
  for(i1 in 1:k2)
  {
    m <- 3*k1+i1
    m2 <- c(m2,m)
  }
  m2
  vec2 <- c(a,b,d,m2)
  data <- list(vec,vec1,vec2)
  data
}


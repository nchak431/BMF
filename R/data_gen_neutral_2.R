library(MASS)
library(Matrix)

Gen_A = function(k1,k2)
{
val <- runif((3*(k1+k2))^2,0.7,2)*ifelse(rbinom((3*(k1+k2))^2,1,0.7),-1,1)
mat <- matrix(val,nrow=(3*(k1+k2)),ncol=(3*(k1+k2)))
eig=max(abs(eigen(mat)$values))
while (eig > 0.81){
  mat = mat*0.95
  eig = max(abs(eigen(mat)$values)) 
}  
mat
}

# permutation
permutation <- function(k1,k2)
{


b<- NULL
for(i1 in 1:(3*k1))
{
  m <- i1
  b <- c(b,m)
}
b

a<- NULL
for(i1 in 1:k2)
{
  m <- ((3*k1)+(2*k2)+i1)
  a <- c(a,m)
}
a
vech <- c(b,a)

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
data <- list(vech,vec1,vec2)
data
}

Gen_Sigma = function(k1,k2)
{
  p <- 3*(k1 +k2)
  A <- genPositiveDefMat(p, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,1), lambdaLow=1, ratioLambda=2)$Sigma
  sigma_true <- A %*% (t(A))
  sigma_true
}

#==========================================
k1 <- 3
k2 <- 30
nobs <- 100
theta_true <- 0.5
horizon <- 10
#==========================================
vech <- permutation(k1,k2)[[1]]  
vec1 <- permutation(k1,k2)[[2]]

p <- (3*k1) + k2
p1 <- (k1+k2)

A <- Gen_A(k1,k2)
sigma_true <- Gen_Sigma (k1,k2)
desired_SNR <- 2
k <- norm(A, type="F")/norm(desired_SNR*sigma_true, type="F")
sigma_true <- k * sigma_true

#=====================Generate data=====================
y_true <- list()
y <- list()
y_mod <- list()
err <- list()
y_true[[1]] <- mvrnorm(1,rep(0,(3*(k1+k2))),sigma_true)
y[[1]]<- y_true[[1]][vech]
y_mod[[1]]<- y[[1]][vec1]
for ( t in 2:(nobs + horizon + 1))
{
  err[[t]] <- mvrnorm(1,rep(0,(3*(k1+k2))),sigma_true)
  y_true[[t]] <- as.vector((A %*% y_true[[t-1]]) + err[[t]])
  y[[t]]<-y_true[[t]][vech]
  y_mod[[t]]<-y[[t]][vec1]
  
}
sigma_true <- sigma_true[vech,]
sigma_true <- sigma_true[,vech]

y_full <- as.matrix(do.call(cbind,y))
y_mod_full <- do.call(cbind,y_mod)
data <-y_full[,1:(nobs + 1)]

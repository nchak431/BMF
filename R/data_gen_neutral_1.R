#========Generate VAR transition matrix using a neutral setting=======
Gen_A = function(k1,k2,rho_h,rho_l)
{
hh <- matrix(0,nrow=(3*k1),ncol=(2*k1))
h1 <- rho_h*diag(k1)
h2 <- (rho_h^2)*diag(k1)
h3 <- (rho_h^3)*diag(k1)
hh_block <- cbind(hh,rbind(h1,h2,h3))
ll_block <- rho_l*diag(k2)
temp1 <- runif((k1*k2),0.7,2)
temp2 <- runif((k1*k2),0.7,2)
hl_main <- matrix(temp1,nrow=k1,ncol=k2)
hl1 <- (1 + rho_h)* hl_main
hl2 <- (1+ rho_h + (rho_h^2))* hl_main
lh_main <- matrix(temp2,nrow=k2,ncol=k1)
hl_block <- rbind(hl_main,hl1,hl2)
lh1 <- lh_main*0.3
lh2 <- lh_main*0.4
lh3 <- lh_main*0.3
lh_block <- cbind(lh_main,lh1,lh2)
mat <-rbind(cbind(hh_block,hl_block),cbind(lh_block,ll_block))
eig=max(abs(eigen(mat)$values))
while (eig > 0.81){
  mat = mat*0.95
  eig = max(abs(eigen(mat)$values)) 
}  
mat
}

permutation <- function(k1,k2)
{
a<- NULL
for(i1 in 1:k1)
{
  m <- (2*k1+i1)
  a <- c(a,m)
}
a

b<- NULL
for(i1 in 1:k1)
{
  m <- i1
  b <- c(b,m)
}
b

d<- NULL
for(i1 in 1:k1)
{
  m <- (k1+i1)
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
vech <- c(a,b,d,m2)

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

Gen_Sigma = function(k1,k2,rho_h)
{
  A <- genPositiveDefMat((3*k1+k2), covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),  alphad=1, eta=1, rangeVar=c(0.1,0.9), lambdaLow=1, ratioLambda=2)$Sigma
  sigma <- A %*% (t(A))
  sigma
}

#===========================================
k1 <- 3
k2 <- 30
nobs <- 100
theta_true <- 0.5
horizon <- 20
rho_h <- 0.5
rho_l <- 0.6
#===========================================
vech <- permutation(k1,k2)[[1]]
vec1 <- permutation(k1,k2)[[2]]

p <- (3*k1) + k2
p1 <- (k1+k2)

A <- Gen_A(k1,k2,rho_h,rho_l)
sigma_true <- Gen_Sigma (k1,k2,rho_h)
desired_SNR <- 2
k <- norm(A, type="F")/norm(desired_SNR*sigma_true, type="F")
sigma_true <- k * sigma_true
#=========================Generate data=============================
y_true <- list()
y <- list()
y_mod <- list()
err <- list()
y_true[[1]] <- mvrnorm(1,rep(0,p),sigma_true)
y[[1]]<-y_true[[1]][vech]
y_mod[[1]]<- y[[1]][vec1]
for ( t in 2:(nobs + horizon + 1))
{
  err[[t]] <- mvrnorm(1,rep(0,p),sigma_true)
  y_true[[t]] <- as.vector((A %*% y_true[[t-1]]) + err[[t]])
  y[[t]]<-y_true[[t]][vech]
  y_mod[[t]]<-y[[t]][vec1]
  
}
y_full <- as.matrix(do.call(cbind,y))
y_mod_full <- do.call(cbind,y_mod)
data <-y_full[,1:(nobs + 1)]


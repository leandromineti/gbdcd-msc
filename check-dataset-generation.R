rm(list=ls())

# Import functions
source("scripts/gbdcd.R") 
source("scripts/partition.R") 

load("data/cluster.RData")  # Load basic cluster structure

# Testando os dados
plot(y ~ x, data=grid, pch=15, col="light blue", cex=2)

# Desenha as conexoes do grid (ok!)
for(cont in 1:nrow(viz)){
  segments(x0=grid$x[ viz[cont,1] ],
           y0=grid$y[ viz[cont,1] ],
           x1=grid$x[ viz[cont,2] ],
           y1=grid$y[ viz[cont,2] ],
           col="orange", lwd=2)
}

points(y ~ x, data=grid, pch=15, col="light blue", cex=2.5)
text(x=grid$x, y=grid$y, labels=grid$idx)

centros   <- c(36, 37)
#centros   <- c(1, 64, 8, 57)
#centros   <- c(1, 64, 8, 57, 28)

neigh     <- viz
partition <- rcpp_partition(neigh, centros)
partition <- as.numeric( as.factor(partition) )

table(partition)
image(x=unique(grid$x), y=unique(grid$y),
      z=matrix(partition, nrow=8), col=c("gray","dark gray","green","yellow"),
      xlab="x grid", ylab="y grid")
text(x=grid$x, y=grid$y, labels=grid$idx)

## - - - - - - - - - - - - - - - - - - - - - - - 
## Creates a data set

# Define the number of simulations
n_simulations <- 10000

# Define simulation parameters
sigma <- 0.5 # Noise standard deviation
Z     <- abs(qnorm(0.95))

xa <- seq(0, 1, length = 32)  
xb <- seq(-1, 0, length = 32)
k  <- sqrt( 1/sum(xa^2) )

# Betas - for the linear regression equations
beta_a <- + Z*k*sigma*sqrt(2)
beta_b <- - Z*k*sigma*sqrt(2)

ya <- 0 + beta_a*xa + rnorm(length(xa), sd=sigma)
yb <- 0 + beta_b*xb + rnorm(length(xb), sd=sigma)

plot(c(ya,yb) ~ c(xa,xb), pch=19, xlab="[xb, xa]", ylab="[yb, ya]"); grid()
lines(beta_a*xa ~ xa, col="red", lwd=2)
lines(beta_b*xb ~ xb, col="red", lwd=2)

## Insere os valores ajustados
bA     <- coef( lm(ya ~ -1 + xa)  )
bB     <- coef( lm(yb ~ -1 + xb)  )
lines(bA*xa ~ xa, col="orange", lwd=2, lty=2)
lines(bB*xb ~ xb, col="orange", lwd=2, lty=2)
print(bA - bB)

## Colocando os dados de volta no grid...
y <- rep(NA, 64)
y[partition == 1] <- ya
y[partition == 2] <- yb
image(x=unique(grid$x), y=unique(grid$y),
      z=matrix(y, nrow=8),
      xlab="x grid", ylab="y grid")
text(x=grid$x, y=grid$y, labels=grid$idx)

## Testando a presenca dos clusters
y <- rep(NA, 64)
y[partition == 1] <- ya
y[partition == 2] <- yb

x <- rep(NA, 64)
x[partition == 1] <- xa
x[partition == 2] <- xb

summary( lm(y ~ -1 + x + x:I(partition-1)) )

# save(centros, viz, grid, y, x, file="TwoClustersMarcelo.RData")


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Aplicando um teste estatistico se os clusters sao conhecidos
## a priori...
threshold <- Z*k*sigma*sqrt(2)
vec       <- rep(NA, n_simulations)

for(s in 1:n_simulations){
  ya     <- 0 + beta_a*xa + rnorm(length(xa), sd=sigma)
  yb     <- 0 + beta_b*xb + rnorm(length(xb), sd=sigma)
  bA     <- coef( lm(ya ~ -1 + xa)  )
  bB     <- coef( lm(yb ~ -1 + xb)  )
  vec[s] <- bA - bB 
}

hist(abs(vec)); rug(abs(vec))
abline( v=threshold, lwd=2, col="red" )
sum(vec < threshold )/length(vec)


# Create datasets with one one cluster ------------------------------------
for(s_index in 1:n_simulations){
  y <- rnorm(64, mean=0, sd=1)
  save(centros, viz, grid, y, file=sprintf("datasets/one-cluster/c_one_%d.RData", s_index))
}

# Create datasets with two diagonal clusters ------------------------------

save(centros, viz, grid, y, file=sprintf("datasets/two-diagonal-clusters/c_diagonal_%d.RData", s_index))

# Create datasets with four square clusters -------------------------------

# Second script -----------------------------------------------------------

# Testando os dados
plot(y ~ x, data=grid, pch=15, col="light blue", cex=2)
# Desenha as conexoes do grid (ok!)
for(cont in 1:nrow(viz)){
  segments(x0=grid$x[ viz[cont,1] ],
           y0=grid$y[ viz[cont,1] ],
           x1=grid$x[ viz[cont,2] ],
           y1=grid$y[ viz[cont,2] ],
           col="orange", lwd=2)
}
points(y ~ x, data=grid, pch=15, col="light blue", cex=2.5)
text(x=grid$x, y=grid$y, labels=grid$idx)

centros   <- c(36, 37)

neigh     <- viz
partition <- RcppPartition(neigh, centros)
partition <- as.numeric( as.factor(partition) )

table(partition)
image(x=unique(grid$x), y=unique(grid$y),
      z=matrix(partition, nrow=8), col=c("gray","dark gray"),
      xlab="x grid", ylab="y grid")
text(x=grid$x, y=grid$y, labels=grid$idx)

## - - - - - - - - - - - - - - - - - - - - - - - 
## Creates a data set

# Define simulation parameters
sigma <- 0.5    # Noise standard deviation
Z     <- abs(qnorm(0.95))

xa <- seq(0, 1, length = 32)  
xb <- seq(-1, 0, length = 32)
k  <- sqrt( 1/sum(xa^2) )    #k  <- sqrt( 1/sum(xb^2) )  # Redundant

# Betas - for the linear regression equations
beta_a <- + Z*k*sigma*sqrt(2)
beta_b <- - Z*k*sigma*sqrt(2)

ya <- 0 + beta_a*xa + rnorm(length(xa), sd=sigma)
yb <- 0 + beta_b*xb + rnorm(length(xb), sd=sigma)

plot(c(ya,yb) ~ c(xa,xb), pch=19, xlab="[xb, xa]", ylab="[yb, ya]"); grid()
lines(beta_a*xa ~ xa, col="red", lwd=2)
lines(beta_b*xb ~ xb, col="red", lwd=2)

## Insere os valores ajustados
bA     <- coef( lm(ya ~ -1 + xa)  )
bB     <- coef( lm(yb ~ -1 + xb)  )
lines(bA*xa ~ xa, col="orange", lwd=2, lty=2)
lines(bB*xb ~ xb, col="orange", lwd=2, lty=2)
print(bA - bB)

## Colocando os dados de volta no grid...
y <- rep(NA, 64)
y[partition == 1] <- ya
y[partition == 2] <- yb
image(x=unique(grid$x), y=unique(grid$y),
      z=matrix(y, nrow=8),
      xlab="x grid", ylab="y grid")
text(x=grid$x, y=grid$y, labels=grid$idx)

## Testando a presenca dos clusters
y <- rep(NA, 64)
y[partition == 1] <- ya
y[partition == 2] <- yb

x <- rep(NA, 64)
x[partition == 1] <- xa
x[partition == 2] <- xb

summary( lm(y ~ -1 + x + x:I(partition-1)) )

# save(centros, viz, grid, y, x, file="TwoClustersMarcelo.RData")


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## TESTANDO O PASSO DE UPDATE
plot(ya ~ xa, col="blue", pch=19); abline( lm(ya ~ xa), lty=1, col="red")
summary( lm(ya ~ xa) )
sig.hat <- sqrt( anova( lm(ya ~ xa) )["Residuals", "Mean Sq"] )

X       <- cbind(1, xa)
Y       <- as.matrix(ya)
lambda0 <- 1e-9
a0      <- 2.1
b0      <- 1.1
sigma2  <- 0.5^2


require(mvtnorm) # The Multivariate Normal Distribution


steps     <- 10000
vec.sigma <- rep(NA, steps)
mat.beta  <- matrix(NA, nrow=steps, ncol=2) 

for(st in 1:steps){
  # Nova proposta para o vetor Beta
  auX     <- solve( crossprod(X,X) + lambda0*diag(rep(1,2)) )
  MuBeta  <- auX%*%crossprod(X,ya)
  CovBeta <- sigma2*auX
  
  newBeta <- rmvnorm(n=1, mean=MuBeta, sigma=CovBeta)
  
  fit     <- X%*%t(newBeta)
  #lines(fit ~ xa, col="black", lty=3)
  
  an <- a0 + nrow(X)/2
  bn <- b0 + sum( (ya - fit)^2 )/2
  
  sigma2 <- 1/rgamma(n=1, shape=an, rate=bn)
  
  vec.sigma[st] <- sqrt(sigma2)
  mat.beta[st,] <- newBeta
}

summary(vec.sigma^2)
names(mat.beta) <- c("beta0", "beta1")
summary(mat.beta)
print( c( anova( lm(ya ~ xa) )["Residuals", "Mean Sq"],
          coef(lm(ya ~ xa)) )  )

windows()
par(mfrow=c(2,2))
coeficientes <- coef( lm(ya ~ xa) )
hist(vec.sigma); abline(v=sig.hat, col="red")
hist(mat.beta[,1]); abline(v=coeficientes[1], col="red")
hist(mat.beta[,2]); abline(v=coeficientes[2], col="red")



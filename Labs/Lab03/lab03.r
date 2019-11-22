library(gclus)

data(bank)
bank = as.matrix(bank)
y <- bank[,1]
X <- cbind(1,bank[,2:5])
q <- dim(X)[2] 
n <- dim(X)[1]

#likelihood

lik = function(X, y, beta){
  a = pnorm(X%*%matrix(beta,ncol = 1))
  b = 1-a
  return(prod(a^y*b^(1-y)))
}

#posterior is k*likelihood -> stays like this

posterior = lik

#M-H

#beta_star ~ N_{q+1} (beta_{k-1}, tau^2 sigma_hat), tau^2 is a constant

metro_hastings = function(nsim, tau, y, X){
  beta = matrix(0, nrow=nsim, ncol=ncol(X))
  initializer_model = glm(y~X, family=binomial(link="probit"))
  
  sigma_hat = summary(initializer_model)$cov.unscaled
  
  beta[1,] = summary(initializer_model)$coefficients[,1]
  
  for(s in 2:nsim){
    beta_proposal = mvtnorm::rmvnorm(1, mean = beta[s-1,], sigma = tau^2*sigma_hat)
    accept_prob = min(1, posterior(X, y, beta_proposal)/posterior(X, y, beta[s-1,]))
    u = runif(1)
    if(u<accept_prob){
      beta[s,] = beta_proposal
    }else{
      beta[s,] = beta[s-1,]
    }
  }
  return(beta)
}

mh1 = metro_hastings(1e4, .1, y, X)
mh2 = metro_hastings(1e4, 1, y, X)
mh3 = metro_hastings(1e4, 4, y, X)
mh = list(1 = mh1, mh2 = mh2, mh3 = mh3)


#traceplots
par(mfrow = c(3,5))
for(m in 1:3)
  for(i in 1:5)
    plot(mh[[m]][,i], type='l', ylab = paste0("beta",i-1))

#histogram
par(mfrow = c(3,5))
for(m in 1:3)
  for(i in 1:5)
    hist(mh[[m]][,i], breaks=100)

#ACF
par(mfrow = c(3,5))
for(m in 1:3)
  for(i in 1:5)
    acf(mh[[m]][,i], breaks=100)

#Means
for(m in 1:3){
  for(i in 1:5){
    a <- cumsum(mh[[m]][,i])/1:1e4
    plot(a, type="l",ylab="Cumulative mean plot", xlab="Iteration")
    abline(h=mean(mh[[m]][i]), col="firebrick3", lty=2)
  }
}


#gibbs 

gibbs = function(nsim, X, y){
  p = ncol(X)
  n = nrow(X)
  beta = matrix(0, nrow=nsim, ncol=p)
  
  beta[1,] = summary(glm(y~X, family=binomial(link="probit")))$coefficients[,1]
  z = X%*%beta[1,] + rnorm(n)
  y = ifelse(z>0, 1, 0)
  
  XtX_1 = solve(t(X)%*%X)
  XtX_1Xt = XtX_1 %*% t(X)
  
  
  for(s in 2:nsim){
    z = ifelse(y==1,
               truncnorm::rtruncnorm(1, a = 0, mean = X%*%beta[s-1,], sd = 1),
               truncnorm::rtruncnorm(1, b = 0, mean = X%*%beta[s-1,], sd = 1))
    beta[s,] = mvtnorm::rmvnorm(1, mean = XtX_1Xt %*% z, sigma = XtX_1)
  }
  
  return(beta)
  
}

gibbs(1000, X, y)

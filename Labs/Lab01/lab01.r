#gioia.dicredico@deams.units.it


samples = list()
n = c(20, 100, 1000)
p = .5
index = 1
for (i in n){
  samples[index]=rpois(i, p)
  index = index + 1
}

curve(Vectorize())


g1 = function(x) dgamma(x, 1, .5)
#h1 = Vectorize(g1)
g2 = function(x) dgamma(x, 1, 2)
#h2 = Vectorize(g2)
g3 = function(x) dgamma(x, 1, 10)
#h3 = Vectorize(g3)
g4 = function(x) dgamma(x, 2, 2)
#h4 = Vectorize(g4)
g5 = function(x) 1/sqrt(x)
priors = c(g1, g2, g3, g4, g5)
curve(g1, xlim = c(0, 3), ylim = c(0, 3))
curve(g2, add = T, col = "red")
curve(g3, add = T, col = "orange")
curve(g4, add = T, col = "blue")
curve(g5, add = T, col = "chocolate4", lwd = 2, lty = 5)

sample1 = rpois(20, .5)
sample2 = rpois(100, .5)
sample3 = rpois(1000, .5)

posterior = function(sample, prior_alpha, prior_beta){
  sumy = sum(sample)
  n = length(sample)
  return ( function(x)dgamma(x, prior_alpha + sumy, prior_beta + n))
}

prior_alpha = c(1,1,1,2,0)
prior_beta = c(.5, 2, 10, 2, 0.5)
posteriors = list()#(nrow = length(prior_alpha), ncol = 3, data = NA)

for(i in 1:length(prior_alpha)){
  k = c(posterior(sample1, prior_alpha[i], prior_beta[i]),
        posterior(sample2, prior_alpha[i], prior_beta[i]),
        posterior(sample3, prior_alpha[i], prior_beta[i]))
  posteriors[[i]] = k
  
}

curve(posteriors[[1]][[1]])

for(i in 1:length(prior_alpha)){
  print("MEANS")
  print((prior_alpha[i]+sum(sample1))/(length(sample1) + prior_beta[i]))
  print((prior_alpha[i]+sum(sample2))/(length(sample2) + prior_beta[i]))
  print((prior_alpha[i]+sum(sample3))/(length(sample3) + prior_beta[i]))
  print("MAP")
  print((prior_alpha[i]+sum(sample1)-1)/(length(sample1) + prior_beta[i]))
  print((prior_alpha[i]+sum(sample2)-1)/(length(sample2) + prior_beta[i]))
  print((prior_alpha[i]+sum(sample3)-1)/(length(sample3) + prior_beta[i]))
}


simpost = rgamma(10000, sum(sample1) + .5, 20)

q = quantile(simpost, c(.025,.975))

q_correct = density(simpost)

MAP = max(q_correct$y)
mode = q_correct$x[q_correct$y==MAP]

hist(simpost)
abline(v = q[1], col = "red", lty = 5)
abline(v = q[2], col = "red", lty = 5)



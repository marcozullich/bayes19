library(rstan)
library(dplyr)
library(lubridate)
library(ggplot2)
library(bayesplot)

pest_data = readRDS(file = "pest_data.RDS")


theme_set(bayesplot::theme_default())

# seed for R's pseudo-RNGs, not Stan's
set.seed(1123) 

N_buildings <- length(unique(pest_data$building_id))
N_buildings

pest_data %>% pull(date) %>% as.numeric 


corrplot::corrplot.mixed(cor(pest_data %>% select(-date)))

#simple poisson model
comp_model_P <- stan_model('simple_poisson_regression2.stan')

stan_dat_simple <- list(
  N = nrow(pest_data), 
  x = pest_data$traps,
  y = pest_data$complaints
)
str(stan_dat_simple)

fit_P_real_data <- sampling(comp_model_P, data = stan_dat_simple)
print(fit_P_real_data, pars = c('alpha','beta'))
# n_eff -> number of effective draws that can be considered independent
# rhat -> did parameter converge?

mcmc_hist(as.matrix(fit_P_real_data, pars = c('alpha','beta')))
#correlation is good since it shows there's relation between intercept and beta
#if it weren't like this, it would be an indicator of little correlation


mcmc_scatter(as.matrix(fit_P_real_data, pars = c('alpha','beta')), alpha = 0.2)


## posterior predictive checking
y_rep <- as.matrix(fit_P_real_data, pars = "y_rep")
ppc_dens_overlay(y = stan_dat_simple$y, y_rep[1:200,])
#low no. of complaints -> underestimate
#medium no. of complaints -> largely overestimate
#then -> ~ similar
#--> poisson model not OK

## standardised residuals of the observed vs predicted number of complaints
mean_y_rep <- colMeans(y_rep)
std_resid <- (stan_dat_simple$y - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)
#model underestimates, bc there're a lot of pos. outliers

#################### insert model with one more covariate + offset
ggplot(pest_data, aes(x = log(total_sq_foot), y = log1p(complaints))) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

## add the two variables to the list of the data
stan_dat_simple$offset <- log(pest_data$total_sq_foot/1e4)
stan_dat_simple$x2 <- pest_data$live_in_super
str(stan_dat_simple)

## compile the model
comp_model_P_mult <- stan_model('multiple_poisson_regression.stan')
fit_model_P_mult_real <- sampling(comp_model_P_mult, data = stan_dat_simple)
y_rep <- as.matrix(fit_model_P_mult_real, pars = "y_rep")
ppc_dens_overlay(stan_dat_simple$y, y_rep[1:200,])
#still underestimating
mean_y_rep <- colMeans(y_rep)
std_resid <- (stan_dat_simple$y - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)

ppc_intervals(
  y = stan_dat_simple$y, 
  yrep = y_rep,
  x = stan_dat_simple$x
) + 
  labs(x = "Number of traps", y = "Number of complaints")

################## BINEG MODEL
comp_model_NB <- stan_model('multiple_NB_regression.stan')
fitted_model_NB <- sampling(comp_model_NB, data = stan_dat_simple)

samps_NB <- rstan::extract(fitted_model_NB)
## predictions vs. the data
y_rep <- samps_NB$y_rep
ppc_dens_overlay(stan_dat_simple$y, y_rep[1:200,])

#std res
mean_inv_phi <- mean(samps_NB$inv_phi)
mean_y_rep <- colMeans(y_rep)
std_resid <- (stan_dat_simple$y - mean_y_rep) / sqrt(mean_y_rep + mean_y_rep^2*mean_inv_phi)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)

#posterior predictive check
ppc_stat_grouped(
  y = stan_dat_simple$y, 
  yrep = y_rep, 
  group = pest_data$building_id, 
  stat = 'mean',
  binwidth = 0.2
)
#what does this chart tell us? there's a nice variability in mean of complaints for each building
#so we can sketch a hierarchical model
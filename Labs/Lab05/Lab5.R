library(rstan)
library(dplyr)
library(lubridate)
library(ggplot2)
library(bayesplot)
library(loo)

theme_set(bayesplot::theme_default())
rstan_options(auto_write=TRUE) #during compilation, stan saves file .rds translating the
#model, so that there's no need to recompile it again

pest_data <- readRDS('../Lab04/pest_data.RDS')
str(pest_data)

N_buildings <- length(unique(pest_data$building_id))

N_months <- length(unique(pest_data$date))

# Add some IDs for building and month
pest_data <- pest_data %>%
  mutate(
    building_fac = factor(building_id, levels = unique(building_id)),
    building_idx = as.integer(building_fac),
    ids = rep(1:N_months, N_buildings),
    mo_idx = lubridate::month(date)
  )

# Center and rescale the building specific data
building_data <- pest_data %>%
  select(
    building_idx,
    live_in_super,
    age_of_building,
    total_sq_foot,
    average_tenant_age,
    monthly_average_rent
  ) %>%
  unique() %>%
  arrange(building_idx) %>%
  select(-building_idx) %>%
  scale(scale=FALSE) %>%
  as.data.frame() %>%
  mutate( # scale by constants
    age_of_building = age_of_building / 10,
    total_sq_foot = total_sq_foot / 10000,
    average_tenant_age = average_tenant_age / 10,
    monthly_average_rent = monthly_average_rent / 1000
  ) %>%
  as.matrix()

str(building_data)

stan_dat_hier <-
  with(pest_data,
       list(complaints = complaints,
            traps = traps,
            N = length(traps),
            J = N_buildings,
            log_sq_foot = log(pest_data$total_sq_foot/1e4),
            building_data = building_data[,-3],
            mo_idx = as.integer(as.factor(date)),
            K = 4,
            building_idx = building_idx
       )
  )
str(stan_dat_hier)

comp_model_NB_hier <- stan_model('hier_NB_regression.stan')

fitted_model_NB_hier <- sampling(comp_model_NB_hier, data = stan_dat_hier,
                                 chains = 4, cores = 4, iter = 4000)
# divergence due to difficult sampling in the parameters space

samps_hier_NB <- rstan::extract(fitted_model_NB_hier)
print(fitted_model_NB_hier, pars = c('sigma_mu','beta','alpha','phi','mu'))
# n_eff -> number of independent simulations, if high, inference may be problematic

mcmc_trace(
  as.array(fitted_model_NB_hier,pars = 'sigma_mu'),
  np = nuts_params(fitted_model_NB_hier),
  window = c(500,1000)
)

# DIVERGENCY PLOT
# assign to object so we can compare to another plot later
scatter_with_divs <- mcmc_scatter(
  as.array(fitted_model_NB_hier),
  pars = c("mu[4]", 'sigma_mu'),
  transform = list('sigma_mu' = "log"),
  np = nuts_params(fitted_model_NB_hier)
)
scatter_with_divs
#COncentration on lower part of parameters space

N_sims <- 1000
log_sigma <- rep(NA, N_sims)
theta <- rep(NA, N_sims)
for (j in 1:N_sims) {
  log_sigma[j] <- rnorm(1, mean = 0, sd = 1)
  theta[j] <- rnorm(1, mean = 0, sd = exp(log_sigma[j]))
}
draws <- cbind("mu" = theta, "log(sigma_mu)" = log_sigma)
mcmc_scatter(draws)
#when values for sigma are high, moving in the space is easier
#sol:
#define an auxiliary parameter such that sampling is easier


parcoord_with_divs <- mcmc_parcoord(
  as.array(fitted_model_NB_hier, pars = c("sigma_mu", "mu")),
  np = nuts_params(fitted_model_NB_hier)
)
parcoord_with_divs
#what parameter is creating troubles??

#refit
comp_model_NB_hier_rep <- stan_model('hier_NB_regression_reparam.stan')

fitted_model_NB_hier_rep <- sampling(comp_model_NB_hier_rep, data = stan_dat_hier,
                                 chains = 4, cores = 4, iter = 4000)

print(fitted_model_NB_hier_rep, pars = c('sigma_mu','beta','alpha','phi','mu'))

scatter_no_divs <- mcmc_scatter(
  as.array(fitted_model_NB_hier_rep),
  pars = c("mu[4]", 'sigma_mu'),
  transform = list('sigma_mu' = "log"),
  np = nuts_params(fitted_model_NB_hier_rep)
)
bayesplot_grid(scatter_with_divs, scatter_no_divs,
               grid_args = list(ncol = 2), ylim = c(-11, 1))
#better exploration within the lower parameter space




parcoord_no_divs <- mcmc_parcoord(
  as.array(fitted_model_NB_hier_rep, pars = c("sigma_mu", "mu")),
  np = nuts_params(fitted_model_NB_hier_rep)
)
bayesplot_grid(parcoord_with_divs, parcoord_no_divs,
               ylim = c(-3, 3))

samps_NB_hier_rep <- rstan::extract(fitted_model_NB_hier_rep, pars = c('y_rep','inv_phi'))
y_rep <- as.matrix(fitted_model_NB_hier_rep, pars = "y_rep")
ppc_dens_overlay(stan_dat_hier$complaints, y_rep[1:200,])




ppc_stat_grouped(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  group = pest_data$building_id,
  stat = 'mean',
  binwidth = 0.5
)


## modello gerarchico con pendenza variabile

stan_dat_hier <- readRDS('pest_data_longer_stan_dat.RDS')

comp_model_NB_hier_slopes <- stan_model('hier_NB_regression_ncp_slopes_mod.stan')

fitted_model_NB_hier_slopes <-
  sampling(
    comp_model_NB_hier_slopes,
    data = stan_dat_hier,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )


mcmc_hist(
  as.matrix(fitted_model_NB_hier_slopes, pars = "sigma_kappa"),
  binwidth = 0.005
) #sigma_kappa values small but not concentrated around zero -> effective partial pooling

# 10 - add temporal component - is there a difference between months?

select_vec <- which(stan_dat_hier$mo_idx %in% 1:12)

ppc_stat_grouped(
  y = stan_dat_hier$complaints[select_vec],
  yrep = y_rep[,select_vec],
  group = stan_dat_hier$mo_idx[select_vec],
  stat = 'mean'
) + xlim(0, 11)

comp_model_NB_hier_mos <- stan_model('hier_NB_regression_ncp_slopes_mod_mos.stan')
fitted_model_NB_hier_mos <- sampling(comp_model_NB_hier_mos, data = stan_dat_hier, chains = 4, cores = 4, control = list(adapt_delta = 0.9))

y_rep <- as.matrix(fitted_model_NB_hier_mos, pars = "y_rep")
ppc_dens_overlay(
  y = stan_dat_hier$complaints,
  yrep = y_rep[1:200,]
)

#LOOCV
comp_model_NB_hier_slopes_llik <- stan_model('hier_NB_regression_ncp_slopes_mod.stan')

fitted_model_NB_hier_slopes_llik <-
  sampling(
    comp_model_NB_hier_slopes_llik,
    data = stan_dat_hier,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )

comp_model_NB_hier_mos_llik <- stan_model('hier_NB_regression_ncp_slopes_mod_mos.stan')

fitted_model_NB_hier_mos_llik <-
  sampling(
    comp_model_NB_hier_mos_llik,
    data = stan_dat_hier,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )



log_lik_slopes <- extract_log_lik(fitted_model_NB_hier_slopes_llik)
loo_slopes <- loo(log_lik_slopes)
loo_slopes

log_lik_mos <- extract_log_lik(fitted_model_NB_hier_mos_llik)
loo_mos <- loo(log_lik_mos)
loo_mos

compare(loo_slopes, loo_mos)

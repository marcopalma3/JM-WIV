###How to use: 
###Rscript HPC_Simulation_JMWIV.R





cat("Start simulation.\n\n")
library(pacman)
p_load(rstan, dplyr, magrittr, gridExtra, formula.tools)
p_load(stringr, nlme, survival, JMbayes2, here, rlist, lme4, mvQuad, loo, rstanarm)

library(rstanjmwiv)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(digits = 5)

job_id_string <- Sys.getenv("SLURM_ARRAY_JOB_ID", unset = 1)
array_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = 1)
cat("The job ID string is", job_id_string, ".\n\n")
cat("The array ID string is", array_id_string, ".\n\n")

out_filename_jmbayes2 <- paste0("Simulations/Results/SIM_jmbayes2_", as.numeric(job_id_string), "_", str_pad(as.numeric(array_id_string), 3, pad = "0"), ".txt")
out_filename_jmwiv <- paste0("Simulations/Results/SIM_jmwiv_", as.numeric(job_id_string), "_", str_pad(as.numeric(array_id_string), 3, pad = "0"), ".txt")
out_filename_rstanarm <- paste0("Simulations/Results/SIM_rstanarm_", as.numeric(job_id_string), "_", str_pad(as.numeric(array_id_string), 3, pad = "0"), ".txt")
out_filename_timing <- paste0("Simulations/Results/TIMING_", as.numeric(job_id_string), "_", str_pad(as.numeric(array_id_string), 3, pad = "0"), ".txt")
out_filename_loo_jmwiv <- paste0("Simulations/Results/LOO_jmwiv_", as.numeric(job_id_string), "_", str_pad(as.numeric(array_id_string), 3, pad = "0"), ".txt")
out_filename_loo_jmbayes2 <- paste0("Simulations/Results/LOO_jmbayes2_", as.numeric(job_id_string), "_", str_pad(as.numeric(array_id_string), 3, pad = "0"), ".txt")
out_filename_loo_rstanarm <- paste0("Simulations/Results/LOO_rstanarm_", as.numeric(job_id_string), "_", str_pad(as.numeric(array_id_string), 3, pad = "0"), ".txt")

cat("The JMbayes2 output will be saved in", here(out_filename_jmbayes2), ".\n\n")
cat("The JM-WIV output will be saved in", here(out_filename_jmwiv), ".\n\n")
cat("The rstanarm output will be saved in", here(out_filename_rstanarm), ".\n\n")
cat("The timing results will be saved in", here(out_filename_timing), ".\n\n")

timing_results <- data.frame(
  model = character(),
  user_time = numeric(),
  system_time = numeric(),
  elapsed_time = numeric(),
  stringsAsFactors = FALSE
)

ranef_corr <- matrix(c(1, 0.129, 0.502, 0.019,
                       0.129, 1, -0.006, 0.243,
                       0.502, -0.006, 1, 0.491,
                       0.019, 0.243, 0.491, 1),
                     byrow = F,
                     ncol = 4)
sd_vec <- c(0.81, 0.44, 0.52, 0.16) * diag(4)
ranef_mat <- sd_vec %*% ranef_corr %*% t(sd_vec)

e_beta_input <- c(0.93, -2.3)
a_beta_input <- c(-2.24, 1.9, 0.55, 0.35)
y1_mu_fixed_input <- c(2.19, -0.04)
y2_mu_fixed_input <- c(1.04, 0.35)

input_vec <- c(e_beta_input,
               a_beta_input[c(1,3)],
               y1_mu_fixed_input,
               y2_mu_fixed_input)

cat("Simulate data.\n\n")

sim <- simulate_jmwiv(seed = as.numeric(job_id_string) + as.numeric(array_id_string),
                      nsub = 1000,
                      assoc = 'LP',
                      lambda = 5,
                      g_shape = 0,
                      e_beta = e_beta_input,
                      a_beta = a_beta_input,
                      distObs = 1,
                      max_time = 10,
                      ranef_covmat = ranef_mat,
                      y1mu_fixed = y1_mu_fixed_input,
                      y1sigma_fixed = c(-1.36, -0.02),
                      y2mu_fixed = y2_mu_fixed_input,
                      y2sigma_fixed = c(0.16, 0.01))

head(sim$dataLong)
head(sim$dataEvent)

entry_times <- sim$dataLong %>%
  group_by(ID) %>%
  summarise(entry_time = min(visit_times)) %>%
  ungroup()

sim$dataEvent <- sim$dataEvent %>%
  left_join(entry_times, by = "ID")

sim[[1]]$ID <- as.numeric(sim[[1]]$ID)
sim[[2]]$ID <- as.numeric(sim[[2]]$ID)

fL1 <- list(formula(y1 ~ visit_times + (1|ID)),
            formula(sigma ~ visit_times + (1|ID)))
fL2 <- list(formula(y2 ~ visit_times + (1|ID)),
            formula(sigma ~ visit_times + (1|ID)))
fEvent <- survival::Surv(e_time, e_status) ~ binary_cov + norm_cov

prior_input <- list(
  a_prior_scale = rep(1, 4),
  e_prior_scale_for_aux = rep(3, 6),
  b_prior_scale = rep(3, 4),
  b_prior_df = rep(3, 4),
  b_prior_regularization = 3
)

cat("Start rstanarm.\n\n")

time_rstanarm <- system.time({fit_rstanarm <- rstanarm::stan_jm(
  formulaLong = list(y1 ~ visit_times + (1|ID),
                     y2 ~ visit_times + (1|ID)),
  dataLong = sim$dataLong,
  formulaEvent = survival::Surv(e_time, e_status) ~ binary_cov + norm_cov,
  dataEvent = sim$dataEvent,
  time_var = "visit_times",
  id_var = "ID",
  basehaz = "bs",
  basehaz_ops = list(df = 6),
  qnodes = 15,
  assoc = "etavalue",
  max_treedepth = 15,
  chains = 2,
  cores = 2,
  iter = 6000,
  seed = as.numeric(job_id_string) + as.numeric(array_id_string)
)})

print(summary(fit_rstanarm))

posterior_samples <- as.matrix(fit_rstanarm)

rstanarm_mcmc <- data.frame(
  e_beta_binary = posterior_samples[, "Event|binary_cov"],
  e_beta_norm = posterior_samples[, "Event|norm_cov"],
  a_beta_y1 = posterior_samples[, "Assoc|Long1|etavalue"],
  a_beta_y2 = posterior_samples[, "Assoc|Long2|etavalue"],
  y1mu_Intercept = posterior_samples[, "Long1|(Intercept)"],
  y1mu_visit_times = posterior_samples[, "Long1|visit_times"],
  y2mu_Intercept = posterior_samples[, "Long2|(Intercept)"],
  y2mu_visit_times = posterior_samples[, "Long2|visit_times"],
  tau_y1 = sqrt(posterior_samples[, "Sigma[ID:Long1|(Intercept),Long1|(Intercept)]"]),
  tau_y2 = sqrt(posterior_samples[, "Sigma[ID:Long2|(Intercept),Long2|(Intercept)]"]),
  log_sigma_y1 = log(posterior_samples[, "Long1|sigma"]),
  log_sigma_y2 = log(posterior_samples[, "Long2|sigma"])
)

write.table(round(rstanarm_mcmc, 5), file = here(out_filename_rstanarm), sep = "\t",
            row.names = FALSE, col.names = TRUE)
cat(c(unname(round(colMeans(rstanarm_mcmc), 5)), "rstanarm"), "\n",
    file = here("Simulations/SIM_summary_rstanarm.txt"), append = TRUE)

model_summary <- summary(fit_rstanarm)
CI_rstanarm <- model_summary[c("Event|binary_cov", "Event|norm_cov",
                               "Assoc|Long1|etavalue", "Assoc|Long2|etavalue",
                               "Long1|(Intercept)", "Long1|visit_times",
                               "Long2|(Intercept)", "Long2|visit_times",
                               "Sigma[ID:Long1|(Intercept),Long1|(Intercept)]",
                               "Sigma[ID:Long2|(Intercept),Long2|(Intercept)]",
                               "Long1|sigma", "Long2|sigma"), c('2.5%', '97.5%')]

input_vec_rstanarm <- c(e_beta_input, a_beta_input[c(1,3)], y1_mu_fixed_input, y2_mu_fixed_input,
                        c(0.81^2, 0.52^2), exp(-1.36), exp(0.16))
coverage_rstanarm <- as.numeric(input_vec_rstanarm > CI_rstanarm[,1] & input_vec_rstanarm < CI_rstanarm[,2])
cat(c(coverage_rstanarm, "rstanarm"), "\n",
    file = here("Simulations/SIM_COVERAGE_rstanarm.txt"), append = TRUE)

log_lik_matrix_rstanarm <- rstanarm::log_lik(fit_rstanarm)
loo_rstanarm <- loo(log_lik_matrix_rstanarm)
print(loo_rstanarm)
saveRDS(loo_rstanarm, file = here(out_filename_loo_rstanarm))

timing_results <- rbind(timing_results,
                        data.frame(model = "rstanarm",
                                   user_time = time_rstanarm["user.self"],
                                   system_time = time_rstanarm["sys.self"],
                                   elapsed_time = time_rstanarm["elapsed"]))

cat("Start JM-WIV.\n\n")

time_jmwiv <- system.time({fit_jmwiv <- rstanjmwiv::jmwiv_stan(
  formulaLong1 = fL1,
  formulaLong2 = fL2,
  dataLong = sim$dataLong,
  id_var = "ID",
  time_varLong = "visit_times",
  formulaEvent = fEvent,
  dataEvent = sim$dataEvent,
  a_K = 4,
  assoc = "LP",
  basehaz = "bs",
  basehaz_aux = list(df = 6, knots = NULL, degree = 3),
  qnodes = 15L,
  prior_list = prior_input,
  cores = 2,
  chains = 2,
  warmup = 3000,
  iter = 6000,
  seed = as.numeric(job_id_string) + as.numeric(array_id_string),
  init_r = 1,
  control = list(max_treedepth = 15)
)})

print(fit_jmwiv,
      pars = c("y1mu_Intercept","y1sigma_Intercept","y2mu_Intercept","y2sigma_Intercept",
               "y1mu_beta","y1sigma_beta","y2mu_beta","y2sigma_beta"),
      probs = c(0.025,0.975))
print(fit_jmwiv, pars = "a_beta", probs = c(0.025,0.975))
print(fit_jmwiv, pars = c("e_Intercept", "e_beta"), probs = c(0.025,0.975))
print(fit_jmwiv, pars = "b_sd", probs = c(0.025,0.975))

JMWIV_mcmc <- cbind(
  rstan::extract(fit_jmwiv, pars = c("e_beta"))$e_beta,
  rstan::extract(fit_jmwiv, pars = c("a_beta"))$a_beta[,c(1,3)],
  rstan::extract(fit_jmwiv, pars = c("y1mu_Intercept"))[[1]],
  rstan::extract(fit_jmwiv, pars = c("y1mu_beta"))[[1]],
  rstan::extract(fit_jmwiv, pars = c("y2mu_Intercept"))[[1]],
  rstan::extract(fit_jmwiv, pars = c("y2mu_beta"))[[1]],
  rstan::extract(fit_jmwiv, pars = c("b_sd"))[[1]],
  rstan::extract(fit_jmwiv, pars = c("a_beta"))$a_beta[,c(2,4)],
  rstan::extract(fit_jmwiv, pars = c("y1sigma_Intercept"))[[1]],
  rstan::extract(fit_jmwiv, pars = c("y1sigma_beta"))[[1]],
  rstan::extract(fit_jmwiv, pars = c("y2sigma_Intercept"))[[1]],
  rstan::extract(fit_jmwiv, pars = c("y2sigma_beta"))[[1]],
  rstan::extract(fit_jmwiv, pars = c("e_tau"))[[1]]
)

JMWIV_mcmc <- data.frame(JMWIV_mcmc)

write.table(round(JMWIV_mcmc, 5), file = here(out_filename_jmwiv), sep = "\t",
            row.names = FALSE, col.names = TRUE)
cat(c(unname(round(colMeans(JMWIV_mcmc), 5)), "JM-WIV"), "\n",
    file = here("Simulations/SIM_summary_JMWIV.txt"), append = TRUE)

CI_JMWIV <- summary(fit_jmwiv)$summary[c("e_beta[1]", "e_beta[2]",
                                          "a_beta[1]", "a_beta[3]",
                                          "y1mu_Intercept","y1mu_beta[1]",
                                          "y2mu_Intercept","y2mu_beta[1]",
                                          "b_sd[1]", "b_sd[2]", "b_sd[3]", "b_sd[4]",
                                          "a_beta[2]", "a_beta[4]",
                                          "y1sigma_Intercept","y1sigma_beta[1]",
                                          "y2sigma_Intercept","y2sigma_beta[1]"), c('2.5%', '97.5%')]

input_vec_JMWIV <- c(input_vec, c(0.81, 0.44, 0.52, 0.16), a_beta_input[c(2,4)],
                     c(-1.36, -0.02), c(0.16, 0.01))
coverage_JMWIV <- as.numeric(input_vec_JMWIV > CI_JMWIV[,1] & input_vec_JMWIV < CI_JMWIV[,2])
cat(c(coverage_JMWIV, "JMWIV"), "\n",
    file = here("Simulations/SIM_COVERAGE_JMWIV.txt"), append = TRUE)

log_lik_matrix_jmwiv <- rstan::extract(fit_jmwiv, pars = "log_lik")$log_lik
loo_jmwiv <- loo(log_lik_matrix_jmwiv)
print(loo_jmwiv)
saveRDS(loo_jmwiv, file = here(out_filename_loo_jmwiv))

timing_results <- rbind(timing_results,
                        data.frame(model = "JMWIV",
                                   user_time = time_jmwiv["user.self"],
                                   system_time = time_jmwiv["sys.self"],
                                   elapsed_time = time_jmwiv["elapsed"]))

cat("Start JMbayes2.\n\n")

fm1 <- lme(fixed = y1 ~ visit_times, random = ~ 1 | ID, data = sim$dataLong)
fm2 <- lme(fixed = y2 ~ visit_times, random = ~ 1 | ID, data = sim$dataLong)
Mixed <- list(fm1, fm2)

fCox1 <- coxph(Surv(e_time, e_status) ~ binary_cov + norm_cov, sim$dataEvent)

time_jmbayes2 <- system.time({fit_jmbayes2 <- jm(fCox1, Mixed, time_var = "visit_times",
                                                  n_chains = 2L, n_iter = 6000L, n_burnin = 3000L,
                                                  cores = 2L,
                                                  seed = as.numeric(job_id_string) + as.numeric(array_id_string),
                                                  save_logLik_contributions = TRUE)})

summary(fit_jmbayes2)

JMbayesfit_mcmc <- cbind(
  fit_jmbayes2$mcmc %>% .$gammas %>% .[[1]],
  fit_jmbayes2$mcmc %>% .$alphas %>% .[[1]],
  fit_jmbayes2$mcmc %>% .$betas1 %>% .[[1]],
  fit_jmbayes2$mcmc %>% .$betas2 %>% .[[1]],
  fit_jmbayes2$mcmc %>% .$sigmas %>% .[[1]] %>% log(.),
  fit_jmbayes2$mcmc %>% .$D %>% .[[1]] %>% .[, -2] %>% sqrt()
)

colnames(JMbayesfit_mcmc)[1:4] <- paste0("e.", colnames(JMbayesfit_mcmc)[1:4])
colnames(JMbayesfit_mcmc)[5:6] <- paste0("y1.", colnames(JMbayesfit_mcmc)[5:6])
colnames(JMbayesfit_mcmc)[7:8] <- paste0("y2.", colnames(JMbayesfit_mcmc)[7:8])
colnames(JMbayesfit_mcmc)[9:10] <- paste0("logsigma.", colnames(JMbayesfit_mcmc)[9:10])
colnames(JMbayesfit_mcmc)[11:12] <- c("b_sd1", "b_sd3")

JMbayesfit_mcmc <- data.frame(JMbayesfit_mcmc)

write.table(round(JMbayesfit_mcmc, 5), file = here(out_filename_jmbayes2), sep = "\t",
            row.names = FALSE, col.names = TRUE)
cat(c(unname(round(colMeans(JMbayesfit_mcmc), 5)), "JMbayes2"), "\n",
    file = here("Simulations/SIM_summary_JMbayes2.txt"), append = TRUE)

CI_JMbayes2 <- rbind(
  summary(fit_jmbayes2)$Survival[, c('2.5%', '97.5%')],
  summary(fit_jmbayes2)$Outcome1[1:2, c('2.5%', '97.5%')],
  summary(fit_jmbayes2)$Outcome2[1:2, c('2.5%', '97.5%')],
  summary(fit_jmbayes2)$Outcome1[3, c('2.5%', '97.5%')],
  summary(fit_jmbayes2)$Outcome2[3, c('2.5%', '97.5%')],
  fit_jmbayes2$mcmc %>% .$D %>% .[[1]] %>% .[, -2] %>% sqrt() %>%
    apply(., 2, quantile, probs = c(0.025, 0.975)) %>% t()
)

input_vec_JMbayes2 <- c(input_vec, -1.36, 0.16, 0.81, 0.52)
coverage_JMbayes2 <- as.numeric(input_vec_JMbayes2 > CI_JMbayes2[,1] & input_vec_JMbayes2 < CI_JMbayes2[,2])
cat(c(coverage_JMbayes2, "JMbayes2"), "\n",
    file = here("Simulations/SIM_COVERAGE_JMbayes2.txt"), append = TRUE)

log_lik_matrix_jmbayes2 <- fit_jmbayes2$logLik
loo_jmbayes2 <- loo(log_lik_matrix_jmbayes2)
print(loo_jmbayes2)
saveRDS(loo_jmbayes2, file = here(out_filename_loo_jmbayes2))

timing_results <- rbind(timing_results,
                        data.frame(model = "JMbayes2",
                                   user_time = time_jmbayes2["user.self"],
                                   system_time = time_jmbayes2["sys.self"],
                                   elapsed_time = time_jmbayes2["elapsed"]))

write.table(timing_results, file = here(out_filename_timing), sep = "\t", row.names = FALSE, col.names = TRUE)

cat(paste(c(as.numeric(job_id_string), as.numeric(array_id_string),
            timing_results$elapsed_time, "\n"), collapse = "\t"),
    file = here("Simulations/TIMING_summary_all.txt"), append = TRUE)

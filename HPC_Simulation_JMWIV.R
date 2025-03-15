###How to use: 
###Rscript HPC_Simulation_JMWIV.R





cat("Start simulation.\n\n")


library(pacman)
p_load(rstan, dplyr, magrittr, gridExtra, formula.tools)
p_load(stringr, nlme, survival, JMbayes2, here, rlist, lme4, mvQuad)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(digits = 5)



job_id_string <- Sys.getenv("SLURM_ARRAY_JOB_ID")
array_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
cat("The job ID string is", job_id_string, ".\n\n")
cat("The array ID string is", array_id_string, ".\n\n")


out_filename_jmbayes2 <- paste0("Simulations/Results/SIM_jmbayes2_", as.numeric(job_id_string), "_", str_pad(as.numeric(array_id_string), 3, pad = "0"), ".txt")
out_filename_jmwiv <- paste0("Simulations/Results/SIM_jmwiv_", as.numeric(job_id_string), "_", str_pad(as.numeric(array_id_string), 3, pad = "0"), ".txt")
cat("The JMbayes2 output will be saved file in", here(out_filename_jmbayes2), ".\n\n")
cat("The JM-WIV output will be saved file in", here(out_filename_jmwiv), ".\n\n")


library(rstanjmwiv)

cat("Simulate data.\n\n")





ranef_corr <- matrix(c(1, 0.129, 0.502, 0.019,
                       0.129, 1, -0.006, 0.243,
                       0.502, -0.006, 1, 0.491,
                       0.019, 0.243, 0.491, 1),
                     byrow = F,
                     ncol = 4)  ###check that it is symmetric
sd_vec <- c(0.81, 0.44, 0.52, 0.16) * diag(4)
ranef_mat <- sd_vec %*% ranef_corr %*% t(sd_vec)


e_beta_input <- c(0.93, -2.3)
a_beta_input <- c(-2.24, 1.08, 0.55, -0.12)
y1_mu_fixed_input <- c(2.19, -0.04)
y2_mu_fixed_input <- c(1.04, 0.35)

input_vec <- c(e_beta_input, 
               a_beta_input[c(1,3)],
               y1_mu_fixed_input,
               y2_mu_fixed_input)


simdata_jmwiv <- simulate_jmwiv(seed  = as.numeric(job_id_string) + as.numeric(array_id_string),
                                nsub = 1000,
                                assoc = 'LP',
                                lambda = 5,
                                g_shape = 0,
                                e_beta = e_beta_input,
                                a_beta = a_beta_input,
                                distObs  = 1,
                                max_time = 10,
                                ranef_covmat = ranef_mat,
                                y1mu_fixed = y1_mu_fixed_input,
                                y1sigma_fixed = c(-1.36, -0.02),
                                y2mu_fixed = y2_mu_fixed_input,
                                y2sigma_fixed = c(0.16, 0.01))



head(simdata_jmwiv$dataLong)
head(simdata_jmwiv$dataEvent)


#################### JMbayes2 #################################################
cat("Start JMbayes2.\n\n")


library(JMbayes2)

fm1 <- lme(fixed = y1 ~ visit_times, random = ~ 1 | ID, data = simdata_jmwiv$dataLong)
fm2 <- lme(fixed = y2 ~ visit_times, random = ~ 1 | ID, data = simdata_jmwiv$dataLong)
Mixed <- list(fm1, fm2)

fCox1 <- coxph(Surv(e_time, e_status) ~ binary_cov + norm_cov, simdata_jmwiv$dataEvent)
JMbayes2_fit <- jm(fCox1, Mixed, time_var = "visit_times",
                        n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)



summary(JMbayes2_fit)

JMbayesfit_mcmc <- cbind(
JMbayes2_fit$mcmc %>% 
  .$gammas %>% .[[1]],

JMbayes2_fit$mcmc %>% 
  .$alphas %>% .[[1]],

JMbayes2_fit$mcmc %>% 
  .$betas1 %>% .[[1]],

JMbayes2_fit$mcmc %>% 
  .$betas2 %>% .[[1]],

JMbayes2_fit$mcmc %>% 
  .$sigmas %>% .[[1]] %>% log(.),

JMbayes2_fit$mcmc %>% 
  .$D %>% .[[1]] %>% .[, -2] %>% sqrt()

)

colnames(JMbayesfit_mcmc)[1:4] <- paste0("e.", colnames(JMbayesfit_mcmc)[1:4])
colnames(JMbayesfit_mcmc)[5:6] <- paste0("y1.", colnames(JMbayesfit_mcmc)[5:6])
colnames(JMbayesfit_mcmc)[7:8] <- paste0("y2.", colnames(JMbayesfit_mcmc)[7:8])
colnames(JMbayesfit_mcmc)[9:10] <- paste0("logsigma.", colnames(JMbayesfit_mcmc)[9:10])
colnames(JMbayesfit_mcmc)[11:12] <- c("b_sd1", "b_sd3")


JMbayesfit_mcmc <-  data.frame(JMbayesfit_mcmc)

head(JMbayesfit_mcmc)

write.table(round(JMbayesfit_mcmc, 5), file = here(out_filename_jmbayes2), sep = "\t", row.names = FALSE, col.names = TRUE)


cat(c(unname(round(colMeans(JMbayesfit_mcmc), 5)), "JMbayes2"),"\n", file = here("Simulations/SIM_summary_JMbayes2.txt"), append = TRUE)

CI_JMbayes2 <- rbind(
  summary(JMbayes2_fit)$Survival[, c('2.5%', '97.5%')],
  summary(JMbayes2_fit)$Outcome1[1:2, c('2.5%', '97.5%')],
  summary(JMbayes2_fit)$Outcome2[1:2, c('2.5%', '97.5%')],
  summary(JMbayes2_fit)$Outcome1[3, c('2.5%', '97.5%')],
  summary(JMbayes2_fit)$Outcome2[3, c('2.5%', '97.5%')],
  JMbayes2_fit$mcmc %>%.$D %>% .[[1]] %>% .[, -2] %>% sqrt() %>% apply(., 2, quantile, probs = c(0.025,0.975)) %>% t()
)

input_vec_JMbayes2 <- c(input_vec, -1.36, 0.18, 0.81, 0.52)
coverage_JMbayes2 <- as.numeric(input_vec_JMbayes2 > CI_JMbayes2[,1] & input_vec_JMbayes2 < CI_JMbayes2[,2]) 
cat(c(coverage_JMbayes2, "JMbayes2"),"\n", file = here("Simulations/SIM_COVERAGE_JMbayes2.txt"), append = TRUE)




################## JM-WIV ########################################

cat("Start JM-WIV.\n\n")



fL1 <- list(formula(y1 ~ visit_times + (1|ID)),  
            formula(sigma ~ visit_times + (1|ID)))
fL2 <- list(formula(y2 ~ visit_times + (1|ID)),
            formula(sigma ~ visit_times + (1|ID)))
fEvent <- survival::Surv(e_time, e_status) ~ binary_cov + norm_cov




prior_input <- list(
  #y1mu_prior_mean = NULL,
  # y2mu_prior_mean = NULL,
  # y1sigma_prior_mean= NULL,
  # y2sigma_prior_mean= NULL,
  # e_prior_mean = NULL,
  # a_prior_mean = rep(0, 4),
  # ymu_prior_mean_for_intercept = NULL,
  # ysigma_prior_mean_for_intercept = NULL,
  # e_prior_mean_for_aux = NULL,  
  # y1mu_prior_scale = 5,
  # y2mu_prior_scale = 5,
  # y1sigma_prior_scale = 5,
  # y2sigma_prior_scale = 5,
  # e_prior_scale = 5,
  a_prior_scale = rep(1, 4),
  # ymu_prior_scale_for_intercept = rep(0.5, 2),
  # ysigma_prior_scale_for_intercept = rep(0.5, 2),
  e_prior_scale_for_aux = rep(3, 6),     
  b_prior_scale = rep(3, 4),
  b_prior_df = rep(3, 4),
  b_prior_regularization = 3
)


simdata_jmwiv[[1]]$ID <- as.numeric(simdata_jmwiv[[1]]$ID)
simdata_jmwiv[[2]]$ID <- as.numeric(simdata_jmwiv[[2]]$ID)





simData_jmwiv_LPassoc <- rstanjmwiv::jmwiv_stan(formulaLong1 = fL1,
	formulaLong2 = fL2,
           dataLong = simdata_jmwiv$dataLong,
	id_var = "ID",
           time_varLong = "visit_times",
           formulaEvent = fEvent,
           dataEvent = simdata_jmwiv$dataEvent,
           a_K = 4,
           assoc = "LP",
           basehaz = "bs",
           basehaz_aux = list(df = 6,  
                              knots = NULL, 
                              degree = 3),
           qnodes = 15L,
           prior_list = prior_input,
  	cores = 2,            
  	chains = 2,             
  	warmup = 3000,          
  	iter = 6000,            
  	#refresh = 0.2,           
  	seed = as.numeric(job_id_string) + as.numeric(array_id_string), 
  	control = list(max_treedepth = 15) 
)









print(simData_jmwiv_LPassoc,  
      pars = c("y1mu_Intercept","y1sigma_Intercept","y2mu_Intercept","y2sigma_Intercept",
               "y1mu_beta","y1sigma_beta","y2mu_beta","y2sigma_beta"),
      probs = c(0.025,0.975))
print(simData_jmwiv_LPassoc,  pars = "a_beta", probs = c(0.025,0.975))
print(simData_jmwiv_LPassoc,  pars = c("e_Intercept", "e_beta"), probs = c(0.025,0.975))
print(simData_jmwiv_LPassoc,  pars = "b_sd", probs = c(0.025,0.975))


JMWIV_mcmc <- cbind(
  rstan::extract(simData_jmwiv_LPassoc, pars = c("e_beta"))$e_beta,
  rstan::extract(simData_jmwiv_LPassoc, pars = c("a_beta"))$a_beta[,c(1,3)],
  rstan::extract(simData_jmwiv_LPassoc, pars = c("y1mu_Intercept"))[[1]],
  rstan::extract(simData_jmwiv_LPassoc, pars = c("y1mu_beta"))[[1]],
  rstan::extract(simData_jmwiv_LPassoc, pars = c("y2mu_Intercept"))[[1]],
  rstan::extract(simData_jmwiv_LPassoc, pars = c("y2mu_beta"))[[1]],
  rstan::extract(simData_jmwiv_LPassoc, pars = c("b_sd"))[[1]],
  rstan::extract(simData_jmwiv_LPassoc, pars = c("a_beta"))$a_beta[,c(2,4)],
  rstan::extract(simData_jmwiv_LPassoc, pars = c("y1sigma_Intercept"))[[1]],
  rstan::extract(simData_jmwiv_LPassoc, pars = c("y1sigma_beta"))[[1]],
  rstan::extract(simData_jmwiv_LPassoc, pars = c("y2sigma_Intercept"))[[1]],
  rstan::extract(simData_jmwiv_LPassoc, pars = c("y2sigma_beta"))[[1]]
)

JMWIV_mcmc <- data.frame(JMWIV_mcmc)

write.table(round(JMWIV_mcmc, 5), file = here(out_filename_jmwiv), sep = "\t",
            row.names = FALSE, col.names = TRUE)

cat(c(unname(round(colMeans(JMWIV_mcmc), 5)), "JM-WIV"),"\n", file = here("Simulations/SIM_summary_JMWIV.txt"), append = TRUE)


CI_JMWIV <- summary(simData_jmwiv_LPassoc)$summary[c("e_beta[1]", "e_beta[2]", 
                                     "a_beta[1]", "a_beta[3]",
                                     "y1mu_Intercept","y1mu_beta[1]", 
                                     "y2mu_Intercept","y2mu_beta[1]",
			  	    "b_sd[1]", "b_sd[2]", "b_sd[3]", "b_sd[4]",
			  	    "a_beta[2]", "a_beta[4]",
			  	    "y1sigma_Intercept","y1sigma_beta[1]", 
                                     "y2sigma_Intercept","y2sigma_beta[1]"), c('2.5%', '97.5%')]



input_vec_JMWIV <- c(input_vec, c(0.81, 0.44, 0.52, 0.16), a_beta_input[c(2,4)], c(-1.36, -0.02), c(0.16, 0.01))
coverage_JMWIV <- as.numeric(input_vec_JMWIV > CI_JMWIV[,1] & input_vec_JMWIV < CI_JMWIV[,2]) 
cat(c(coverage_JMWIV, "JMWIV"),"\n", file = here("Simulations/SIM_COVERAGE_JMWIV.txt"), append = TRUE)

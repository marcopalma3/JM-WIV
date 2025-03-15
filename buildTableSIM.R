library(pacman)
p_load(dplyr, tidyr)

data_JMbayes2 <- read.table("Simulations/SIM_summary_JMbayes2.txt", header = TRUE) 
data_JMWIV <- read.table("Simulations/SIM_summary_JMWIV.txt", header = TRUE) 

summary(data_JMbayes2)
summary(data_JMWIV)

names(data_JMbayes2)[9:12] <- names(data_JMWIV)[c(15, 17, 9, 11)]
data_JMbayes2_withNA <- data.frame(data_JMbayes2[, 1:12], matrix(NA, ncol = 6), data_JMbayes2$model) 
names(data_JMbayes2_withNA)[13:18] <-  names(data_JMWIV)[c(10, 12:14, 16, 18)]


coverage_JMbayes2 <- read.table("Simulations/SIM_COVERAGE_JMbayes2.txt", header = TRUE) %>% tail(200)
coverage_JMWIV <- read.table("Simulations/SIM_COVERAGE_JMWIV.txt", header = TRUE) %>% tail(200)

colMeans(coverage_JMbayes2[, -19])
colMeans(coverage_JMWIV[, -19])


names(coverage_JMbayes2)[9:12] <- names(coverage_JMWIV)[c(15, 17, 9, 11)]
coverage_JMbayes2_withNA <- data.frame(coverage_JMbayes2[, 1:12], matrix(NA, ncol = 6), coverage_JMbayes2$model) 
names(coverage_JMbayes2_withNA)[13:18] <-  names(coverage_JMWIV)[c(10, 12:14, 16, 18)]





####################
e_beta_input <- c(0.93, -2.3)
a_beta_input <- c(-2.24, 1.08, 0.55, -0.12)
y1_mu_fixed_input <- c(2.19, -0.04)
y2_mu_fixed_input <- c(1.04, 0.35)
####################
sd_vec_input <- c(0.81, 0.44, 0.52, 0.16)
y1sigma_fixed_input <-  c(-1.36, -0.02)
y2sigma_fixed_input = c(0.16, 0.01)


input_vec <- c(y1_mu_fixed_input, sd_vec_input[1],
	    y1sigma_fixed_input, sd_vec_input[2], 
	    y2_mu_fixed_input, sd_vec_input[3], 
	    y2sigma_fixed_input, sd_vec_input[4], 
		e_beta_input, 
               	a_beta_input)


vec_names <- c("y1..Intercept.",  "y1.visit_times", "b_sd1",
		"y1sigma..Intercept.", "y1sigma.visit_times", "b_sd2",
		"y2..Intercept.", "y2.visit_times", "b_sd3",
		"y2sigma..Intercept.", "y2sigma.visit_times", "b_sd4",
		"e.binary_cov", "e.norm_cov", 
		"e.value.y1.",  "e.value.y1sigma", "e.value.y2.", "e.value.y2sigma" 
)



options(knitr.kable.NA = '-')
data.frame("input" = input_vec,
	"Mean_JMWIV" = colMeans(data_JMWIV[, vec_names]),
	"SEMean_JMWIV" = apply(data_JMWIV[, vec_names], 2, sd),
	"Coverage_JMWIV" = colMeans(coverage_JMWIV[, vec_names]),
	"Mean_JMbayes2" = colMeans(data_JMbayes2_withNA[, vec_names]),
	"SEMean_JMbayes2" = apply(data_JMbayes2_withNA[, vec_names], 2, sd),
	"Coverage_JMbayes2" = colMeans(coverage_JMbayes2_withNA[, vec_names])
) %>% round(3) %>% knitr::kable(format = "latex")






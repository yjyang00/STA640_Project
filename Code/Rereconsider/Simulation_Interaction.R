library(lme4)
library(tidyverse)
library(foreach)
library(doParallel)
library(magrittr)

set.seed(19991109, kind = "L'Ecuyer-CMRG")

# Initialize parallel backend, adjust as needed
cl <- makeCluster(8)

# Register the parallel backend
registerDoParallel(cl)


panel_bias_sim = function(S = 100, n, t, t_treat, delta, gamma, rho, phi, confound_treatment){
  
  result = foreach (i = 1:S, .combine = 'rbind', .errorhandling='remove') %dopar% {
    library(lme4)
    library(tidyverse)
    
    X <- runif(n, 0, 10)
    X_it = c()
    for (j in 1:n){
      X_it = c(X_it, rnorm(n = t, mean = X[j], sd = 1))
    }
    
    if (confound_treatment == "Small"){
      A <- rnorm(n)
      p = ifelse(A >= 0, 0.55, 0.45)
      D <- rbinom(n, 1, p)
    }
    
    if (confound_treatment == "Strong"){
      A <- rnorm(n)
      p = ifelse(A >= 0, 0.7, 0.3)
      D <- rbinom(n, 1, p)
    }
    
    if (confound_treatment == "Very Strong"){
      A <- runif(n)
      p = A
      D <- rbinom(n, 1, A)
    }
    
    dat = tidyr::expand_grid(data.frame(id = 1:n, A = A, D = D), Time = 1:t) %>% 
      mutate(X_it = X_it) %>% 
      mutate(time = factor(Time), id = factor(id)) %>% 
      mutate(Time_to_treat = ifelse(Time < t_treat, 0, 1)) %>% 
      mutate(D_it = D * Time_to_treat) %>% 
      mutate(epsilon = rnorm(n*t, mean = 0, sd = 1)) %>% 
      mutate(Y = (Time)^2 + delta*X_it + gamma*A + rho*D_it + phi*(A)^3*D_it + epsilon)
    
    # dat %>% ggplot(aes(x = Time, y = Y, group = factor(id), color = factor(D))) + geom_point() + geom_line() + facet_wrap(~factor(D))
    
    # DID estimator
    bar_Y_1_t2 = dat %>% filter(Time == t_treat, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
    bar_Y_1_t1 = dat %>% filter(Time == t_treat-1, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
    
    bar_Y_0_t2 = dat %>% filter(Time == t_treat, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
    bar_Y_0_t1 = dat %>% filter(Time == t_treat-1, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
    
    DID_est = (bar_Y_1_t2 - bar_Y_1_t1) - (bar_Y_0_t2 - bar_Y_0_t1)
    DID_bias = DID_est - rho
    
    # OLS
    mod0 = lm(Y ~ X_it + time + D_it, dat)
    # summary(mod0)
    OLS_est = tail(mod0$coefficients, 1)
    OLS_bias = (OLS_est - rho) %>% as.numeric()
    OLS_CI = tail(confint(mod0),1)[1] < rho & rho < tail(confint(mod0),1)[2]
    
    # fixed-effects model
    mod1 = lm(Y ~ id + X_it + time + D_it - 1, dat)
    # summary(mod1)
    FE_est = tail(mod1$coefficients, 1)
    FE_bias = (FE_est - rho) %>% as.numeric()
    FE_CI = tail(confint(mod1),1)[1] < rho & rho < tail(confint(mod1),1)[2]
    
    # random-effects model (with random intercept)
    mod2 = lme4::lmer(Y ~ X_it + time + D_it + (1 | id) - 1, dat)
    summ_lmer = summary(mod2)
    # summ_lmer
    RE_est = tail(summ_lmer$coefficients, 1)[1]
    RE_bias = (RE_est - rho) %>% as.numeric()
    RE_ci =  tail(confint(mod2, method="Wald"),1)
    RE_CI = RE_ci[1] < rho & rho < RE_ci[2]
    
    
    res_est = c(DID_est, OLS_est, FE_est, RE_est)
    res_bias = c(DID_bias, OLS_bias, FE_bias, RE_bias)
    CI_capture = c(OLS_CI, FE_CI, RE_CI)
    # c(res_est, res_bias, CI_capture)
    return(c(res_est, res_bias, CI_capture))
    
  }
  
  average_bias = colMeans(result)
  paste(average_bias, collapse=" ")
}



n <- 50 
t <- 10 
t_treat <- 5
delta <- 10
gamma <- 1
rho <- 1
phi = (0:5)*2
confound_treatment = c("Small","Strong","Very Strong")


# Define a function to be applied in parallel
bias_mutate <- function(df) {
  df %>%  mutate(result = panel_bias_sim(S = 500, n, t, t_treat, delta, gamma, rho, phi, confound_treatment))
}

bias_result_interaction <- tidyr::expand_grid(n, t, t_treat, delta, gamma, rho, phi, confound_treatment) %>% 
  group_by(n, t, t_treat, delta, gamma, rho, phi, confound_treatment) %>% 
  do(bias_mutate(.)) %>% separate(result, c("DID_est", "OLS_est", "FE_est", "RE_est",
                                            "DID_bias", "OLS_bias", "FE_bias", "RE_bias",
                                            "OLS_CI", "FE_CI", "RE_CI"), " ", convert = TRUE)


# save(bias_result_interaction, file = "bias_result_interaction.RData")

# no interaction

bias_result_interaction %>% pivot_longer(cols = c("OLS_bias","FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(phi == 0, gamma == 8)%>% 
  ggplot(aes(x = rho, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)

bias_result_interaction %>% pivot_longer(cols = c("OLS_bias","FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(phi == 0, rho == 8)%>% 
  ggplot(aes(x = gamma, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)

bias_result_interaction %>% pivot_longer(cols = c("FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(phi == 0, gamma == 8)%>% 
  ggplot(aes(x = rho, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)

bias_result_interaction %>% pivot_longer(cols = c("FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(phi == 0, rho == 8)%>% 
  ggplot(aes(x = gamma, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)


# some interaction

bias_result_interaction %>% pivot_longer(cols = c("OLS_bias","FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(phi == 2, gamma == 8)%>% 
  ggplot(aes(x = rho, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)

bias_result_interaction %>% pivot_longer(cols = c("OLS_bias","FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(phi == 2, rho == 8)%>% 
  ggplot(aes(x = gamma, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)

bias_result_interaction %>% pivot_longer(cols = c("FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(phi == 2, gamma == 8)%>% 
  ggplot(aes(x = rho, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)

bias_result_interaction %>% pivot_longer(cols = c("FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(phi == 2, rho == 8)%>% 
  ggplot(aes(x = gamma, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)




# interaction 

bias_result_interaction %>% pivot_longer(cols = c("OLS_bias","FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(rho == 2, gamma == 2)%>% 
  ggplot(aes(x = phi, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)

bias_result_interaction %>% pivot_longer(cols = c("FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(rho == 2, gamma == 2)%>% 
  ggplot(aes(x = phi, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)

bias_result_interaction %>% pivot_longer(cols = c("OLS_CI", "FE_CI", "RE_CI"), names_to = "Bias_Type", values_to = "Bias") %>%
  filter(rho == 2, gamma == 2)%>% 
  ggplot(aes(x = phi, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)


bias_result_interaction %>% pivot_longer(cols = c("OLS_bias","FE_bias", "RE_bias"), names_to = "Bias_Type", values_to = "Bias") %>%
  ggplot(aes(x = phi, y = Bias, color = Bias_Type)) + geom_point() + geom_line() + facet_wrap(~confound_treatment)

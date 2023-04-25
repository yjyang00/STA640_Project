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
      p = ifelse(A >= 1, 0.55, 0.45)
      D <- rbinom(n, 1, p)
    }
    
    if (confound_treatment == "Strong"){
      A <- rnorm(n)
      p = ifelse(A >= 1, 0.7, 0.3)
      D <- rbinom(n, 1, p)
    }
    
    if (confound_treatment == "Uniform"){
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
      mutate(Y = (Time)^2 + delta*X_it + gamma*A + rho*D_it + epsilon)
    
    # dat %>% ggplot(aes(x = Time, y = Y, group = factor(id), color = factor(D))) + geom_point() + geom_line() + facet_wrap(~factor(D))
    
    # DID estimator
    bar_Y_1_t2 = dat %>% filter(Time == t_treat, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
    bar_Y_1_t1 = dat %>% filter(Time == t_treat-1, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
    
    bar_Y_0_t2 = dat %>% filter(Time == t_treat, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
    bar_Y_0_t1 = dat %>% filter(Time == t_treat-1, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
    
    DID_est = (bar_Y_1_t2 - bar_Y_1_t1) - (bar_Y_0_t2 - bar_Y_0_t1)
    DID_bias = abs(DID_est - rho)
    
    # OLS
    mod0 = lm(Y ~ X_it + time + D_it, dat)
    # summary(mod0)
    # tail(confint(mod0),1)
    OLS_est = tail(mod0$coefficients, 1)
    OLS_bias = abs(OLS_est - rho) %>% as.numeric()
    
    # fixed-effects model
    mod1 = lm(Y ~ id + X_it + time + D_it - 1, dat)
    # summary(mod1)
    # tail(confint(mod1),1)
    FE_est = tail(mod1$coefficients, 1)
    FE_bias = abs(FE_est - rho) %>% as.numeric()
    
    # random-effects model (with random intercept)
    mod2 = lme4::lmer(Y ~ X_it + time + D_it + (1 | id) - 1, dat)
    summ_lmer = summary(mod2)
    # summ_lmer
    # tail(confint(mod2),1)
    RE_est = tail(summ_lmer$coefficients,1)[1]
    RE_bias = abs(RE_est - rho) %>% as.numeric()
    
    
    res_est = c(DID_est, OLS_est, FE_est, RE_est)
    res_bias = c(DID_bias, OLS_bias, FE_bias, RE_bias)
    # c(res_est, res_bias)
    return(c(res_est, res_bias))
    
  }
  
  average_bias = colMeans(result)
  paste(average_bias, collapse=" ")
}



n <- 100 
t <- 10 
t_treat <- 5
delta <- 5
gamma <- (1:10)*2
rho <- (1:10)*2
phi = 0
confound_treatment = c("Small","Strong","Uniform")


# Define a function to be applied in parallel
bias_mutate <- function(df) {
  df %>%  mutate(result = panel_bias_sim(S = 100, n, t, t_treat, delta, gamma, rho, phi, confound_treatment))
}

bias_result_no_interaction <- tidyr::expand_grid(n, t, t_treat, delta, gamma, rho, phi, confound_treatment) %>% 
  group_by(n, t, t_treat, delta, gamma, rho, phi, confound_treatment) %>% 
  do(bias_mutate(.)) %>% separate(result, c("DID_est", "OLS_est", "FE_est", "RE_est",
                                            "DID_bias", "OLS_bias", "FE_bias", "RE_bias"), " ", convert = TRUE)


# save(bias_result_no_interaction, file = "bias_result_no_interaction.RData")
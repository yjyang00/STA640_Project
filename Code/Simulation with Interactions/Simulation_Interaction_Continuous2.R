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


panel_bias_sim = function(S = 100, n, t, t_treat, delta, gamma, rho, phi, beta, beta_a){
  
  result = foreach (i = 1:S, .combine = 'rbind', .errorhandling='remove') %dopar% {
    library(lme4)
    library(tidyverse)
    
    X <- runif(n, 0, 10)
    X_it = c()
    for (j in 1:n){
      X_it = c(X_it, rnorm(n = t, mean = X[j], sd = 1))
    }
    
    A <- rnorm(n, 0, 3)
    alpha_a <- rnorm(n, 0, 1)
    
    dat = tidyr::expand_grid(data.frame(id = 1:n, A, alpha_a), Time = 1:t) %>% 
      mutate(epsilon_a = rnorm(n*t, mean = 0, sd = 1)) %>% 
      mutate(D_it = alpha_a + beta_a*A + epsilon_a) %>% 
      mutate(X_it = X_it) %>% 
      mutate(time = factor(Time), id = factor(id)) %>% 
      mutate(epsilon_y = rnorm(n*t, mean = 0, sd = 1)) %>% 
      mutate(Y = (Time)^2 + gamma*A + 0*X_it + rho*D_it + phi*A*D_it + beta*A*Time + epsilon_y)
    
    
    # OLS
    mod0 = lm(Y ~  time + D_it, dat)
    # summary(mod0)
    OLS_est = tail(mod0$coefficients, 1)
    OLS_bias = (OLS_est - rho) %>% as.numeric()
    OLS_CI = tail(confint(mod0),1)[1] < rho & rho < tail(confint(mod0),1)[2]
    
    # fixed-effects model
    mod1 = lm(Y ~ id + time + D_it - 1, dat)
    # summary(mod1)
    FE_est = tail(mod1$coefficients, 1)
    FE_bias = (FE_est - rho) %>% as.numeric()
    FE_CI = tail(confint(mod1),1)[1] < rho & rho < tail(confint(mod1),1)[2]
    
    # random-effects model (with random intercept)
    mod2 = lme4::lmer(Y ~  time + D_it + (1 | id) - 1, dat)
    summ_lmer = summary(mod2)
    # summ_lmer
    RE_est = tail(summ_lmer$coefficients, 1)[1]
    RE_bias = (RE_est - rho) %>% as.numeric()
    RE_ci =  tail(confint(mod2, method="Wald"),1)
    RE_CI = RE_ci[1] < rho & rho < RE_ci[2]
    
    
    res_est = c(OLS_est, FE_est, RE_est)
    res_bias = c(OLS_bias, FE_bias, RE_bias)
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
delta <- 1
gamma <- c(0,1,3,5,7)
rho <- 1
phi =  0
beta = 0
beta_a = c(0,1,3,5,7)

# Define a function to be applied in parallel
bias_mutate <- function(df) {
  df %>%  mutate(result = panel_bias_sim(S = 10000, n, t, t_treat, delta, gamma, rho, phi, beta, beta_a))
}

bias_result_continuous <- tidyr::expand_grid(n, t, t_treat, delta, gamma, rho, phi, beta, beta_a) %>% 
  group_by(n, t, t_treat, delta, gamma, rho, phi, beta, beta_a) %>% 
  do(bias_mutate(.)) %>% separate(result, c("OLS_est", "FE_est", "RE_est",
                                            "OLS_bias", "FE_bias", "RE_bias",
                                            "OLS_CI", "FE_CI", "RE_CI"), " ", convert = TRUE)


# save(bias_result_interaction_phi, file = "result_data/bias_result_interaction_phi.RData")


bias_result_continuous %>% pivot_longer(cols = c("OLS_bias","FE_bias", "RE_bias"), names_to = "Model", values_to = "Bias") %>% 
  filter(beta_a == 3) %>% 
  ggplot(aes(x = gamma, y = Bias, color = Model)) + geom_point() + geom_line() 



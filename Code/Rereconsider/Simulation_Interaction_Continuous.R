library(tidyverse)
library(foreach)
library(doParallel)

set.seed(19991109, kind = "L'Ecuyer-CMRG")

# Initialize parallel backend
cl <- makeCluster(8)

# Register the parallel backend
registerDoParallel(cl)


S = 1000
result = foreach (i = 1:S, .combine = 'rbind', .errorhandling='remove') %dopar% {
  library(tidyverse)
  library(lme4)
  n <- 50 
  t <- 10
  t_treat <- 5
  
  delta <- 1
  gamma <- 1
  rho <- 1
  phi <- 0
  beta <- 0
  
  X <- runif(n, 0, 10)
  
  X_it = c()
  for (j in 1:n){
    X_it = c(X_it, rnorm(n = t, mean = X[j], sd = 1))
  }
  
  A <- rnorm(n, 0, 1)
  alpha_a <- rnorm(n, 0, 1)
  beta_a <- 1
  
  dat = tidyr::expand_grid(data.frame(id = 1:n, A, alpha_a), Time = 1:t) %>% 
    mutate(epsilon_a = rnorm(n*t, mean = 0, sd = 1)) %>% 
    mutate(D_it = alpha_a + beta_a*A + epsilon_a) %>% 
    mutate(X_it = X_it) %>% 
    mutate(time = factor(Time), id = factor(id)) %>% 
    mutate(epsilon_y = rnorm(n*t, mean = 0, sd = 1)) %>% 
    mutate(Y = (Time)^2 + gamma*A + delta*X_it + rho*D_it + phi*A*D_it + beta*A*Time + epsilon_y)
    
  # dat %>% ggplot(aes(x = Time, y = Y, group = factor(id), color = factor(D))) + geom_line() + facet_wrap(~factor(D))
  
  # OLS
  mod0 = lm(Y ~ X_it + time + D_it, dat)
  # summary(mod0)
  # tail(confint(mod0),1)
  OLS_est = tail(mod0$coefficients, 1)
  OLS_bias = (OLS_est - rho) %>% as.numeric()
  OLS_CI = tail(confint(mod0),1)[1] < rho & rho < tail(confint(mod0),1)[2]

  # fixed-effects model
  mod1 = lm(Y ~ id + X_it + time + D_it - 1, dat)
  # summary(mod1)
  # tail(confint(mod1),1)
  FE_est = tail(mod1$coefficients, 1)
  FE_bias = (FE_est - rho) %>% as.numeric()
  FE_CI = tail(confint(mod1),1)[1] < rho & rho < tail(confint(mod1),1)[2]
  
  # random-effects model (with random intercept)
  mod2 = lme4::lmer(Y ~ X_it + time + D_it + (1 | id) - 1, dat)
  summ_lmer = summary(mod2)
  # summ_lmer
  # tail(confint(mod2),1)
  RE_est = tail(summ_lmer$coefficients,1)[1]
  RE_bias = (RE_est - rho) %>% as.numeric()
  RE_ci =  tail(confint(mod2, method="Wald"),1)
  RE_CI = RE_ci[1] < rho & rho < RE_ci[2]
  
  
  res_est = c(DID_est, OLS_est, FE_est, RE_est)
  res_bias = c(DID_bias, OLS_bias, FE_bias, RE_bias)
  CI_capture = c(OLS_CI, FE_CI, RE_CI)
  # c(res_est, res_bias, CI_capture)
  return(c(res_est, res_bias, CI_capture))
}

colnames(result) = c("DID_est", "OLS_est", "FE_est", "RE_est", 
                     "DID_bias", "OLS_bias", "FE_bias", "RE_bias",
                     "OLS_CI", "FE_CI", "RE_CI")
colMeans(result)

result_10_5_gamma1_beta0 = result

save(result_10_5_gamma1_beta0, file = "result_data/result_10_5_gamma1_beta0.RData")


# result %>% as_tibble() %>%
#   pivot_longer(cols = c("OLS_bias","FE_bias", "RE_bias"), names_to = "Type", values_to = "Value")%>% 
#   mutate(Type = factor(Type))%>% 
#   mutate(Type = factor(Type, levels = c("OLS_bias","FE_bias", "RE_bias")))%>% 
#   ggplot() + geom_boxplot(aes(x = Type, y = Value))
# 
# result %>% as_tibble() %>%
#   pivot_longer(cols = c("FE_bias", "RE_bias"), names_to = "Type", values_to = "Value")%>%
#   mutate(Type = factor(Type))%>%
#   mutate(Type = factor(Type, levels = c("FE_bias", "RE_bias")))%>%
#   ggplot() + geom_boxplot(aes(x = Type, y = Value))

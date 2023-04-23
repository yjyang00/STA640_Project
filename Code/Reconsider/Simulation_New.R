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
  n <- 100 
  t <- 10
  t_treat <- 5
  
  delta <- 5
  gamma <- 3
  rho <- 8
  
  A <- rnorm(n)
  X <- rnorm(n)
  
  p = ifelse(A >= 1, 0.75, 0.45)
  D <- rbinom(n, 1, p)
  
  dat = tidyr::expand_grid(data.frame(id = 1:n, A = A, X = X, D = D), Time = 1:t) %>% 
    mutate(id = factor(id)) %>% 
    mutate(epsilon = rnorm(n*t, mean = 0, sd = 1)) %>% 
    mutate(Y  = ifelse(Time < t_treat, 
                       (Time)^2 + delta*X + gamma*A + epsilon,
                       (Time)^2 + delta*X + gamma*A + rho*D + epsilon))
  
  # dat %>% ggplot(aes(x = Time, y = Y, group = factor(id), color = factor(D))) + geom_point() + geom_line() + facet_wrap(~factor(D))

  
  # DID estimator
  bar_Y_1_t2 = dat %>% filter(Time == t_treat, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
  bar_Y_1_t1 = dat %>% filter(Time == t_treat-1, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
  
  bar_Y_0_t2 = dat %>% filter(Time == t_treat, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
  bar_Y_0_t1 = dat %>% filter(Time == t_treat-1, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
  
  DID_est = (bar_Y_1_t2 - bar_Y_1_t1) - (bar_Y_0_t2 - bar_Y_0_t1)
  DID_bias = DID_est - rho
  
  
  # create D_it
  dat = dat %>% mutate(time = factor(Time)) %>% 
    mutate(Time_to_treat = ifelse(Time < t_treat, 0, 1)) %>% 
    mutate(D_it = D * Time_to_treat)
  
  # OLS
  mod0 = lm(Y ~ X + time + D_it, dat)
  # summary(mod0)
  # tail(confint(mod0),1)
  OLS_est = tail(mod0$coefficients, 1)
  OLS_bias = (OLS_est - rho) %>% as.numeric()

  # fixed-effects model
  mod1 = lm(Y ~ id + X + time + D_it, dat)
  # summary(mod1)
  # tail(confint(mod1),1)
  FE_est = tail(mod1$coefficients, 1)
  FE_bias = (FE_est - rho) %>% as.numeric()
  
  # random-effects model (with random intercept)
  mod2 = lme4::lmer(Y ~ X + time + D_it + (1 | id), dat)
  summ_lmer = summary(mod2)
  # summ_lmer
  # tail(confint(mod2),1)
  RE_est = tail(summ_lmer$coefficients,1)[1]
  RE_bias = (RE_est - rho) %>% as.numeric()
  
  
  res_est = c(DID_est, OLS_est, FE_est, RE_est)
  res_bias = c(DID_bias, OLS_bias, FE_bias, RE_bias)
  
  return(c(res_est, res_bias))
}

colnames(result) = c("DID_est", "OLS_est", "FE_est", "RE_est", "DID_bias", "OLS_bias", "FE_bias", "RE_bias")
colMeans(result)


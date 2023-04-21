library(foreach)
library(doParallel)

set.seed(19991109, kind = "L'Ecuyer-CMRG")

# Initialize parallel backend
cl <- makeCluster(8)

# Register the parallel backend
registerDoParallel(cl)


S = 5000
result = foreach (i = 1:S, .combine = 'rbind', .errorhandling='remove') %dopar% {
  library(tidyverse)
  library(lme4)
  n <- 100 
  t <- 2 
  delta <- 5
  gamma <- 3
  rho <- 8
  
  A <- rnorm(n)
  X <- rnorm(n)

  p = ifelse(A >= 1, 0.75, 0.45)
  D <- rbinom(n, 1, p)
  
  dat = tidyr::expand_grid(data.frame(id = 1:n, A = A, X = X, D = D), Time = 1:t) %>% 
    mutate(epsilon = rnorm(n*t, mean = 0, sd = 1)) %>% 
    mutate(Y  = ifelse(Time < 2, 
                       (Time) + delta*X + gamma*A + epsilon,
                       (Time)^2 + delta*X + gamma*A + rho*D + epsilon))
  
  
  # DID estimator
  bar_Y_1_t2 = dat %>% filter(Time == 2, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
  bar_Y_1_t1 = dat %>% filter(Time == 1, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
  
  bar_Y_0_t2 = dat %>% filter(Time == 2, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
  bar_Y_0_t1 = dat %>% filter(Time == 1, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
  
  DID_est = (bar_Y_1_t2 - bar_Y_1_t1) - (bar_Y_0_t2 - bar_Y_0_t1)
  DID_bias = DID_est - rho


  # create D_it
  dat = dat %>% mutate(time = factor(Time)) %>% 
    mutate(Time_to_treat = ifelse(Time < 2, 0, 1)) %>% 
    mutate(D_it = D * Time_to_treat)
  
  # OLS
  mod0 = lm(Y ~ X + time + D_it, dat)
  # summary(mod0)
  OLS_bias = (mod0$coefficients[4] - rho) %>% as.numeric()
  
  
  # FE: allow different intercepts for each id
  mod1 = lme4::lmer(Y ~ time + X + D_it + (1 | id), dat)
  summ_lmer = summary(mod1)
  # summ_lmer
  FE_bias = (summ_lmer$coefficients[4,1] - rho) %>% as.numeric()
  
  res = c(DID_bias, OLS_bias, FE_bias)
  
  return(res)
}

colnames(result) = c("DID_bias", "OLS_bias", "FE_bias")
colMeans(result)

# save(result, file = "result.RData")



### simulation 2

S = 5000
result2 = foreach (i = 1:S, .combine = 'rbind', .errorhandling='remove') %dopar% {
  library(tidyverse)
  library(lme4)
  n <- 100 
  t <- 2 
  delta <- 5
  gamma <- 3
  rho <- 8
  
  A <- rnorm(n)
  X <- rnorm(n)
  
  p = ifelse(A >= 1, 0.85, 0.35)
  D <- rbinom(n, 1, p)
  
  dat = tidyr::expand_grid(data.frame(id = 1:n, A = A, X = X, D = D), Time = 1:t) %>% 
    mutate(epsilon = rnorm(n*t, mean = 0, sd = 1)) %>% 
    mutate(Y  = ifelse(Time < 2, 
                       (Time) + delta*X + gamma*A + epsilon,
                       (Time)^2 + delta*X + gamma*A + rho*D + epsilon))
  
  
  # DID estimator
  bar_Y_1_t2 = dat %>% filter(Time == 2, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
  bar_Y_1_t1 = dat %>% filter(Time == 1, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
  
  bar_Y_0_t2 = dat %>% filter(Time == 2, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
  bar_Y_0_t1 = dat %>% filter(Time == 1, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
  
  DID_est = (bar_Y_1_t2 - bar_Y_1_t1) - (bar_Y_0_t2 - bar_Y_0_t1)
  DID_bias = DID_est - rho
  
  
  # create D_it
  dat = dat %>% mutate(time = factor(Time)) %>% 
    mutate(Time_to_treat = ifelse(Time < 2, 0, 1)) %>% 
    mutate(D_it = D * Time_to_treat)
  
  # OLS
  mod0 = lm(Y ~ X + time + D_it, dat)
  # summary(mod0)
  OLS_bias = (mod0$coefficients[4] - rho) %>% as.numeric()
  
  
  # FE: allow different intercepts for each id
  mod1 = lme4::lmer(Y ~ time + X + D_it + (1 | id), dat)
  summ_lmer = summary(mod1)
  # summ_lmer
  FE_bias = (summ_lmer$coefficients[4,1] - rho) %>% as.numeric()
  
  res = c(DID_bias, OLS_bias, FE_bias)
  
  return(res)
}

colMeans(result2)
colnames(result2) = c("DID_bias", "OLS_bias", "FE_bias")

# save(result2, file = "result2.RData")




### simulation 3

S = 5000
result3 = foreach (i = 1:S, .combine = 'rbind', .errorhandling='remove') %dopar% {
  library(tidyverse)
  library(lme4)
  n <- 100 
  t <- 2 
  delta <- 5
  gamma <- 10
  rho <- 8
  
  A <- rnorm(n)
  X <- rnorm(n)
  
  p = ifelse(A >= 1, 0.85, 0.35)
  D <- rbinom(n, 1, p)
  
  dat = tidyr::expand_grid(data.frame(id = 1:n, A = A, X = X, D = D), Time = 1:t) %>% 
    mutate(epsilon = rnorm(n*t, mean = 0, sd = 1)) %>% 
    mutate(Y  = ifelse(Time < 2, 
                       (Time) + delta*X + gamma*A + epsilon,
                       (Time)^2 + delta*X + gamma*A + rho*D + epsilon))
  
  
  # DID estimator
  bar_Y_1_t2 = dat %>% filter(Time == 2, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
  bar_Y_1_t1 = dat %>% filter(Time == 1, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
  
  bar_Y_0_t2 = dat %>% filter(Time == 2, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
  bar_Y_0_t1 = dat %>% filter(Time == 1, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
  
  DID_est = (bar_Y_1_t2 - bar_Y_1_t1) - (bar_Y_0_t2 - bar_Y_0_t1)
  DID_bias = DID_est - rho
  
  
  # create D_it
  dat = dat %>% mutate(time = factor(Time)) %>% 
    mutate(Time_to_treat = ifelse(Time < 2, 0, 1)) %>% 
    mutate(D_it = D * Time_to_treat)
  
  # OLS
  mod0 = lm(Y ~ X + time + D_it, dat)
  # summary(mod0)
  OLS_bias = (mod0$coefficients[4] - rho) %>% as.numeric()
  
  
  # FE: allow different intercepts for each id
  mod1 = lme4::lmer(Y ~ time + X + D_it + (1 | id), dat)
  summ_lmer = summary(mod1)
  # summ_lmer
  FE_bias = (summ_lmer$coefficients[4,1] - rho) %>% as.numeric()
  
  res = c(DID_bias, OLS_bias, FE_bias)
  
  return(res)
}

colMeans(result3)
colnames(result3) = c("DID_bias", "OLS_bias", "FE_bias")

# save(result3, file = "result3.RData")


load("result.RData")
load("result2.RData")
load("result3.RData")
colMeans(result)
colMeans(result2)
colMeans(result3)

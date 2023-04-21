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


bias_sim = function(S = 100, n, t, delta, gamma, rho, confound_treatment){
  
  result = foreach (i = 1:S, .combine = 'rbind', .errorhandling='remove') %dopar% {
    library(lme4)
    library(tidyverse)
    
    A <- rnorm(n)
    X <- rnorm(n)
    
    if (confound_treatment == 1){
      p = ifelse(A >= 1, 0.55, 0.45)
      D <- rbinom(n, 1, p)
    }
    
    if (confound_treatment == 2){
      p = ifelse(A >= 1, 0.65, 0.35)
      D <- rbinom(n, 1, p)
    }
    
    if (confound_treatment == 3){
      p = ifelse(A >= 1, 0.75, 0.25)
      D <- rbinom(n, 1, p)
    }
    
    
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
  
  # colnames(result) = c("DID_bias", "OLS_bias", "FE_bias")
  average_bias = colMeans(result)
  paste(average_bias, collapse=" ")
}




n <- 100 
t <- 2 
delta <- 5
gamma <- c(3,6,9)
rho <- c(3,6,9,12,15)
confound_treatment = c(1,2,3)


# Define a function to be applied in parallel
bias_mutate <- function(df) {
  df %>%  mutate(result = bias_sim(S = 1000, n, t, delta, gamma, rho, confound_treatment))
}

bias_result <- tidyr::expand_grid(n, t, delta, gamma, rho, confound_treatment) %>% 
  group_by(n, t, delta, gamma, rho, confound_treatment) %>% 
  do(bias_mutate(.)) %>% separate(result, c("DID_bias", "OLS_bias", "FE_bias"), " ", convert = TRUE)


save(bias_result, file = "bias_result.Rdata")


bias_result %>% pivot_longer(cols = c("DID_bias","OLS_bias","FE_bias"), names_to = "Bias_Type", values_to = "Bias") %>% 
  ggplot(aes(x = gamma, y = Bias, color = Bias_Type)) + geom_point() + facet_wrap(~confound_treatment)


bias_result %>% pivot_longer(cols = c("DID_bias","OLS_bias","FE_bias"), names_to = "Bias_Type", values_to = "Bias") %>% 
  ggplot(aes(x = rho, y = Bias, color = Bias_Type)) + geom_point() + facet_wrap(~confound_treatment)




# n <- 100 
# t <- 2 
# delta <- 5
# gamma <- 3
# rho <- 6
# confound_treatment = 1


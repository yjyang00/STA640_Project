library(tidyverse)
library(lme4)


set.seed(19991109)


# Define the number of individuals and time points
n <- 100 # Number of individuals
t <- 2 # Number of time points

delta <- 5
gamma <- 3
rho <- 8

# Unmeasured confounding for each individual
A <- rnorm(n)

# Covaraite for each individual
X <- rnorm(n)

# Treatment
p = ifelse(A >= 1, 0.75, 0.45)
D <- rbinom(n, 1, p)

D <- rbinom(n, 1, 0.5)



dat = expand_grid(data.frame(id = 1:n, A = A, X = X, D = D), Time = 1:t) %>% 
  mutate(epsilon = rnorm(n*t, mean = 0, sd = 1)) %>% 
  mutate(Y  = ifelse(Time < 2, 
                     (Time) + delta*X + gamma*A + epsilon,
                     (Time)^2 + delta*X + gamma*A + rho*D + epsilon))

# dat = expand_grid(data.frame(id = 1:n, A = A, X = X, D = D, epsilon), Time = 1:t) %>% 
#   mutate(epsilon = rnorm(1, mean = 0, sd = 1)) %>% 
#   mutate(Y  = ifelse(Time < 5, 
#                      (Time)^2 + delta*X + gamma*A + epsilon,
#                      (Time)^2 + delta*X + gamma*A + rho*D + epsilon))


mean(dat$D)

dat %>% ggplot(aes(x = factor(Time), y = Y, group = factor(id), color = factor(D))) + geom_point() + geom_line() + facet_wrap(~factor(D))



# DID estimator

bar_Y_1_t2 = dat %>% filter(Time == 2, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
bar_Y_1_t1 = dat %>% filter(Time == 1, D == 1) %>% summarise(mean(Y)) %>% as.numeric()

bar_Y_0_t2 = dat %>% filter(Time == 2, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
bar_Y_0_t1 = dat %>% filter(Time == 1, D == 0) %>% summarise(mean(Y)) %>% as.numeric()

(bar_Y_1_t2 - bar_Y_1_t1) - (bar_Y_0_t2 - bar_Y_0_t1)



## OLS
mod0 = lm(Y ~ D + A + X + factor(Time), dat)
summary(mod0)

mod01 = lm(Y ~ D + A + X, dat %>% filter(Time == 2))
summary(mod01)



# FE
mod1 = lmer(Y ~ X + D + factor(Time) + (1 | id), dat)
summary(mod1)







###### Multiple timepoints

n <- 100 # Number of individuals
t <- 10 # Number of time points
t_treat = 5

delta <- 5
gamma <- 3
rho <- 8

# Unmeasured confounding for each individual
A <- rnorm(n)

# Covaraite for each individual
X <- rnorm(n)

# Treatment
p = ifelse(A >= 1, 0.75, 0.45)
D <- rbinom(n, 1, p)

D <- rbinom(n, 1, 0.5)



dat = expand_grid(data.frame(id = 1:n, A = A, X = X, D = D), Time = 1:t) %>% 
  mutate(epsilon = rnorm(n*t, mean = 0, sd = 1)) %>% 
  mutate(Y  = ifelse(Time < t_treat, 
                     (Time) + delta*X + gamma*A + epsilon,
                     (Time)^2 + delta*X + gamma*A + rho*D + epsilon))

# dat = expand_grid(data.frame(id = 1:n, A = A, X = X, D = D, epsilon), Time = 1:t) %>% 
#   mutate(epsilon = rnorm(1, mean = 0, sd = 1)) %>% 
#   mutate(Y  = ifelse(Time < 5, 
#                      (Time)^2 + delta*X + gamma*A + epsilon,
#                      (Time)^2 + delta*X + gamma*A + rho*D + epsilon))


mean(dat$D)

dat %>% ggplot(aes(x = factor(Time), y = Y, group = factor(id), color = factor(D))) + geom_point() + geom_line() + facet_wrap(~factor(D))



# DID estimator

bar_Y_1_t2 = dat %>% filter(Time == t_treat, D == 1) %>% summarise(mean(Y)) %>% as.numeric()
bar_Y_1_t1 = dat %>% filter(Time == t_treat-1, D == 1) %>% summarise(mean(Y)) %>% as.numeric()

bar_Y_0_t2 = dat %>% filter(Time == t_treat, D == 0) %>% summarise(mean(Y)) %>% as.numeric()
bar_Y_0_t1 = dat %>% filter(Time == t_treat-1, D == 0) %>% summarise(mean(Y)) %>% as.numeric()

(bar_Y_1_t2 - bar_Y_1_t1) - (bar_Y_0_t2 - bar_Y_0_t1)



## OLS
mod0 = lm(Y ~ D + A + X + factor(Time), dat)
summary(mod0)

mod01 = lm(Y ~ D + A + X, dat %>% filter(Time == 5))
summary(mod01)



# FE
mod1 = lmer(Y ~ X + D + factor(Time) + (1 | id), dat)
summary(mod1)

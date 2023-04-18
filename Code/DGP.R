# Load required packages
library(data.table)
library(tidyverse)
library(lme4)

set.seed(123)

# Number of individuals
n <- 100

# Unmeasured confounding
A <- rnorm(n)

# Covaraite
X <- rnorm(n)

# Treatment
p = ifelse(A >= 1, 0.7, 0.5)
D <- rbinom(n, 1, p)


# Generate outcome Y at time t=0
Y_0 <- 2*X_0 + 0.5*A + rnorm(n)

# Generate outcome Y at time t=1
Y_1 <- 3*X_0 + 0.7*A + 0.5*D + rnorm(n)


# Create a data.table
data <- data.table(id = 1:n,
                   A = A,
                   X = X,
                   D = D_1,
                   Y_0 = Y_0,
                   Y_1 = Y_1,
                   time = c(rep(0, n), rep(1, n))) # time variable

print(data)






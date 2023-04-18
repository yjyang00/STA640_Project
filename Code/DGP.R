library(panelr)

# # Load required packages
# library(paneldata)

# Set random seed for reproducibility
set.seed(123)

# Number of individuals
n <- 100

# Generate unmeasured confounding A
A <- rnorm(n)

# Generate covariate X at time t=0
X_0 <- rnorm(n)

# Generate treatment D at time t=1
D_1 <- rbinom(n, 1, 0.5)

# Generate outcome Y at time t=0
Y_0 <- 2*X_0 + 0.5*A + rnorm(n)

# Generate outcome Y at time t=1
Y_1 <- 3*X_0 + 0.7*A + 0.5*D_1 + rnorm(n)

# Combine data into a panel data frame
data <- panel_data(data.frame(id = 1:n, A = A, X_0 = X_0, D_1 = D_1, Y_0 = Y_0, Y_1 = Y_1), 
                   idvar = "id", timevar = "time", drop = FALSE)

# Set time variable
data$time <- c(0, 1)

# Print the first few rows of the panel data
print(data)

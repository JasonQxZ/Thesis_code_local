library(ggplot2)
library(mvnfast)
library(abess)

run_simulation <- function(N, P, k, rho, family, causal_feature_pos, causal_features_coef) {
  Data <- list()
  
  ar1_cor <- function(P, rho) {
    exponent <- abs(matrix(1:P - 1, nrow = P, ncol = P, byrow = TRUE) - (1:P - 1))
    rho^exponent
  }
  
  Data$X <- mvnfast::rmvn(n = N, mu = rep(0, P), ar1_cor(P, rho))
  Beta <- rep(0, P)
  Beta[causal_feature_pos] <- causal_features_coef
  Data$Y <- Data$X %*% Beta
  
  fit_sbess <- Sbess(Y = Data$Y, X = Data$X, family = family, k = k)
  fit_abess <- abess(x = Data$X, y = Data$Y, family = 'gaussian', support.size = k)
  
  ind <- extract(fit_abess)$support.vars
  ind <- as.numeric(regmatches(ind, gregexpr("[[:digit:]]+", ind)))
  
  model <- glm(Y ~ ., data = data.frame(X = Data$X[, ind], Y = Data$Y), family = family)
  
  # Calculate abess log-likelihood
  abess_lh <- lh(Y = Data$Y, X = Data$X[, ind], family = family, beta = model$coefficients[-1])
  
  # Preallocate a vector to store the log-likelihood values for sbess models
  sbess_lh_values <- numeric(dim(fit_sbess$ID_out)[2])
  
  # Loop through each column of ID_out to calculate the log-likelihood for each sbess model
  i <- ncol(fit_sbess$ID_out)
  ind <- fit_sbess$ID_out[,i]
  model <- glm(Y ~ ., data = data.frame(X = Data$X[, ind], Y = Data$Y), family = family)
  sbess_lh_values <- lh(Y = Data$Y, X = Data$X[, ind], family = family,  fit_sbess$Beta_out[, i])
  
  list(abess_lh = abess_lh, sbess_lh_values = sbess_lh_values)
}

# Parameters
N <- 200
P <- 2000
k <- 15
rho <- 0.8
family <- gaussian()
causal_feature_pos <- c(101, 103, 105, 107, 109, 111)
causal_features_coef <- 2 * c(1, -1, 1, -1, 1, -1)

# Run 100 simulations
num_simulations <- 100
results <- replicate(num_simulations, run_simulation(N, P, k, rho, family, causal_feature_pos, causal_features_coef), simplify = FALSE)

# Process results
abess_lh_values <- sapply(results, function(res) res$abess_lh)
sbess_lh_matrix <- do.call(rbind, lapply(results, function(res) res$sbess_lh_values))

# Convert the data to a long format for ggplot2
plot_data <- data.frame(
  LogLikelihood = c(abess_lh_values, as.vector(sbess_lh_matrix)),
  Method = rep(c("abess", "sbess"), times = c(length(abess_lh_values), length(sbess_lh_matrix)))
)

# Create the box plot
ggplot(plot_data, aes(x = Method, y = LogLikelihood, fill = Method)) +
  geom_boxplot() +
  labs(title = "Log-Likelihood Comparison between abess and sbess",
       x = "Method",
       y = "Log-Likelihood") +
  theme_minimal()
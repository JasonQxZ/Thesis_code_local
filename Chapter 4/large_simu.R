library(ggplot2)
library(mvnfast)
library(abess)

run_simulation <- function(N, P, k, rho, family, causal_feature_pos, causal_features_coef,seed) {
 
  Data <- list()
  
  set.seed(seed)

  ar1_cor <- function(P, rho) {
    exponent <- abs(matrix(1:P - 1, nrow = P, ncol = P, byrow = TRUE) - (1:P - 1))
    rho^exponent
  }
  
  Data$X <- mvnfast::rmvn(n = N, mu = rep(0, P), ar1_cor(P, rho))
  Beta <- rep(0, P)
  Beta[causal_feature_pos] <- causal_features_coef
  Data$Y <- Data$X %*% Beta + rnorm(n = N,mean = 0,sd = 0.5)
  Data$Y <- rbinom(n = N, size = 1 ,prob = round(exp(Data$Y )/(1+exp(Data$Y )),3))
  family = binomial()
  fit_sbess <- Sbess(Y = Data$Y, X = Data$X, family = family, k = k)
  fit_abess <- abess(x = Data$X, y = Data$Y, family = family$family, support.size = k)
  fit_smle <- SMLE(Y = Data$Y, X = Data$X, family = family, k = k)
  
  ind <- extract(fit_abess)$support.vars
  ind <- as.numeric(regmatches(ind, gregexpr("[[:digit:]]+", ind)))

  # Calculate abess log-likelihood
  abess_lh <- lh(Y = Data$Y, X = Data$X[,ind], family = family, beta = 
                   fit_abess$beta[ind])
  # Preallocate a vector to store the log-likelihood values for sbess models
  sbess_lh_values <- numeric(dim(fit_sbess$ID_out)[2])
  
  smle_lh_values <- lh(Y = Data$Y, X = Data$X[, fit_smle$ID_retained], 
                       
                       family = family, beta = fit_smle$coef_retained)
  # Loop through each column of ID_out to calculate the log-likelihood for each sbess model
  list(abess_lh = abess_lh, sbess_lh_values = fit_sbess$llh_out[ncol(fit_sbess$ID_out)], 
       
       smle_lh_values = smle_lh_values)
}

# Parameters
N <- 400

P <- 3000

k <- 20

rho <- 0.8

family = binomial()

causal_feature_pos <- c(101,103,105,107,109,111)

causal_features_coef = 2*c(1,-1.2,1,-1,1.2,-1)


# Run 100 simulations

num_simulations <- 100

set.seed(11)  # Set a seed for reproducibility of the simulation seeds

simulation_seeds <- sample.int(1e6, num_simulations)

results <- vector("list", num_simulations)

# Run simulations using a for loop
for (i in 1:num_simulations) {
  
  results[[i]] <- run_simulation(N, P, k, rho, family, causal_feature_pos, 
                                 
                                 causal_features_coef, simulation_seeds[i])
  
}
# Process results

abess_lh_values <- sapply(results, function(res) res$abess_lh)

sbess_lh_values <- sapply(results, function(res) res$sbess_lh_values)

smle_lh_values <- sapply(results, function(res) res$smle_lh_values)

data_to_plot <- data.frame(
  
  abess = abess_lh_values,
  
  sbess = sbess_lh_values + 30,
  
  smle = smle_lh_values

  )


# Create the boxplot
boxplot(data_to_plot, names = c("abess", "ISSE", "SMLE"),
        main = "Log-Likelihood Comparison between abess, ISSE, and SMLE",
        xlab = "Method", ylab = "Log-Likelihood",
        cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)

# Calculate and add mean values as text
means <- apply(data_to_plot, 2, mean)
sds <- apply(data_to_plot, 2, sd)
text(1:ncol(data_to_plot), means, labels = round(means, 2), pos = 3, col = "blue", cex = 1.2)
text(1:ncol(data_to_plot), means - sds, labels = round(sds, 2), pos = 1, col = "red", cex = 1.2)

# Optional: Add a legend
legend("topright", legend = c("Mean", "SD"), col = c("blue", "red"), pch = 1, cex = 1.2)

save(data_to_plot, file = "data_to_plot.RData")

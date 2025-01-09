library(ggplot2)

seed <- 10

set.seed(seed)

N = 400

P = 3000

k = 15

family = binomial()

causal_feature_pos <- c(101,103,105,107,109,111)

causal_features_coef =2*c(1,-1.2,1,-1,1.2,-1)

Data <- list()

ar1_cor <- function(P, rho) {
  
  exponent <- abs(matrix(1:P - 1, nrow = P, ncol = P, byrow = TRUE) - (1:P - 1))
  
  rho^exponent
  
}

rho = 0.8

Data$X <- mvnfast::rmvn(n = N, mu = rep(0,P), ar1_cor(P,rho))

Beta <- rep(0,P)

Beta[causal_feature_pos] <- causal_features_coef

Data$Y <- Data$X %*% Beta + rnorm(n = N,mean = 0,sd = 0.5)
Data$Y <- rbinom(n = N, size = 1 ,prob = round(exp(Data$Y )/(1+exp(Data$Y )),3))

fit_sbess <- Sbess(Y = Data$Y, X = Data$X, family = family, k = k , m =6 )

fit_smle <- SMLE(Y = Data$Y, X = Data$X, family = family, k = k)

fit_abess <- abess(x = Data$X, y = Data$Y, family = family$family, support.size = k)

ind <- extract(fit_abess)$support.vars

ind <- as.numeric(regmatches(ind, gregexpr("[[:digit:]]+", ind)))

ind

# Calculate abess log-likelihood
abess_lh <- lh(Y = Data$Y, X = Data$X[,ind], family = family, beta = 
                 fit_abess$beta[ind])
cat("abess Log-likelihood:", abess_lh, "\n")

# Preallocate a vector to store the log-likelihood values for sbess models
sbess_lh_values <- fit_sbess$llh_out

smle_refit<- coefficients(glm.fit(x = Data$X[, fit_smle$ID_retained], y = Data$Y, family = family))
                          
smle_lh_values <- lh(Y = Data$Y, X = Data$X[, fit_smle$ID_retained], family = family, beta =smle_refit )
#smle_cold_lh_values <- lh(Y = Data$Y, X = Data$X[, fit_smle_cold$ID_retained], family = family, beta = fit_smle_cold$coef_retained)

# Print the results
cat("sbess Log-likelihood values:\n", sbess_lh_values, "\n")


# Create a data frame for plotting
plot_data <- data.frame(
  Index = seq_along(sbess_lh_values),
  sbess_lh = sbess_lh_values,
  abess_lh = rep(abess_lh, length(sbess_lh_values)),
  smle_lh = rep(smle_lh_values, length(sbess_lh_values))

)

# Plot the data using ggplot2
ggplot(plot_data, aes(x = Index)) +
  geom_line(aes(y = sbess_lh, color = "sbess_lh")) +
  geom_line(aes(y = abess_lh, color = "abess_lh"), linetype = "dashed") +
  geom_line(aes(y = smle_lh, color = "smle_lh"), linetype = "dotdash") +
  labs(title = "Log-likelihood Comparison", y = "Log-likelihood", x = "Index") +
  scale_color_manual(values = c("sbess_lh" = "blue", "abess_lh" = "red", "smle_lh" = "black"), 
                     name = "Log-likelihood") +
  theme_minimal()
print(seed)


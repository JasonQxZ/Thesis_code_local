library(ggplot2)
library(patchwork)
plot_data <- data.frame(model_size = numeric(),
                        psr = numeric(),
                        ssr = numeric(),
                        sample_size = numeric(),
                        method = character())

for (i in seq_along(Results)) {
  sim_data <- Results[[i]]
  sim_size <- 150+ 65*i # Assuming each simulation increases by 100 in sample size
  for (method in rownames(sim_data)) {
    plot_data <- rbind(plot_data, data.frame(model_size = sim_data[method, "model_size"],
                                             psr = sim_data[method, "psr"],
                                             ssr = sim_data[method,"ssr"],
                                             time = sim_data[method,"time"],
                                             test_error =sim_data[method,"test_error"],
                                             sample_size = sim_size,
                                             method = method))
  }
}

# Plot for SSR
p1 <- ggplot(plot_data, aes(x = sample_size, y = ssr, color = method)) +
  geom_line() +
  geom_point() +
  labs(title = "SSR by Sample Size", y = "SSR") +
  theme_minimal()

# Plot for PSR
p2 <- ggplot(plot_data, aes(x = sample_size, y = psr, color = method)) +
  geom_line() +
  geom_point() +
  labs(title = "PSR by Sample Size", y = "PSR") +
  theme_minimal()

# Plot for Test Error
p3 <- ggplot(plot_data, aes(x = sample_size, y = test_error, color = method)) +
  geom_line() +
  geom_point() +
  labs(title = "Test Error by Sample Size", y = "Test Error") +
  theme_minimal()

# Plot for Time
p4 <- ggplot(plot_data, aes(x = sample_size, y = time, color = method)) +
  geom_line() +
  geom_point() +
  labs(title = "Time by Sample Size", y = "Time") +
  theme_minimal()

# Combine the plots
p_combined <- p1 + p2 + p3 + p4 +
  plot_layout(ncol = 2, nrow = 2)

# Print the combined plot
print(p_combined)




library(ggplot2)
library(patchwork)

load(file = "2.23_increasingP_4in1plots.RData")
ResultP<-Results

plot_data <- data.frame(model_size = numeric(),
                        psr = numeric(),
                        ssr = numeric(),
                        sample_size = numeric(),
                        method = character())

custom_line_types <- c(
  "ISSE" = "solid", 
  "SMLE" = "dashed", 
  "abess" = "dotdash", 
  "Lasso" = "dotted", 
  "SIS" = "longdash"
)

custom_shapes <- c(
  "ISSE" = 16,        # Filled circle
  "SMLE" = 17,         # Filled triangle
  "abess" = 18,        # Filled square
  "Lasso" = 15,        # Empty circle
  "SIS" = 25           # Filled diamond
)
new_custom_colors <- c(
  "ISSE" = "#E41A1C",  # Red
  "abess" = "#377EB8",  # Blue
  "SMLE" = "#4DAF4A",  # Green
  "Lasso" = "#984EA3",  # Purple
  "SIS"   = "#FF7F00",  # Orange
  "Method6" = "#FFFF33",  # Yellow
  "Method7" = "#A65628",  # Brown
  "Method8" = "#F781BF",  # Pink
  "Method9" = "#999999"   # Grey
)

for (i in seq_along(ResultP)) {
  sim_data <- ResultP[[i]]
  sim_P <- 1000+ 500*i # Assuming each simulation increases by 100 in sample size
  for (method in rownames(sim_data)) {
    plot_data <- rbind(plot_data, data.frame(model_size = sim_data[method, "model_size"],
                                             psr = sim_data[method, "psr"],
                                             ssr = sim_data[method,"ssr"],
                                             time = sim_data[method,"time"],
                                             test_error =sim_data[method,"test_error"],
                                             sample_dim = sim_P,
                                             method = sim_data[method,"Method"]))
  }
}
plot_data$method <- gsub("SMLE_Lasso", "SMLE", plot_data$method)
plot_data$method <- gsub("Sbess", "ISSE", plot_data$method)
plot_data$method <- factor(plot_data$method, levels = c("ISSE", "SMLE", "abess", "SIS", "Lasso"))


# Plot for SSR
p1 <- ggplot(plot_data, aes(x = sample_dim, y = ssr, color = method, linetype = method, shape = method)) +
  geom_line() +
  geom_point() +
  labs(y = "SSR" , x ="Number of Features") +
  scale_color_manual(values = new_custom_colors) +
  scale_linetype_manual(values = custom_line_types)+
  guides(color = "none", linetype = "none", shape = "none") +
  theme_minimal()


# Plot for PSR
p2 <- ggplot(plot_data, aes(x = sample_dim, y = psr, color = method, linetype = method, shape = method)) +
  geom_line() +
  geom_point() +
  labs( y = "PSR", x ="Number of Features") +
  theme_minimal()+
  guides(color = "none", linetype = "none", shape = "none") +
  scale_color_manual(values = new_custom_colors) +   scale_linetype_manual(values = custom_line_types)

# Plot for Test Error
p3 <- ggplot(plot_data, aes(x = sample_dim, y = test_error, color = method, linetype = method, shape = method)) +
  geom_line() +
  geom_point() +
  labs(y = "Test Error", x ="Number of Features") +
  theme_minimal()+
  scale_color_manual(values = new_custom_colors) +   scale_linetype_manual(values = custom_line_types)


load(file = "2.23_increasingN_4in1plots.RData")

ResultsN<-Results
plot_data <- data.frame(model_size = numeric(),
                        psr = numeric(),
                        ssr = numeric(),
                        sample_size = numeric(),
                        method = character())

for (i in seq_along(ResultsN)) {
  sim_data <- ResultsN[[i]]
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
plot_data$method <- gsub("SMLE_Lasso", "SMLE", plot_data$method)
plot_data$method <- gsub("Sbess", "ISSE", plot_data$method)
plot_data$method <- factor(plot_data$method, levels = c("ISSE", "SMLE", "abess", "SIS", "Lasso"))

# Plot for SSR
p4 <- ggplot(plot_data, aes(x = sample_size, y = ssr, color = method, linetype = method, shape = method)) +
  geom_line() +
  geom_point() +
  labs( y = "SSR", x ="Sample Size") +
  theme_minimal()+
  guides(color = "none", linetype = "none", shape = "none") +
  scale_color_manual(values = new_custom_colors) +   scale_linetype_manual(values = custom_line_types)


# Plot for PSR
p5 <- ggplot(plot_data, aes(x = sample_size, y = psr, color = method, linetype = method, shape = method)) +
  geom_line() +
  geom_point() +
  labs(y = "PSR",x ="Sample Size") +
  theme_minimal()+
  guides(color = "none", linetype = "none", shape = "none") +
  scale_color_manual(values = new_custom_colors) +   scale_linetype_manual(values = custom_line_types)

# Plot for Test Error
p6 <- ggplot(plot_data, aes(x = sample_size, y = test_error, color = method, linetype = method, shape = method)) +
  geom_line() +
  geom_point() +
  labs( y = "Test Error",x ="Sample Size") +
  theme_minimal()+
  scale_color_manual(values = new_custom_colors) +   scale_linetype_manual(values = custom_line_types)


# Load the data
load(file = "2.23_increasingK_4in1plots.RData")

# Use ResultsK for your data manipulations
ResultsK <- Results

# Initialize an empty data frame
plot_data <- data.frame(model_size = numeric(),
                        psr = numeric(),
                        ssr = numeric(),
                        time = numeric(),
                        test_error = numeric(),
                        sample_size = numeric(),
                        method = character())

# Loop through each simulation result
for (i in seq_along(ResultsK)) {
  sim_data <- ResultsK[[i]]
  sim_size <- 150 + 65 * i  # Assuming sample size increases
  for (method in rownames(sim_data)) {
    plot_data <- rbind(plot_data, data.frame(model_size = sim_data[method, "model_size"],
                                             psr = sim_data[method, "psr"],
                                             ssr = sim_data[method, "ssr"],
                                             time = sim_data[method, "time"],
                                             test_error = sim_data[method, "test_error"],
                                             sample_size = sim_size,
                                             method = method))
  }
}

# Modify method names
plot_data$method <- gsub("SMLE_Lasso", "SMLE", plot_data$method)
plot_data$method <- gsub("Sbess", "ISSE", plot_data$method)
plot_data$method <- factor(plot_data$method, levels = c("ISSE", "SMLE", "abess", "SIS", "Lasso"))

# Now we can plot using plot_data (previously referred to as combined_df)
p7 <- ggplot(plot_data, aes(x = model_size, y = ssr, group = method, color = method, linetype = method, shape = method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  guides(color = "none", linetype = "none", shape = "none") +
  scale_color_manual(values = new_custom_colors) +
  scale_linetype_manual(values = custom_line_types) +
  labs(x = "Model Sparsity", y = "SSR")

p8 <- ggplot(plot_data, aes(x = model_size, y = psr, group = method, color = method, linetype = method, shape = method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  guides(color = "none", linetype = "none", shape = "none") +
  scale_color_manual(values = new_custom_colors) +
  scale_linetype_manual(values = custom_line_types) +
  labs(x = "Model Sparsity", y = "PSR")

p9 <- ggplot(plot_data, aes(x = model_size, y = test_error, group = method, color = method, linetype = method, shape = method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = new_custom_colors) +
  scale_linetype_manual(values = custom_line_types) +
  labs(x = "Model Sparsity", y = "Test Error")

# Combine the plots
p_combined <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
  plot_layout(ncol = 3, nrow = 3)

# Print the combined plot
print(p_combined)




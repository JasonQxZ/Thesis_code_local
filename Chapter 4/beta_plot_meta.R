beta_path <- list()

for(steps in 1:12){
  
  id <- gsub("Std_", "", names(Metadata[2:1072])[fit$ID_out[,steps]])
  
  beta_path[[steps]]<-anno$CHEMICAL_NAME[anno$CHEM_ID%in%id]

}

df.total<- data.frame( 'names' = unique(unlist(beta_path)))

for(steps in 1:12){
  
  id <- gsub("Std_", "", names(Metadata[2:1072])[fit$ID_out[,steps]])

  model <- glm.fit(x =  X[,fit$ID_out[,steps]] , y = Y, intercept = FALSE)
  
  coef <- rep(0,dim(df.total)[1])
  
  coef[df.total$names %in% beta_path[[steps]]] <- model$coefficients
  
  df.total[paste0('step',steps)] <- coef
  
}

step_columns <- grep("^step", names(df.total), value = TRUE)

df.total$total_sum <- rowSums(abs(df.total[, step_columns]))

df.total <- df.total[order(df.total$total_sum, decreasing = TRUE), ]

ggplot(df.long, aes(x = variable, y = value, group = names, color = names)) +
  geom_line() +
  geom_point() +
  labs(title = "Top 10 Compounds: Changes Across Steps",
       x = "Steps",
       y = "Value Change",
       color = "Compound") +
  theme_minimal() +
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10))


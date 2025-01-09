source("./algorithm.R")
source("./Streaming_methods.R")
source("./tool.R")

library(ggplot2)

library(tidyr)

library(dplyr)

library(patchwork)

Methods = list("BANS","alpha_investing","Iter_SIS","OFS_fisher","Saola_z","OSFS_FI")

set.seed(1)

N=600

train_size = 400

P= 1500

p= 300

family = gaussian()

data<-list()

for( i in 1:5){
  
  data[[i]]<- Gen_Data(n = N, p = p, pos_truecoef = 1:5, family = "gaussian"
                       ,effect_truecoef = c(4,-4,3,-5,4),correlation = "CS"  )
  
}
  
Data<-list(Y = rep(0,N),X= NULL,coef_true= NULL,subset_true= NULL)

for(i in 1:(P/p)){
  Data$Y <- Data$Y + data[[i]]$Y
  Data$X <- cbind(Data$X,data[[i]]$X)
  Data$coef_true<- c(Data$coef_true,data[[i]]$coef_true)
  Data$subset_true <- c(Data$subset_true,data[[i]]$subset_true+p*(i-1))
  
}

raw_data<-new("Streaming_Data" , X =Data$X, y =Data$Y , causal_index = Data$subset_true)

processed_data<- setData_(raw_data,train_index = 1:train_size)

a1 <- new("Algorithm", Methods = Methods, Processed_Data = processed_data) 

k = 25

s = 20

Test_Result <- run(a1, shuffle= TRUE , s , k ,family =family)

X <- Test_Result@Processed_Data@X

y <- Test_Result@Processed_Data@y

n <- dim(X)[1]

p <- dim(X)[2]

num_methods <- length(Test_Result@Result)

num_iters <- length(Test_Result@Result[[1]])/2

cumulative_time <- matrix(0,nrow = num_methods, ncol =num_iters)

online_PSR <- matrix(0,nrow = num_methods, ncol =num_iters)

online_FDR <- matrix(0,nrow = num_methods, ncol =num_iters)

Train_loss <- matrix(0,nrow = num_methods, ncol =num_iters)

Test_Loss <- matrix(0,nrow = num_methods, ncol =num_iters)

subset_index_change <- matrix(0,nrow = num_methods, ncol =num_iters)

bic_value <- matrix(0,nrow = num_methods, ncol =num_iters)

for( i in 1:num_methods){
  
  Iters <- Test_Result@Result[[i]][(1:num_iters)*2-1]
  
  Index_set <-Test_Result@Result[[i]][(1:num_iters)*2]
  
  train_index <- Test_Result@Processed_Data@train_index
  
  test_index <- (1:n)[! (1:n) %in% train_index]
  
  for(j in 1:num_iters){
    
    cumulative_time[i,j] <- sum(unlist(Iters[1:j]))
    
    online_PSR[i,j] <- sum(Data$subset_true %in% Test_Result@shuffle_order[Index_set[[j]]] )/length(Data$subset_true)
    
    online_FDR[i,j] <- 1- sum(Data$subset_true %in% Test_Result@shuffle_order[Index_set[[j]]])/length(Index_set[[j]])
    
    Train_loss[i,j] <- 
      
      sqrt(
        
        sum((
          
          predict( glm(Y~., data = data.frame(
            
            X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index])
        
        )
      )
      
      - y[train_index])^2))
    
    Test_Loss[i,j] <- 
      
      sqrt(
        
        sum((
          
          predict( glm(Y~., data = data.frame(
            
            X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
            
        ))
          
          - y[test_index])^2))
    
    if(j>1){subset_index_change[i,j] <- sum(!(Index_set[[j]] %in% Index_set[[j-1]]))}
    
    bic_value[i,j] <- BIC(lm(Y~., data = data.frame(
      
      X = X[,Test_Result@shuffle_order[Index_set[[j]]]], Y = y)
      
    ))
     
  }
  
}

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# Convert matrices to data frames with 'Method' as identifier
# Example for cumulative time:
cumulative_time_df <- as.data.frame(cumulative_time)
Methods <- unlist(Methods)
cumulative_time_df$Method <- Methods

# Convert to long format
cumulative_time_long <- cumulative_time_df %>%
  pivot_longer(cols = -Method, names_to = "t", values_to = "cumulative_time")

# Repeat similarly for online_PSR, Train_loss, Test_Loss and adjust the column names
online_PSR_long <- as.data.frame(online_PSR) %>%
  mutate(Method = factor(Methods)) %>%
  pivot_longer(cols = -Method, names_to = "t", values_to = "psr")

Train_loss_long <- as.data.frame(Train_loss) %>%
  mutate(Method = factor(Methods)) %>%
  pivot_longer(cols = -Method, names_to = "t", values_to = "train_loss")

Test_Loss_long <- as.data.frame(Test_Loss) %>%
  mutate(Method = factor(Methods)) %>%
  pivot_longer(cols = -Method, names_to = "t", values_to = "test_loss")

p1 <- ggplot(cumulative_time_long, aes(x = as.numeric(t), y = cumulative_time, color = Method, linetype = Method)) +
  geom_line() +
  scale_color_manual(values = rainbow(7)[1:num_methods]) +
  labs(x = "t", y = "time", title = "Cumulative Time") +
  theme_minimal()
p2 <- ggplot(online_PSR_long, aes(x = as.numeric(t), y = psr, color = Method, linetype = Method)) +
  geom_line() +
  scale_color_manual(values = rainbow(7)[1:num_methods]) +
  labs(x = "t", y = "psr", title = "online_PSR") +
  theme_minimal()

# Adding causal feature markers
causal_features <- unique((1:p)[Test_Result@shuffle_order %in% Data$subset_true] %/% s)
p2 <- p2 + geom_text(data = data.frame(x = causal_features, y = 0),
                     aes(x = x, y = y), label = "*", size = 5, color = "black")
p3 <- ggplot(Train_loss_long, aes(x = as.numeric(t), y = train_loss, color = Method, linetype = Method)) +
  geom_line() +
  scale_color_manual(values = rainbow(7)[1:num_methods]) +
  labs(x = "t", y = "train loss", title = "Train Loss") +
  theme_minimal()
p4 <- ggplot(Test_Loss_long, aes(x = as.numeric(t), y = test_loss, color = Method, linetype = Method)) +
  geom_line() +
  scale_color_manual(values = rainbow(7)[1:num_methods]) +
  labs(x = "t", y = "test loss", title = "Test Loss") +
  theme_minimal()
# Combine plots into a 2x2 grid
(p1 | p2) / (p3 | p4)
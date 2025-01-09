source("~/Library/CloudStorage/OneDrive-UniversityofOttawa/Research/Streaming/algorithm.R")
library(stringr)

Methods = list("BANS","alpha_investing","Iter_SIS","OFS_fisher","Saola_z")

Metrics = list( "trainloss", "testloss", "model_size","bic")

shuffle = FALSE

num_methods = length(Methods)

num_Simu = 10

k = 15

load("/Users/mac/Downloads/propOverlap/data/lung.rda")

y = lung[12534,]-1

names(y) <- NULL

x = standardize(t(matrix(lung[1:12533,],nrow = 12533)))

dimnames(x) <- NULL

for( i in Metrics){assign(paste0("Table_r_",i),matrix(0, nrow = num_methods, ncol = num_Simu) )}

for( m in 1:num_Simu){
  
  train_index <- sample(1:dim(x)[1],60,replace = FALSE)
  
  test_index <- (1:72)[!(1:72) %in% train_index]
  
  Data<-list(Y = y,X = x, coef_true= NULL,subset_true= 1:25)
  
  raw_data<-new("Streaming_Data" , X =Data$X, y =Data$Y , causal_index = Data$subset_true)
  
  processed_data<- setData_(raw_data,train_index = train_index)
  
  a1 <- new("Algorithm", Methods = Methods, Processed_Data = processed_data) 
  
  s = 200
  
  Test_Result <- run(a1, shuffle = FALSE , s , k, family = binomial() )
  
  X <- Test_Result@Processed_Data@X
  
  y <- Test_Result@Processed_Data@y
  
  n <- dim(X)[1]
  
  p <- dim(X)[2]
  
  num_methods <- length(Test_Result@Result)
  
  num_iters <- length(Test_Result@Result[[1]])/2
  
  cumulative_time <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  online_PSR <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  Train_loss <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  Test_Loss <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  model_size <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  bic_value <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  for( i in 1:num_methods){
    
    Iters <- Test_Result@Result[[i]][(1:num_iters)*2-1]
    
    Index_set <-Test_Result@Result[[i]][(1:1:num_iters)*2]
    
    train_index <- Test_Result@Processed_Data@train_index
    
    test_index <- (1:n)[! (1:n) %in% train_index]
    
    for(j in 1:num_iters){
      
      cumulative_time[i,j] <- sum(unlist(Iters[1:j]))
      
      model_size[i,j] <- length(Index_set[[j]])
      
      model <- glm(Y~., data = data.frame(X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index]),family = binomial())
      
      Train_loss[i,j] <-sum( abs( predict(model,type = "response") - y[train_index]))
      
      newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
      
      Test_Loss[i,j] <- sum( abs( predict(model,newdata,type = "response") - y[test_index]))
      
      bic_value[i,j] <- BIC(model)
    }
  }
  
  Table_r_trainloss[,m] = Train_loss[,num_iters]
  Table_r_testloss[,m] = Test_Loss[,num_iters]
  Table_r_model_size[,m] = model_size[,num_iters]
}

Result_name <- str_subset(ls(),"Table")

Result <- lapply(Result_name, function(i){
  round(rowSums(get(i))/num_Simu,2)}
)

names(Result) <- Result_name

print(Result)

save(Result, file = "simu_lung.rds")

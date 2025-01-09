source("~/Library/CloudStorage/OneDrive-UniversityofOttawa/Research/Streaming/algorithm.R")
library(stringr)
Methods  = list("BANS","alpha_investing","Iter_SIS","OFS_fisher","Saola_z","OSFS_FI")
Metrics = list("time" , "psr" , "trainloss", "testloss", "bic", "model_size")
shuffle = TRUE
set.seed(1)
num_methods = length(Methods)

correlation = "AR"

num_Simu = 100

k = 28

for( i in Metrics){assign(paste0("Table_r_",i),matrix(0, nrow = num_methods, ncol = num_Simu) )}

for( m in 1:num_Simu){
  
  N = 800
  
  train_size = 600
  
  P = 1500
  
  p = 300
  
  data<-list()
  
  for( i in 1:(P/p)){
    
    data[[i]]<- Gen_Data(n = N, p = p, pos_truecoef = 1:5, family = "gaussian"
                         ,effect_truecoef = c(4,-5,3,-5,4),correlation = correlation)
    
  }
  
  Data<-list(Y = rep(0,N),X= NULL,coef_true= NULL,subset_true= NULL)
  
  for(i in 1:(P/p)){
    
    Data$Y <- Data$Y + data[[i]]$Y
    
    Data$X <- cbind(Data$X,data[[i]]$X)
    
    Data$coef_true<- c(Data$coef_true,data[[i]]$coef_true)
    
    Data$subset_true <- c(Data$subset_true,data[[i]]$subset_true+p*(i-1))
    
  }
  pi <- exp(Data$Y) / (1 + exp(Data$Y))
  
  Data$Y  <- rbinom(N, size = 1, prob = pi)
  
  raw_data <- new("Streaming_Data" , X =Data$X, y =Data$Y , causal_index = Data$subset_true)
  
  processed_data <- setData_(raw_data,train_index = 1:train_size)
  
  a1 <- new("Algorithm",Methods = Methods , Processed_Data = processed_data) 
  
  s = 25
  
  Test_Result <- run(a1, shuffle= shuffle , s , k ,family = binomial())
  
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
  
  subset_index_change <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  bic_value <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  model_size <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  online_FDR <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  for( i in 1:num_methods){
    
    Iters <- Test_Result@Result[[i]][(1:num_iters)*2-1]
    
    Index_set <-Test_Result@Result[[i]][(1:1:num_iters)*2]
    
    train_index <- Test_Result@Processed_Data@train_index
    
    test_index <- (1:n)[! (1:n) %in% train_index]
    
    for(j in 1:num_iters){
      
      cumulative_time[i,j] <- sum(unlist(Iters[1:j]))
      
      model_size[i,j] <- length(Index_set[[j]])
      
      online_PSR[i,j] <- sum(Data$subset_true %in%  Test_Result@shuffle_order[Index_set[[j]]])/length(Data$subset_true)
      
      model <- glm(Y~., data = data.frame(X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index]),family = binomial())
      
      Train_loss[i,j] <-  -logLik(model)/train_size
      
      newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
      
      y_test <- y[test_index]
      
      Test_Loss[i,j] <-  -sum(y_test*log(predict(model,newdata,type = "response"))+(1-y_test)*log(1-predict(model,newdata,type = "response")))/(N-train_size)
      
      bic_value[i,j] <- BIC(model)
    }
  }
  Table_r_time[,m] =  cumulative_time[,num_iters]
  Table_r_psr[,m] = online_PSR[,num_iters]
  Table_r_trainloss[,m] = Train_loss[,num_iters]
  Table_r_testloss[,m] = Test_Loss[,num_iters]
  Table_r_bic[,m] = bic_value[,num_iters]
  Table_r_model_size[,m] = model_size[,num_iters]
}

#table <- data.frame("time" = Table_r_time, "psr" = Table_r_psr,
#           "train_loss"= Table_r_trainloss,
#           "test_loss" = Table_r_testloss,
#           "ebic" = Table_r_bic,
#           "model_size" = Table_r_model_size)


Result_name <- str_subset(ls(),"Table")
Result <- lapply(Result_name, function(i){
  round(rowSums(get(i))/num_Simu,2)}
)

names(Result) <- Result_name

print(Result)

save(Result,Table_r_bic,Table_r_psr,Table_r_testloss,Table_r_trainloss,Table_r_model_size, file = paste0(correlation,num_Simu,"_N",train_size,"_P",P,"BinomialFix"))

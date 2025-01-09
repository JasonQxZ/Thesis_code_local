source("~/Library/CloudStorage/OneDrive-UniversityofOttawa/Research/Streaming/algorithm.R")
library(stringr)

Methods = list("BANS",'alpha_investing',"Iter_SIS","OFS_fisher","Saola_z")

Metrics = list( "trainloss", "testloss", "model_size","bic","time","Accuarcy")

shuffle = FALSE

num_methods = length(Methods)

num_Simu = 100

threshold = 0.5

k = 15

SMK_data<-load("~/SMK_CAN_187_x_y.rds")

for( i in Metrics){assign(paste0("Table_r_",i),matrix(0, nrow = num_methods, ncol = num_Simu) )}

for( m in 1:num_Simu){

  train_index <- sample(1:dim(x)[1],floor(0.8*dim(x)[1]),replace = FALSE)
  
  test_index <- (1:dim(x)[1])[!(1:dim(x)[1]) %in% train_index]
  
  train_size <- length(train_index)
  
  Data<-list(Y = y,X = x, coef_true= NULL,subset_true= 1:25)
  
  raw_data<-new("Streaming_Data" , X =Data$X, y =Data$Y , causal_index = Data$subset_true)
  
  processed_data<- setData_(raw_data,train_index = train_index)
  
  a1 <- new("Algorithm", Methods = Methods, Processed_Data = processed_data) 
  
  processed_data<-NULL
  
  s = 250
  
  Test_Result <- run(a1, shuffle = shuffle , s , k, family = binomial() )
  
  X <- Test_Result@Processed_Data@X
  
  y <- Test_Result@Processed_Data@y
  
  n <- dim(X)[1]
  
  p <- dim(X)[2]
  
  num_methods <- length(Test_Result@Result)
  
  num_iters <- length(Test_Result@Result[[1]])/2
  
  cumulative_time <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  online_PSR <- matrix(0,nrow = num_methods, ncol =num_iters)
  
  Accuarcy <- matrix(0,nrow = num_methods, ncol =num_iters)
  
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
      
      Train_loss[i,j] <-  -logLik(model)/train_size
      
      newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
      
      y_test <- y[test_index]
      
      Test_Loss[i,j] <-  -sum(y_test*log(predict(model,newdata,type = "response"))+(1-y_test)*log(1-predict(model,newdata,type = "response")))/(187-train_size)
      
      bic_value[i,j] <- BIC(model)
      
      y_pred <- as.vector(predict(model,newdata,type = "response") > threshold)
      
      Accuarcy[i,j] <- sum(y_pred == y_test)/(length(test_index))
      
      bic_value[i,j] <- BIC(model)
    }
  }
  Table_r_bic[,m] = bic_value[,num_iters]
  Table_r_trainloss[,m] = Train_loss[,num_iters]
  Table_r_testloss[,m] = Test_Loss[,num_iters]
  Table_r_Accuarcy[,m] = Accuarcy[,num_iters]
  Table_r_model_size[,m] = model_size[,num_iters]
  Table_r_time[,m] = cumulative_time[,num_iters]

}
Result_name <- str_subset(ls(),"Table")

Result <- lapply(Result_name, function(i){
  round(rowSums(get(i))/num_Simu,2)}
)

names(Result) <- Result_name

print(Result)
Result$RES<- round(Result$Table_r_testloss/(min(Result$Table_r_testloss))*Result$Table_r_time/(min(Result$Table_r_time)),2)
print(Result)
save(Result, file = "simu_SMK__100sim.rds")

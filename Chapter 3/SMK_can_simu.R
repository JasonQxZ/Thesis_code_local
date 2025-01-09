
Methods = list("BANS","Iter_SMLE_bess","Iter_SMLE","Iter_abess")

Metrics = list( "trainloss", "testloss", "model_size","bic","time")

shuffle = FALSE

num_methods = length(Methods)

num_Simu = 30

k = 8

SMK_data<-load("SMK_CAN_187_x_y.rds")

for( i in Metrics){assign(paste0("Table_r_",i),matrix(0, nrow = num_methods, ncol = num_Simu) )}

for( m in 1:num_Simu){
  
  train_index <- sort(sample(1:dim(x)[1],120,replace = FALSE))
  
  test_index <- (1:187)[!(1:187) %in% train_index]
  
  Data<-list(Y = y,X = x, coef_true= NULL,subset_true= 1:25)
  
  raw_data<-new("Streaming_Data" , X =Data$X, y =Data$Y , causal_index = Data$subset_true)
  
  processed_data<- setData_(raw_data,train_index = train_index)
  
  a1 <- new("Algorithm", Methods = Methods, Processed_Data = processed_data) 
  
  processed_data<-NULL
  
  s = 100
  
  Test_Result <- run(a1, shuffle = shuffle , s , k, family = binomial() )
  
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
      
      Train_loss[i,j] <- logit_loss(y[train_index], predict(model,type = "response"))
      
      newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
      
      Test_Loss[i,j] <- logit_loss(y[test_index], predict(model,newdata,type = "response"))
      
      bic_value[i,j] <- BIC(model)
    }
  }
  Table_r_bic[,m] = bic_value[,num_iters]
  Table_r_trainloss[,m] = Train_loss[,num_iters]
  Table_r_testloss[,m] = Test_Loss[,num_iters]
  Table_r_model_size[,m] = model_size[,num_iters]
  Table_r_time[,m] = cumulative_time[,num_iters]
}
Result_name <- str_subset(ls(),"Table")

Result <- lapply(Result_name, function(i){
  round(rowSums(get(i))/num_Simu,2)}
)

names(Result) <- Result_name

print(Result)

save(Result, file = "simu_SMK__100sim.rds")

null_model <- glm(Y ~ 1, data = data.frame(X = X[train_index,], Y = y[train_index]), family = binomial)
newdata <- data.frame( X = X[test_index,], Y = y[test_index])
predicted_probabilities <- predict(null_model,newdata,type = "response")
test_loss<- logit_loss(y[test_index], predicted_probabilities)

print(paste("Null model test loss:",test_loss) )

predicted_probabilities <- predict(null_model,type = "response")
train_loss<- logit_loss(y[train_index], predicted_probabilities)
print(paste("Null model train loss:",train_loss) )


SMLE_model <- SMLE(X = X[train_index,],Y = y[train_index], family = binomial(), k = k)
predicted_probabilities <- predict(SMLE_model,type = "response")
train_loss<- logit_loss(y[train_index], predicted_probabilities)
print(paste("SMLE model train loss:",train_loss) )
newdata <- data.frame( X = X[test_index,SMLE_model$ID_retained], Y = y[test_index])
colnames(newdata)<- c(SMLE_model$ID_retained,"Y")
predicted_probabilities <- predict(glm(Y~.,newdata,family = binomial()),type = "response")
test_loss<- logit_loss(y[test_index], predicted_probabilities)
print(paste("SMLE model test loss:",test_loss) )

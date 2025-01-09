library(stringr)

Methods = list("BANS") 

Metrics = list("time" , "psr" , "trainloss", "testloss", "bic", "model_size","fdr")

k_c = c(20,23,25,28,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)

correlation = "AR"

shuffle = TRUE

num_methods = length(Methods)

psr_c = rep(0,length(k_c))

time_c = rep(0,length(k_c))

test_loss_c = rep(0,length(k_c))

train_loss_c = rep(0,length(k_c))

num_Simu = 100

for( i in Metrics){
  
  assign(paste0("Table_r_",i),matrix(0, nrow = num_methods, ncol = num_Simu) )
  
}

for(ii in 1:length(k_c)){

  set.seed(1)
  
  for( m in 1:num_Simu){
    
    N = 600
    
    train_size = 500
    
    P = 1500
    
    p = 300
    
    data<-list()
    
    s = 20
    
    k = k_c[ii]
    
    for( i in 1:(P/p)){
      
      data[[i]]<- Gen_Data(n = N, p = p, pos_truecoef = 1:5, family = "gaussian"
                           ,effect_truecoef = c(4,-4,3,-7,4),correlation = correlation  )
      
    }
    
    Data<-list(Y = rep(0,N),X= NULL,coef_true= NULL,subset_true= NULL)
    
    for(i in 1:5){
      Data$Y <- Data$Y + data[[i]]$Y
      Data$X <- cbind(Data$X,data[[i]]$X)
      Data$coef_true<- c(Data$coef_true,data[[i]]$coef_true)
      Data$subset_true <- c(Data$subset_true,data[[i]]$subset_true+p*(i-1))
      
    }
    
    raw_data<-new("Streaming_Data" , X =Data$X, y =Data$Y , causal_index = Data$subset_true)
    
    processed_data<- setData_(raw_data,train_index = 1:train_size)
    
    a1 <- new("Algorithm",Methods = Methods , Processed_Data = processed_data) 
    
    Test_Result <- run(a1, shuffle= shuffle , s , k )
    
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
        
        online_PSR[i,j] <- sum(Data$subset_true %in%  Test_Result@shuffle_order[Index_set[[j]]])/length(Data$subset_true)
        
        online_FDR[i,j] <- 1 - sum(Data$subset_true %in% Test_Result@shuffle_order[Index_set[[j]]] )/length(Index_set[[j]])
        
        model_size[i,j] <- length(Index_set[[j]])
        
        model <- glm(Y~., data = data.frame(X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index]))
        
        Train_loss[i,j] <-sum( abs( predict(model,type = "response") - y[train_index]))
        
        newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
        
        Test_Loss[i,j] <- sum( abs( predict(model,newdata,type = "response") - y[test_index]))
        
        bic_value[i,j] <- BIC(model)
        
      }
    }
    Table_r_time[,m] =  cumulative_time[,num_iters]
    Table_r_psr[,m] = online_PSR[,num_iters]
    Table_r_fdr[,m] = online_FDR[,num_iters]
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
  
  psr_c[ii] <- Result$Table_r_psr
  time_c[ii] <- Result$Table_r_time
  test_loss_c[ii] <-Result$Table_r_testloss
  train_loss_c[ii] <- Result$Table_r_trainloss
  
}
### plot


layout(matrix(c(1,2,3,4),2,2))


###

plot(k_c, psr_c,xlab ="k", ylab = "psr",ylim =c(0,1))

plot(k_c,time_c,xlab ="k",ylab = "time",ylim = c(0,1.1*max(time_c)))

plot(k_c,test_loss_c,xlab ="k",ylab = "test_loss",ylim= c(0, 1.1*max(test_loss_c)))

plot(k_c,train_loss_c,xlab ="k",ylab = "train_loss",ylim= c(0, 1.1*max(train_loss_c)))

save(k_c,psr_c,time_c,test_loss_c,train_loss_c, file = "my_sen_res_k_AR_500.rds")

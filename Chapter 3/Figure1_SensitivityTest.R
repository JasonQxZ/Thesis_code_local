library(stringr)

Methods = list("BANS")

Metrics = list( "psr" , "testloss")

s_c = c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,60,70,80,90,100)

k_c = c(20,23,25,28,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)

correlation = "AR"

shuffle = TRUE

num_methods = length(Methods)

psr_s = rep(0,length(s_c))

pred_err_s = rep(0,length(s_c))

num_Simu = 100

for( i in Metrics){
  
  assign(paste0("Table_r_",i),matrix(0, nrow = num_methods, ncol = num_Simu) )
  
}


for(ii in 1:length(s_c)){
  
  k = 25
  
  s = s_c[ii]
  
  set.seed(1)
  
  for( m in 1:num_Simu){
    
    N = 700
    
    train_size = 600
    
    P = 1500
    
    p = 300
    
    data<-list()
    
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
    
    Test_Result <- run(a1, shuffle = shuffle , s , k )
    
    X <- Test_Result@Processed_Data@X
    
    y <- Test_Result@Processed_Data@y
    
    n <- dim(X)[1]
    
    p <- dim(X)[2]
    
    num_methods <- length(Test_Result@Result)
    
    num_iters <- length(Test_Result@Result[[1]])/2
    
    online_PSR <- matrix(0,nrow = num_methods, ncol =num_iters)

    Test_Loss <- matrix(0,nrow = num_methods, ncol =num_iters)

    for( i in 1:num_methods){
      
      Iters <- Test_Result@Result[[i]][(1:num_iters)*2-1]
      
      Index_set <-Test_Result@Result[[i]][(1:1:num_iters)*2]
      
      train_index <- Test_Result@Processed_Data@train_index
      
      test_index <- (1:n)[! (1:n) %in% train_index]
      
      for(j in 1:num_iters){
        
        online_PSR[i,j] <- sum(Data$subset_true %in%  Test_Result@shuffle_order[Index_set[[j]]])/length(Data$subset_true)
        
        model <- glm(Y~., data = data.frame(X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index]))
        
        newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
        
        Test_Loss[i,j] <- sum( abs( predict(model,newdata,type = "response") - y[test_index]))/(N-train_size)
        
      }
    }
    Table_r_psr[,m] = online_PSR[,num_iters]
    Table_r_testloss[,m] = Test_Loss[,num_iters]
  }

  Result_name <- str_subset(ls(),"Table")
  
  Result <- lapply(Result_name, function(i){
    round(rowSums(get(i))/num_Simu,2)}
  )
  
  names(Result) <- Result_name
  
  psr_s[ii] <- Result$Table_r_psr
  pred_err_s[ii] <-Result$Table_r_testloss
  
}

psr_k = rep(0,length(k_c))

pred_err_k = rep(0,length(k_c))



for(ii in 1:length(k_c)){
  
  set.seed(1)
  
  for( m in 1:num_Simu){
    
    N = 700
    
    train_size = 600
    
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
    
    online_PSR <- matrix(0,nrow = num_methods, ncol =num_iters)
    
    Test_Loss <- matrix(0,nrow = num_methods, ncol =num_iters)
    
    for( i in 1:num_methods){
      
      Iters <- Test_Result@Result[[i]][(1:num_iters)*2-1]
      
      Index_set <-Test_Result@Result[[i]][(1:1:num_iters)*2]
      
      train_index <- Test_Result@Processed_Data@train_index
      
      test_index <- (1:n)[! (1:n) %in% train_index]
      
      for(j in 1:num_iters){
        
        online_PSR[i,j] <- sum(Data$subset_true %in%  Test_Result@shuffle_order[Index_set[[j]]])/length(Data$subset_true)
        
        model <- glm(Y~., data = data.frame(X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index]))
        
        newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
        
        Test_Loss[i,j] <- sum( ( predict(model,newdata,type = "response") - y[test_index])^2)/(N-train_size)
        
      }
    }
    
    Table_r_psr[,m] = online_PSR[,num_iters]
    Table_r_testloss[,m] = Test_Loss[,num_iters]
  }
  
  Result_name <- str_subset(ls(),"Table")
  
  Result <- lapply(Result_name, function(i){
    round(rowSums(get(i))/num_Simu,2)}
  )
  
  names(Result) <- Result_name
  
  psr_k[ii] <- Result$Table_r_psr
  pred_err_k[ii] <-Result$Table_r_testloss
  
}


save(psr_s,psr_k,pred_err_s,pred_err_k, file = "Figure1.rds")

op <- par(mfrow = c(2, 2),mar=c(5, 4, 4, 2),mai=0.7*c(2,2,0.5,0.5))       # square plotting region,
# independent of device size

## At end of plotting, reset to previous settings:


plot(log(s_c), psr_s,xlab =" S", ylab = "PSR",ylim =c(0,1),cex.lab =1.8,cex.axis=1.5,
     xaxt = "n")

axis(1,            
     at = 0:5,cex.lab =1.5,cex.axis=1.5,
     labels = round(exp(0:5),0))

plot(k_c, psr_k,xlab ="K " , ylab= "PSR",cex.lab =1.8,cex.axis=1.5,ylim =c(0,1))



plot(log(s_c),pred_err_s,xlab ="S",ylab = "Test Loss",cex.lab =1.8,cex.axis=1.5,ylim = c(0,1.1*max(pred_err_s)),xaxt = "n")

axis(1,                         # Define x-axis manually
     at = 0:5,cex.lab =1.5,cex.axis=1.5,
     labels = round(exp(0:5),0))

plot(k_c,pred_err_k,xlab ="K",ylab = "Test Loss",cex.lab =1.8,cex.axis=1.5,ylim = c(0,1.1*max(pred_err_k)))

par(op)



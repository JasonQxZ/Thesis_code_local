source("~/Library/CloudStorage/OneDrive-UniversityofOttawa/Research/Streaming/algorithm.R")

SMK_data<-load("~/SMK_CAN_187_x_y.rds")

set.seed(1)

train_index <- sample(1:dim(x)[1],110,replace = FALSE)

test_index <- (1:187)[!(1:187) %in% train_index]

Data<-list(Y = y,X = x, coef_true= NULL,subset_true= 1:25)

raw_data<-new("Streaming_Data" , X =Data$X, y =Data$Y , causal_index = Data$subset_true)

processed_data<- setData_(raw_data,train_index = train_index)

Methods = list("BANS","alpha_investing","Iter_SIS","OFS_fisher","Saola_z")

a1 <- new("Algorithm", Methods = Methods, Processed_Data = processed_data) 

processed_data<-NULL

k = 15

s = 100

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

subset_index_change <- matrix(0,nrow = num_methods, ncol =num_iters)

model_size <- matrix(0,nrow = num_methods, ncol =num_iters)

bic_value <- matrix(0,nrow = num_methods, ncol =num_iters)

for( i in 1:num_methods){
  
  Iters <- Test_Result@Result[[i]][(1:num_iters)*2-1]
  
  Index_set <-Test_Result@Result[[i]][(1:num_iters)*2]
  
  train_index <- Test_Result@Processed_Data@train_index
  
  test_index <- (1:n)[! (1:n) %in% train_index]
  
  for(j in 1:num_iters){
    
    model_size[i,j] <- length(Index_set[[j]])
    
    cumulative_time[i,j] <- sum(unlist(Iters[1:j]))
    
    model <- glm(Y~., data = data.frame(X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index]),family = binomial())
    
    Train_loss[i,j] <-  -logLik(model)/train_size
    
    newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
    
    y_test <- y[test_index]
    
    Test_Loss[i,j] <-  -(y_test%*%log(predict(model,newdata,type = "response"))+(1-y_test)%*%log(1-predict(model,newdata,type = "response")))/(187-train_size)
    
    bic_value[i,j] <- BIC(model)
  }
}

### plot
layout(matrix(c(1,2,3,4),2))
Test_Result@Methods[[2]] <- "AlphaInvesting"
Test_Result@Methods[[3]] <- "SIS"
Test_Result@Methods[[4]] <- "OSFS"
Test_Result@Methods[[5]] <- "Saola"
op <- par(mfrow = c(2, 2),mar=c(5, 4, 4, 2),mai=0.7*c(2,2,0.5,0.5))  
###

plot(NULL,xlab = "step",ylab ="Train Loss", xlim=c(0,num_iters), ylim =c(min(Train_loss),max(Train_loss)), cex.aixs = 1.5,  cex.lab=1.5)

for(i in 1:num_methods){
  
  lines(Train_loss[i,],lty=i,col = rainbow(7)[i])
  
}
legend("topright",legend = paste0(Test_Result@Methods),cex =1 ,lty =1:num_methods,col=rainbow(7)[1:num_methods])



plot(NULL,ylab ="Test Loss",xlab = "step", xlim=c(0,num_iters), ylim =c(min(Test_Loss),max(Test_Loss)), cex.aixs = 1.5,  cex.lab=1.5)

for(i in 1:num_methods){
  
  lines(Test_Loss[i,],lty=i,col = rainbow(7)[i])
  
}


###

plot(NULL,ylab ="BIC",xlab = "step", xlim=c(0,num_iters), ylim =c(min(bic_value),max(bic_value)), cex.aixs = 1.5,  cex.lab=1.5)

for(i in 1:num_methods){
  
  lines(bic_value[i,],lty=i,col = rainbow(7)[i])
  
}

plot(NULL,ylab ="Time",xlab = "step", xlim=c(0,num_iters), ylim =c(min(cumulative_time),max(cumulative_time)), cex.aixs = 1.5,  cex.lab=1.5)

for(i in 1:num_methods){
  
  lines(cumulative_time[i,],lty=i,col = rainbow(7)[i])
  
}

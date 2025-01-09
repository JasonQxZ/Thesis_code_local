source("~/Library/CloudStorage/OneDrive-UniversityofOttawa/Research/Streaming/algorithm.R")

SMK_data<-load("./SMK_CAN_187_x_y.rds")

train_index <- sample(1:dim(x)[1],170,replace = FALSE)

test_index <- (1:187)[!(1:187) %in% train_index]

Data<-list(Y = y,X = x, coef_true= NULL,subset_true= 1:25)

raw_data<-new("Streaming_Data" , X =Data$X, y =Data$Y , causal_index = Data$subset_true)

processed_data<- setData_(raw_data,train_index = train_index)

Methods = list("BANS","alpha_investing","Iter_SIS","OFS_fisher","Saola_z","BANS_l15","BANS_l20")

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
    
    Train_loss[i,j] <-sum( abs( predict(model,type = "response") - y[train_index]))
    
    newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
    
    Test_Loss[i,j] <- sum( abs( predict(model,newdata,type = "response") - y[test_index]))
    
    bic_value[i,j] <- BIC(model)
  }
}

### plot
layout(matrix(c(1,2,3,4),2))

###

plot(NULL,xlab = "t",ylab ="", xlim=c(0,num_iters), ylim =c(min(Train_loss),max(Train_loss)))

for(i in 1:num_methods){
  
  lines(Train_loss[i,],lty=i,col = rainbow(7)[i])
  
}
legend("topright",legend = paste0(Test_Result@Methods),cex =1 ,lty =1:num_methods,col=rainbow(7)[1:num_methods])

title("Train_Loss")


plot(NULL,ylab ="",xlab = "t", xlim=c(0,num_iters), ylim =c(min(Test_Loss),max(Test_Loss)))

for(i in 1:num_methods){
  
  lines(Test_Loss[i,],lty=i,col = rainbow(7)[i])
  
}

title("Test_Loss")

###

plot(NULL,ylab ="",xlab = "t", xlim=c(0,num_iters), ylim =c(min(bic_value),max(bic_value)))

for(i in 1:num_methods){
  
  lines(bic_value[i,],lty=i,col = rainbow(7)[i])
  
}

title("bic_value")

plot(NULL,ylab ="",xlab = "t", xlim=c(0,num_iters), ylim =c(min(cumulative_time),max(cumulative_time)))

for(i in 1:num_methods){
  
  lines(cumulative_time[i,],lty=i,col = rainbow(7)[i])
  
}

title("time")

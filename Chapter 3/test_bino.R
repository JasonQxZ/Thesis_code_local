source("~/Library/CloudStorage/OneDrive-UniversityofOttawa/Research/Streaming/algorithm.R")

Methods = list("BANS","alpha_investing","Iter_SIS","OFS_fisher","Saola_z","OSFS_FI")

correlation = "ID"

#set.seed(6)

N = 600

train_size = 400

P= 1500

p= 300

data<-list()

for( i in 1:(P/p)){
  
  data[[i]]<- Gen_Data(n = N, p = p, pos_truecoef = 1:5, family = "gaussian"
                       ,effect_truecoef = c(4,-4,3,-7,4),correlation = correlation )
  
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

raw_data<-new("Streaming_Data" , X =Data$X, y =Data$Y , causal_index = Data$subset_true)

processed_data<- setData_(raw_data,train_index = 1:train_size)

a1 <- new("Algorithm", Methods = Methods, Processed_Data = processed_data) 

k = 25

s = 20

Test_Result <- run(a1, shuffle= TRUE , s , k ,family = binomial())

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
    
    online_FDR[i,j] <- 1 - sum(Data$subset_true %in% Test_Result@shuffle_order[Index_set[[j]]])/length(Index_set[[j]])
    
    #model_size[i,j] <- length(Index_set[[j]])
    
    model <- glm(Y~., data = data.frame(X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index]),family = binomial())
    
    Train_loss[i,j] <-sum( abs( predict(model,type = "response") - y[train_index]))
    
    newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
    
    Test_Loss[i,j] <- sum( abs( predict(model,newdata,type = "response") - y[test_index]))
    
    bic_value[i,j] <- BIC(model)
  }
}
### plot
layout(matrix(c(1,2,3,4),2,2))
Test_Result@Methods[[2]] = "AlphaInvesting"
Test_Result@Methods[[3]] = "SIS"
###

plot(NULL,xlab = "t",ylab="s", xlim=c(0,num_iters), ylim =c(0,max(cumulative_time)))
for(i in 1:num_methods){
  
  lines(cumulative_time[i,],lty=i,col = rainbow(7)[i])
  
}
legend("topleft",legend = paste0(Test_Result@Methods),bty = "n",cex =1,lty =1:num_methods,col=rainbow(7)[1:num_methods])
title("Cumulative Time")

###
plot(NULL,xlab = "t", ylab="", xlim=c(0,num_iters), ylim =c(0,max(online_PSR)))

for(i in 1:num_methods){
  
  lines(online_PSR[i,],lty=i,col = rainbow(7)[i])
  
}
causl_features <- unique((1:p)[Test_Result@shuffle_order %in% Data$subset_true] %/% s)

for( i in 1: length(causl_features)){
  
  text( x = causl_features[i] , y = 0, expression('*'), cex = 1.5)
  
}

title("online_PSR")


####
plot(NULL,xlab = "t", xlim=c(0,num_iters), ylim =c(min(Train_loss),max(Train_loss)))
for(i in 1:num_methods){
  
  lines(Train_loss[i,],lty=i,col = rainbow(7)[i])
  
}


title("Train_Loss")


###

plot(NULL,xlab = "t",ylab="", xlim=c(0,num_iters), ylim =c(min(Test_Loss),max(Test_Loss)))

for(i in 1:num_methods){
  
  lines(Test_Loss[i,],lty=i,col = rainbow(7)[i])
  
}

title("Test_Loss")

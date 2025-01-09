source("~/Library/CloudStorage/OneDrive-UniversityofOttawa/Research/Streaming/algorithm.R")

Methods = list("BANS","alpha_investing","Iter_SIS","OFS_fisher","Saola_z","OSFS_FI")

correlation = "CS"

#set.seed(6)

N = 800

train_size = 600

P= 1500

p= 300

data<-list()

for( i in 1:(P/p)){
  
  data[[i]]<- Gen_Data(n = N, p = p, pos_truecoef = 1:5, family = "gaussian"
                       ,effect_truecoef = c(4,-7,3,-7,4),correlation = correlation )
  
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

s = 25

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
    
    online_PSR[i,j] <- sum(Data$subset_true %in%  Test_Result@shuffle_order[Index_set[[j]]])/length(Data$subset_true)
    
    model <- glm(Y~., data = data.frame(X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index]),family = binomial())
    
    Train_loss[i,j] <-  -logLik(model)/train_size
    
    newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
    
    y_test <- y[test_index]
    
    Test_Loss[i,j] <-  -sum(y_test*log(predict(model,newdata,type = "response"))+(1-y_test)*log(1-predict(model,newdata,type = "response")))/(N-train_size)
    
    bic_value[i,j] <- BIC(model)
  }
}

layout(matrix(c(1,2,3,4),2,2))
Test_Result@Methods[[2]] <- "AlphaInvesting"
Test_Result@Methods[[3]] <- "SIS"
Test_Result@Methods[[4]] <- "OSFS"
Test_Result@Methods[[5]] <- "Saola"
###
op <- par(mfrow = c(2, 2),mar=c(5, 4, 4, 2),mai=0.7*c(2,2,0.5,0.5))  
plot(NULL,xlab = "step",ylab="Time", xlim=c(0,num_iters), ylim =c(0,max(cumulative_time)), cex.aixs = 1.5,  cex.lab = 2)
for(i in 1:num_methods){
  
  lines(cumulative_time[i,],lty=i,col = rainbow(7)[i])
  
}
legend("topleft",legend = paste0(Test_Result@Methods),bty = "n",cex = 1.2,lty =1:num_methods,col=rainbow(7)[1:num_methods])

###
plot(NULL,xlab = "step", ylab="PSR", xlim=c(0,num_iters), ylim =c(0,max(online_PSR)), cex.aixs = 1.5,  cex.lab = 2)

for(i in 1:num_methods){
  
  lines(online_PSR[i,],lty=i,col = rainbow(7)[i])
  
}
causl_features <- unique((1:p)[Test_Result@shuffle_order %in% Data$subset_true] %/% s)

for( i in 1: length(causl_features)){
  
  text( x = causl_features[i] , y = 0, expression('*'), cex = 1.5)
  
}


####
plot(NULL,xlab = "step", xlim=c(0,num_iters),ylab = "Train Loss", ylim =c(min(Train_loss),max(Train_loss)), cex.aixs = 1.5,  cex.lab = 2)
for(i in 1:num_methods){
  
  lines(Train_loss[i,],lty=i,col = rainbow(7)[i])
  
}



###

plot(NULL,xlab = "step", xlim=c(0,num_iters),ylab = "Test Loss", ylim =c(min(Test_Loss),max(Test_Loss)), cex.aixs = 1.5,  cex.lab = 2)

for(i in 1:num_methods){
  
  lines(Test_Loss[i,],lty=i,col = rainbow(7)[i])
  
}



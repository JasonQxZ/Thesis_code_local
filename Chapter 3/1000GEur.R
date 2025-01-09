snp <- read.table("1000G_Euro_168000SNPs.dat",   header=TRUE)
snp_y <- as.double(snp[,1])
snp_x <- standardize(matrix(as.double(as.matrix(snp[,4:168798])),nrow = 292))           
snp_y<- as.numeric(snp_y)                
                  
train_size =200
Methods = list("BANS","alpha_investing","Iter_SIS")
Metrics = list("time"  , "trainloss", "bic", "model_size")
raw_data <- new("Streaming_Data" , X = snp_x , y =snp_y , causal_index = 1:5)

processed_data <- setData_(raw_data,train_index = 1:200)

a1 <- new("Algorithm",Methods = Methods , Processed_Data = processed_data) 

N = 292

k = 10

s = 1000          

Test_Result <- run(a1, shuffle= FALSE , s , k)

X <- Test_Result@Processed_Data@X

y <- Test_Result@Processed_Data@y

n <- dim(X)[1]

p <- dim(X)[2]

num_methods <- length(Test_Result@Result)

num_iters <- length(Test_Result@Result[[1]])/2

cumulative_time <- matrix(0,nrow = num_methods, ncol = num_iters)

online_PSR <- matrix(0,nrow = num_methods, ncol = num_iters)

Train_loss <- matrix(0,nrow = num_methods, ncol = num_iters)

Test_Loss <- matrix(0,nrow = num_methods, ncol = num_iters)

model_size <- matrix(0,nrow = num_methods, ncol = num_iters)

bic_value <- matrix(0,nrow = num_methods, ncol = num_iters)

for( i in 1:num_methods){
  
  Iters <- Test_Result@Result[[i]][(1:num_iters)*2-1]
  
  Index_set <-Test_Result@Result[[i]][(1:1:num_iters)*2]
  
  train_index <- Test_Result@Processed_Data@train_index
  
  test_index <- (1:n)[! (1:n) %in% train_index]
  
  for(j in 1:num_iters){
    
    cumulative_time[i,j] <- sum(unlist(Iters[1:j]))
    
    model_size[i,j] <- length(Index_set[[j]])
    
    online_PSR[i,j] <-  sum(Data$subset_true %in%  Test_Result@shuffle_order[Index_set[[j]]])/length(Data$subset_true)
    
    model <- glm(Y~., data = data.frame(X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index]))
    
    Train_loss[i,j] <-sum( ( predict(model,type = "response") - y[train_index])^2)/train_size
    
    newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
    
    Test_Loss[i,j] <- sum( ( predict(model,newdata,type = "response") - y[test_index])^2)/(N-train_size)
    
    bic_value[i,j] <- BIC(model)
    
  }
}
### plot
op <- par(mfrow = c(2, 2),mar=c(5, 4, 4, 2),mai=0.7*c(2,2,0.5,0.5))  
plot(NULL,xlab = "step",ylab="Time", xlim=c(0,num_iters), ylim =c(0,max(cumulative_time)),   cex.lab = 2)
for(i in 1:num_methods){
  
  lines(cumulative_time[i,],lty=i,col = rainbow(7)[i])
  
}
legend("topleft",legend = paste0(Test_Result@Methods),bty = "n",cex = 1.2,lty =1:num_methods,col=rainbow(7)[1:num_methods])

###
plot(NULL,xlab = "step", ylab="PSR", xlim=c(0,num_iters), ylim =c(0,max(online_PSR)),   cex.lab = 2)

for(i in 1:num_methods){
  
  lines(online_PSR[i,],lty=i,col = rainbow(7)[i])
  
}
causl_features <- unique((1:p)[Test_Result@shuffle_order %in% Data$subset_true] %/% s)

for( i in 1: length(causl_features)){
  
  text( x = causl_features[i] , y = 0, expression('*'), cex = 1.5)
  
}


####
plot(NULL,xlab = "step", xlim=c(0,num_iters),ylab = "Train Loss", ylim =c(min(Train_loss),max(Train_loss)), cex.lab = 2)
for(i in 1:num_methods){
  
  lines(Train_loss[i,],lty=i,col = rainbow(7)[i])
  
}



###

plot(NULL,xlab = "step", xlim=c(0,num_iters),ylab = "Test Loss", ylim =c(min(Test_Loss),max(Test_Loss)),  cex.lab = 2)

for(i in 1:num_methods){
  
  lines(Test_Loss[i,],lty=i,col = rainbow(7)[i])
  
}



source("~/Library/CloudStorage/OneDrive-UniversityofOttawa/Research/Streaming/algorithm.R")
library(stringr)
Methods = list("BANS","alpha_investing","Iter_SIS","OFS_fisher","Saola_z","OSFS_FI")
Metrics = list("time" , "psr" , "trainloss", "testloss", "bic", "model_size")
correlation = "CS"
shuffle = TRUE
num_methods = length(Methods)
num_Simu = 100
set.seed(1)

for( i in Metrics){assign(paste0("Table_r_",i),matrix(0, nrow = num_methods, ncol = num_Simu) )}

for( m in 1:num_Simu){
  
  N = 600
  
  train_size = 400
  
  P = 1500
  
  p = 300
  
  data<-list()
  
  k = 25
  
  for( i in 1:(P/p)){
    
    data[[i]]<- Gen_Data(n = N, p = p, pos_truecoef = 1:5, family = "gaussian"
                         ,effect_truecoef = c(4,-5,3,-5,4),correlation = correlation  )
    
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
  
  s = 25
  
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
  
  for( i in 1:num_methods){
    
    Iters <- Test_Result@Result[[i]][(1:num_iters)*2-1]
    
    Index_set <-Test_Result@Result[[i]][(1:1:num_iters)*2]
    
    train_index <- Test_Result@Processed_Data@train_index
    
    test_index <- (1:n)[! (1:n) %in% train_index]
    
    for(j in 1:num_iters){
      
      cumulative_time[i,j] <- sum(unlist(Iters[1:j]))
      
      online_PSR[i,j] <- sum(Data$subset_true %in%  Test_Result@shuffle_order[Index_set[[j]]])/length(Data$subset_true)
      
      model_size[i,j] <- length(Index_set[[j]])
      
      model <- glm(Y~., data = data.frame(X = X[train_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[train_index]))
      
      Train_loss[i,j] <-sum( ( predict(model,type = "response") - y[train_index])^2)/train_size
      
      newdata <- data.frame( X = X[test_index,Test_Result@shuffle_order[Index_set[[j]]]], Y = y[test_index])
      
      Test_Loss[i,j] <- sum( ( predict(model,newdata,type = "response") - y[test_index])^2)/(N-train_size)
      
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
Test_Result@Methods[[2]] <- "AlphaInvesting"
Test_Result@Methods[[3]] <- "SIS"
Test_Result@Methods[[4]] <- "OSFS"
Test_Result@Methods[[5]] <- "Saola"

names(Result) <- Result_name

PSR_DIFF <- matrix(0,ncol = 5,nrow = num_Simu)
print(Result)

for(i in 1:5){PSR_DIFF[,i] <- Table_r_psr[1,]-Table_r_psr[i+1,]}
colnames(PSR_DIFF) = unlist(Methods)[2:6]
TEST_DIFF <- matrix(0,ncol = 5,nrow = num_Simu)
for(i in 1:5){TEST_DIFF[,i] <- -Table_r_testloss[1,]+Table_r_testloss[i+1,]}

colnames(TEST_DIFF) = unlist(Methods)[2:6]
par(mar=c(6,5,4,2))

colnames(PSR_DIFF) <- c("AlphaInv.","SIS" ,"OSFS","Saola","FI")
colnames(TEST_DIFF)<-c("AlphaInv.","SIS" ,"OSFS","Saola","FI")

boxplot(PSR_DIFF,cex.axis = 1.5,cex.lab = 2, ylab ="Differece in PSR", boxfill = NA, border = NA)
boxplot(PSR_DIFF, xaxt = "n", add = TRUE, boxfill="red", cex.axis = 1.5,
        boxwex=0.25, at = 1:ncol(PSR_DIFF) - 0.15) 
abline(h=0,lty=2)

#shift these left by -0.15
boxplot(TEST_DIFF,cex.axis = 1.5,cex.lab=2, ylab ="Difference in Test Loss", boxfill = NA, border = NA,cex.axis = 1.2)
boxplot(TEST_DIFF, xaxt = "n", add = T, boxfill="blue",cex.axis = 1.5,
        boxwex=0.25, at = 1:ncol(TEST_DIFF) + 0.15) 
abline(h=0,lty=2)

save(Result,Table_r_bic,Table_r_psr,Table_r_testloss,Table_r_trainloss,Table_r_model_size, file = paste(correlation,num_Simu,"_N",train_size,"_P",P))

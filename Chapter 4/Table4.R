library(stringr)
library(SMLE)
library(SIS)
library(abess)
library(glmnet)
Root <- "C://Users/Gamer PC/Desktop/SMLE_bess/"
set.seed(1)
source(paste0(Root,"/tool.R"))
source(paste0(Root,"/Simu_algo.R"))
source(paste0(Root,"/Sets_define.R"))
source(paste0(Root,"/Sbess.R"))
source(paste0(Root,"/aSBess.R"))

Methods <- list("aSbess","abess_select","SMLE_select","SIS_s","lasso_cv")

Metrics <- list( "TPR" ,"FDR", "SLE", "time","model_size","test_error")

Num_methods <- length(Methods)

Num_simu <- 1000

Results <- list()

N = 500

P = 7000

family = binomial()

causal_feature_pos <- c(101,103,105,107,109,111)

causal_features_coef = 1.5*c(1,-1,1,-1,1,-1)

for( i in Metrics ){assign(paste0("Table_",i),matrix(0,nrow=Num_methods,ncol = Num_simu))}

for(j in 1:Num_simu){
  
  N_test <- floor(0.2*N)
    
  Data <- list()
    
  ar1_cor <- function(P, rho) {
      
    exponent <- abs(matrix(1:P - 1, nrow = P, ncol = P, byrow = TRUE) - (1:P - 1))
      
    rho^exponent
      
  }
    
  rho = 0.8
    
  Data$X <- mvnfast::rmvn(n = N+N_test, mu = rep(0,P), ar1_cor(P,rho))
    
  Beta <- rep(0,P)
    
  Beta[causal_feature_pos] <- causal_features_coef
    
  theta <- Data$X %*% Beta 
    
  prob <- exp(theta)/(1+exp(theta)) 
    
  Data$Y  <- rbinom(N+N_test, 1, prob = prob)
    
  Data_test_X <- Data$X[(N+1):(N+N_test),]
  
  Data_test_y <- Data$Y[(N+1):(N+N_test)]
  
  Data$X <- Data$X[1:N,]
  
  Data$Y <- Data$Y[1:N]
  
  for(i in 1:length(Methods)){
      
    time1 <- proc.time()
      
    id <- Run_method(Methods[[i]], k = NULL , family, Data)
     
    newdata <- data.frame( X = Data_test_X[,id], Y = Data_test_y)
    
    model <- glm(Y~., data = data.frame(X = Data$X[,id], Y = Data$Y),family = family)
    
    theta <- predict(model,newdata,type = "response")
    
    Table_test_error[i,j] <- logit_loss(Data_test_y,theta)
      
    Table_time[i,j] <- (proc.time()-time1)[3]
      
    Table_TPR[i,j] <- sum(causal_feature_pos %in% id)/length(causal_feature_pos)
      
    non_causal <-  sub_off(1:P,causal_feature_pos)
      
    non_select <- sub_off(1:P,id)
      
    if(length(id)==0){Table_FDR[i,j] = 0}else{Table_FDR[i,j] <- sum(non_causal %in% id)/length(id)}
    
    Table_SLE[i,j] <- abs(length(causal_feature_pos)-length(id))
      
    Table_model_size[i,j] <- length(id)
  }
}
  
Result_name <- str_subset(ls(),"Table_")
Result <- lapply(Result_name,function(i){
  round(rowSums(get(i))/Num_simu,2)
})

names(Result) <- str_replace(Result_name,"Table_", "")
  
df.Result <- data.frame(Result,row.names = Methods)
  
print(df.Result)

save(df.Result, file = "Table4.RData")

#Table 5

Table_FDR_transposed <- t(Table_FDR)
Table_SLE_transposed <- t(Table_TPR)
Table_test_error_transposed <- t(Table_test_error)
Table_time_transposed <- t(Table_time)


metrics_list_transposed <- list(
  FDR = Table_FDR_transposed,
  TPR = Table_SLE_transposed,
  Test_Error = Table_test_error_transposed,
  Time = Table_time_transposed
)

# Method names to be used in the boxplots
method_names <- c("aISSE", "Abess", "SMLE-EBIC", "SIS-S", "Lasso-CV")

# Create boxplots for each transposed metric
par(mfrow = c(2, 2)) # To plot all 4 boxplots in a 2x2 grid

# Boxplot for SLE (transposed) with method names
boxplot(metrics_list_transposed$TPR, 
        main = "TPR",
        xlab = "Methods",
        ylab = "TPR",
        col = "lightgreen",
        names = method_names,   # Adding method names as labels
        outline = FALSE)

# Boxplot for FDR (transposed) with method names
boxplot(metrics_list_transposed$FDR, 
        main = "FDR",
        xlab = "Methods",
        ylab = "FDR",
        col = "lightblue",
        names = method_names,   # Adding method names as labels
        outline = FALSE)



# Boxplot for Test Error (transposed) with method names
boxplot(metrics_list_transposed$Test_Error, 
        main = "Test Error",
        xlab = "Methods",
        ylab = "Test Error",
        col = "lightcoral",
        names = method_names,   # Adding method names as labels
        outline = FALSE)

# Boxplot for Time (transposed) with method names
boxplot(metrics_list_transposed$Time, 
        main = "Time",
        xlab = "Methods",
        ylab = "Time (seconds)",
        col = "lightyellow",
        names = method_names,   # Adding method names as labels
        outline = FALSE)

# Reset the plotting layout to the default
par(mfrow = c(1, 1))




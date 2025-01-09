# install.packages(c("MASS","stringr"))

library(stringr)
library(SMLE)
library(SIS)
library(abess)

Root <- "~/Library/CloudStorage/OneDrive-UniversityofOttawa/SMLE_bess/SMLE_bess"

source(paste0(Root,"/tool.R"))
source(paste0(Root,"/Simu_algo.R"))
source(paste0(Root,"/Sets_define.R"))
source(paste0(Root,"/Sbess.R"))

Methods <- list( "SIS", "Lasso","abess","SMLE_Lasso","Sbess","Sbess_S1")

#Methods <- list("SMLE_zero","SMLE_Lasso","Sbess_S3","Sbess_S4")

Metrics <- list("time", "psr" ,"model_size", "ssr","test_error")

Num_methods <- length(Methods)

Num_simu <- 10

for( i in Metrics ){assign(paste0("Table_",i),matrix(0,nrow=Num_methods,ncol = Num_simu))}

family = binomial()

for(j in 1:Num_simu){
  
  N = 200
  
  P = 3000
  
  p = 1000
  
  causal_feature_pos <- c()
  
  data<-list()
  
  k = 150  
  
  N_test <- floor(0.2*N)
  
  for( i in 1:(P/p)){
    
    data[[i]]<- Gen_Data(n = N+N_test , p = p, pos_truecoef = 101:106, family = "gaussian"
                         ,effect_truecoef = c(2,-2,3,-3,-3,3),correlation = "CS"  )
  }
  
  Data<-list(Y = rep(0,N+N_test),X= NULL,coef_true= NULL,subset_true= NULL)
  
  for(i in 1:(P/p)){
    Data$Y <- Data$Y + data[[i]]$Y
    Data$X <- cbind(Data$X,data[[i]]$X)
    causal_feature_pos <- c(causal_feature_pos,data[[i]]$subset_true+p*(i-1))
  }
  
  Data$Y <- rbinom(n = N+N_test, size = 1 ,prob = exp(Data$Y )/(1+exp(Data$Y )))
  
  Data_test_X <- Data$X[(N+1):(N+N_test),]
  Data_test_y <- Data$Y[(N+1):(N+N_test)]
  
  
  Data$X <- Data$X[1:N,]
  Data$Y <- Data$Y[1:N]
  
  for(i in 1:length(Methods)){
    
    time1 <- proc.time()
    
    id <- Run_method(Methods[[i]], k ,family ,Data)
    
    model <- glm(Y~., data = data.frame(X = Data$X[,id], Y = Data$Y),family = family)
    
    newdata <- data.frame( X = Data_test_X[,id], Y = Data_test_y)
    
    Table_test_error[i,j] <- sum(abs(Data_test_y-predict(model,newdata,type = "response")))/N_test
    
    Table_time[i,j] <- (proc.time()-time1)[3]
    
    Table_psr[i,j] <- sum(causal_feature_pos %in% id)/length(causal_feature_pos)
    
    Table_model_size[i,j] <- length(id)
    
    Table_ssr[i,j] <- floor(sum(causal_feature_pos %in% id)/length(causal_feature_pos))
    
  }
  
}

Result_name <- str_subset(ls(),"Table_")

Result <- lapply(Result_name,function(i){
  round(rowSums(get(i))/Num_simu,2)
  
})
names(Result) <- str_replace(Result_name,"Table_", "")

df.Result <- data.frame(Result,row.names = Methods)

filename <- paste0()
write.table(df.Result, file = "model_comparison.txt", sep = "\t")
print(df.Result)

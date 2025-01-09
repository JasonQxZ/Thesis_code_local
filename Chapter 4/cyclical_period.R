# install.packages(c("MASS","stringr"))
set.seed(123)
library(stringr)
library(SMLE)
library(SIS)
library(abess)
Root <- "~/Desktop/SMLE_bess/"
source(paste0(Root,"/tool.R"))
source(paste0(Root,"/Simu_algo.R"))
source(paste0(Root,"/Sets_define.R"))
source(paste0(Root,"/Sbess.R"))

Methods <- list( "SIS", "Lasso","abess","SMLE","Sbess")

#Methods <- list("SMLE_zero","SMLE_Lasso","Sbess_S3","Sbess_S4")

Metrics <- list("time", "psr" ,"model_size", "ssr","test_error")

Num_methods <- length(Methods)

Data<-list()

Num_simu <- 10

k = 15

for( i in Metrics ){assign(paste0("Table_",i),matrix(0,nrow=Num_methods,ncol = Num_simu))}

family = gaussian()

for(j in 1:Num_simu){
  
  N <- 200 # number of observations
  
  N_test <- floor(0.2*N)
  
  P <- 2000 # number of features
  
  cyclical_period <- 1000 # Define cyclical period for correlation
  
  amplitude <- 0.8 # Increase for stronger correlation
  
  # Generate cyclical correlation matrix
  
  cor_matrix <- matrix(nrow = P, ncol = P)
  
  for (m in 1:P) {
    
    for (n in 1:P) {
      
      # Cyclical correlation pattern with increased amplitude
      
      cor_matrix[m, n] <- amplitude * cos(2 * pi * abs(m - n) / cyclical_period)
    }
  }

  diag(cor_matrix) <- 1
  
  X <- mvrnorm(N+N_test, mu = rep(0, P), Sigma = cor_matrix)
  
  beta <- rep(0,P)# Coefficients for each feature, others are 0
  
  causal_feature_coef<-  0.5*c(1,-1.6,1,-1.2,1.2)
      
  causal_feature_pos <- c(1, 500, 1000, 1500, 2000)
      
  beta[causal_feature_pos]<- causal_feature_coef
  
  # Calculate the linear predictor
  
  linear_predictor <- X %*% beta
  
  sigma <- 0.5 # Standard deviation of the noise
  
  y <- linear_predictor + rnorm(N+N_test, mean = 0, sd = sigma)
  
  Data$X <- X
  
  Data$Y <- y
  
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

save(df.Result, file = "Gaussain_cyclical.RData")

print(df.Result)

library(stringr)
library(SMLE)
library(SIS)
library(abess)
library(parallel)
library(mvnfast)

Root <- "~/Desktop/SMLE_bess"

source(paste0(Root,"/tool.R"))
source(paste0(Root,"/Simu_algo.R"))
source(paste0(Root,"/Sets_define.R"))
source(paste0(Root,"/Sbess.R"))

Methods <- list("SIS","ISIS","Lasso", "abess", "SMLE", "Sbess")

Metrics <- list("time", "psr" ,"model_size", "ssr","test_error")

Num_methods <- length(Methods)

Num_simu <- 100

Results <- list()

N = 120

P = 2000

k = 20

cyclical_period <- 2000 # Define cyclical period for correlation

amplitude <- 0.9 # Increase for stronger correlation

causal_feature_coef<-  c(2,-2,2,-2,1)

causal_feature_pos <- c(1,500, 1000,1500,2000)

for( i in Metrics ){assign(paste0("Table_",i),matrix(0, nrow=Num_methods,ncol = Num_simu))}

Gen_Data_Binomial_Simulation <- function(idx , N, P , k , causal_feature_pos,causal_feature_coef,Num_simu, Methods,Metrics){
  
  # Generate cyclical correlation matrix
  
  for(j in 1:Num_simu){
    
    Data<-list()
    
    sigma <- c(0.5, 1)
    
    N_test <- floor(0.2 * N)
    
    cor_matrix <- matrix(nrow = P, ncol = P)
    
    for (m in 1:P) {
      
      for (n in 1:P) {
        
        # Cyclical correlation pattern with increased amplitude
        
        cor_matrix[m, n] <- amplitude * cos(2 * pi * abs(m - n) / cyclical_period)
        
      }
      
    }
    
    diag(cor_matrix) <- 1
    
    X <- mvrnorm(N+N_test, mu = rep(0, P), Sigma = cor_matrix)
    
    beta <- rep(0 , P)# Coefficients for each feature, others are 0
    
    beta[causal_feature_pos]<- causal_feature_coef
    
    # Calculate the linear predictor
    
    Data$X <- X
    
    linear_predictor <- X %*% beta
    
    theta <- linear_predictor + rnorm(N+N_test, mean = 0, sd = sigma[idx])
    
    prob <- exp(theta)/(1+exp(theta))
    
    Data$Y  <- rbinom(N+N_test, 1,prob = prob)
    
    Data_test_X <- Data$X[(N + 1):(N + N_test),]
    
    Data_test_y <- Data$Y[(N + 1):(N + N_test)]
    
    Data$X <- Data$X[1:N,]
    
    Data$Y <- Data$Y[1:N]
    
    for(i in 1:length(Methods)){
      
      time1 <- proc.time()
      
      id <- Run_method(Methods[[i]], k ,binomial() ,Data)
      
      model <- glm(Y~., data = data.frame(X = Data$X[,id], Y = Data$Y),family = binomial())
      
      newdata <- data.frame( X = Data_test_X[,id], Y = Data_test_y)
      
      theta <- predict(model,newdata,type = "response")
      
      Table_test_error[i,j] <- logit_loss(Data_test_y,theta)
      
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
  
  df.Result 
  
}

simulation_times <- 1:2

# Parallel execution

Results <- lapply(simulation_times, function(x) Gen_Data_Simulation(x,N, P, k, causal_feature_pos, causal_feature_coef, Num_simu, Methods, Metrics))

names(Results) <- paste("Simulation", simulation_times)

print(Results)

save(Results, file = "Table1_S2.RData")

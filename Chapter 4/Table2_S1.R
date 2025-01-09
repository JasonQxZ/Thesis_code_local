library(stringr)
library(SMLE)
library(SIS)
library(abess)
library(parallel)

Root <- "C:/Users/Gamer PC/Desktop/SMLE_bess"

source(paste0(Root,"/tool.R"))
source(paste0(Root,"/Simu_algo.R"))
source(paste0(Root,"/Sets_define.R"))
source(paste0(Root,"/Sbess.R"))

Methods <- list( "SIS", "Lasso","abess","SMLE","Sbess")

Metrics <- list("time", "psr" ,"model_size", "ssr","test_error")

Num_methods <- length(Methods)

Num_simu <- 100

Results <- list()

N = 200

P = 4000

k = 20

family = gaussian()

causal_feature_pos <- c(101,103,105,107,109,111)

causal_features_coef = 2*c(1,-1,1,-1,1,-1)

Gen_Data_Simulation <- function(idx , N, P , k , family, causal_feature_pos,causal_features_coef,Num_simu, Methods,Metrics){
  
  N_test <- floor(0.2*N)
  
  for( i in Metrics ){assign(paste0("Table_",i),matrix(0,nrow=Num_methods,ncol = Num_simu))}
  
  for(j in 1:Num_simu){
    
    sigma <- c(0.1, 1 )
    
    Data <- list()
    
    ar1_cor <- function(P, rho) {
      
      exponent <- abs(matrix(1:P - 1, nrow = P, ncol = P, byrow = TRUE) - (1:P - 1))
      
      rho^exponent
      
    }
    
    rho = 0.8
    
    Data$X <- mvnfast::rmvn(n = N+N_test, mu = rep(0,P), ar1_cor(P,rho))
    
    Beta <- rep(0,P)
    
    Beta[causal_feature_pos] <- causal_features_coef
    
    Data$Y <- Data$X %*% Beta + rnorm( N+N_test , 0 ,sigma[idx])
    
    Data_test_X <- Data$X[(N+1):(N+N_test),]
    Data_test_y <- Data$Y[(N+1):(N+N_test)]
    
    Data$X <- Data$X[1:N,]
    Data$Y <- Data$Y[1:N]
    
    for(i in 1:length(Methods)){
      
      time1 <- proc.time()
      
      test_Data <- Data
      
      id <- Run_method(Methods[[i]], k ,family ,Data)
      
      model <- glm(Y~., data = data.frame(X = Data$X[,id], Y = Data$Y),family = family)
      
      newdata <- data.frame( X = Data_test_X[,id], Y = Data_test_y)
      
      Table_test_error[i,j] <- crossprod(Data_test_y-predict(model,newdata,type = "response"))/N_test
      
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

simulation_times <- 1:2 # Or any range/vector you need


# Parallel execution
Results <- lapply(simulation_times, function(x) Gen_Data_Simulation(x,N, P, k, 
                                                                    
  family, causal_feature_pos, causal_features_coef, Num_simu, Methods, Metrics))

# Optional: Name the list elements based on simulation times for clarity
names(Results) <- paste("Simulation", simulation_times)

print(Results)

save(Results, file = "Table2_S1.RData")

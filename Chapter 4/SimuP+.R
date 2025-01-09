# install.packages(c("MASS","stringr"))

library(stringr)
library(SMLE)
library(SIS)
library(abess)
library(parallel)

Root <- "~/Desktop/SMLE_bess/"

source(paste0(Root,"/tool.R"))
source(paste0(Root,"/Simu_algo.R"))
source(paste0(Root,"/Sets_define.R"))
source(paste0(Root,"/Sbess.R"))

Methods <- list( "SIS", "Lasso","abess","SMLE_Lasso","Sbess")

Metrics <- list("time", "psr" ,"model_size", "ssr","test_error")

Num_methods <- length(Methods)

Num_simu <- 50

Results <- list()

N = 300

k = 20

family = binomial()

causal_feature_pos <- c(101,103,105,107,109,111)

causal_features_coef = 2.5*c(2,-2,3,-3,-3,3)

Gen_Data_Simulation <- function(simulation_time, N , k , family, causal_feature_pos,causal_features_coef,Num_simu, Methods,Metrics){
  
  N = 300
  
  N_test <- floor(0.2*N)
  
  P <- 1000+ 500*simulation_time
  
  for( i in Metrics ){assign(paste0("Table_",i),matrix(0,nrow=Num_methods,ncol = Num_simu))}
  
  for(j in 1:Num_simu){
    
    Data<- Gen_Data(n = N+N_test, p = P, family = family$family,  correlation = "AR", rho = 0.8,
                    pos_truecoef = causal_feature_pos, effect_truecoef = causal_features_coef)
    
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
  
  df.Result 
  
}

# Define the range of simulation times
simulation_times <- 1:10  # Or any range/vector you need

# Sequential execution using lapply
Results <- lapply(simulation_times, function(x) {
  Gen_Data_Simulation(x, N, k, family, causal_feature_pos, causal_features_coef, Num_simu, Methods, Metrics)
})

# Optional: Name the list elements based on simulation times for clarity
names(Results) <- paste("Simulation", simulation_times)

# Save the results to a file
save(Results, file = "2.23_increasingP_4in1plots.RData")

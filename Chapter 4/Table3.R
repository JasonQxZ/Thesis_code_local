library(stringr)
library(SMLE)
library(SIS)
library(abess)
library(parallel)

Root <- '/home/jasonz/scratch/Sbess/SMLE_bess'

source(paste0(Root,"/tool.R"))
source(paste0(Root,"/Simu_algo.R"))
source(paste0(Root,"/Sets_define.R"))
source(paste0(Root,"/Sbess.R"))
source(paste0(Root,"/aSbess.R"))

Methods <- list("SIS", "ISIS", "Lasso", "abess", "SMLE", "Sbess")

Metrics <- list("time", "psr", "model_size", "ssr", "test_error")

Num_methods <- length(Methods)

Num_simu <- 1

Results <- list()

for( i in Metrics ){assign(paste0("Table_",i),matrix(0,nrow=Num_methods,ncol = Num_simu))}

Gen_Data_Gaussian_Simulation <- function(idx , N, P , k , causal_feature_pos,causal_features_coef,Num_simu, Methods,Metrics){
  
  for(j in 1:Num_simu){
    
    sigma <- c(0.2, 1)

    p = P/2
    
    causal_feature_pos <- c()
    
    data<-list()
    
    N_test <- floor(0.2*N)
    
    for( i in 1:(P/p)){
      
      data[[i]]<- Gen_Data(n = N+N_test , p = p, pos_truecoef = 100+seq(1,11,by=2), family = "gaussian"
                           ,effect_truecoef = 2*c(1,-1,1,-1,1,-1),correlation = "AR",rho = 0.7+0.1*i ,sigma = 0)
    }
    
    Data<-list(Y = rep(0,N+N_test),X= NULL)
    
    for(i in 1:(P/p)){
      Data$Y <- Data$Y + data[[i]]$Y
      Data$X <- cbind(Data$X,data[[i]]$X)
      causal_feature_pos <- c(causal_feature_pos,data[[i]]$subset_true+p*(i-1))
    }
    
    Data$Y <- Data$Y + rnorm( N+N_test , 0 ,sigma[idx])
    
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

Gen_Data_Binomial_Simulation <- function(idx , N, P , k , causal_feature_pos,causal_features_coef,Num_simu, Methods,Metrics){
  
  for(j in 1:Num_simu){
    
    sigma <- c(0.01, 0.5)
    
    p = P/2
    
    causal_feature_pos <- c()
    
    data<-list()
    
    N_test <- floor(0.2*N)
    
    for( i in 1:(P/p)){
      
      data[[i]]<- Gen_Data(n = N+N_test , p = p, pos_truecoef = 100+seq(1,11,by=2), family = "gaussian"
                           ,effect_truecoef = 4*c(1,-1,1,-1,1,-1),correlation = "AR",rho = 0.7+0.1*i ,sigma=0)
    }
    
    Data<-list(Y = rep(0,N+N_test),X= NULL)
    
    for(i in 1:(P/p)){
      Data$Y <- Data$Y + data[[i]]$Y
      
      Data$X <- cbind(Data$X,data[[i]]$X)
      
      causal_feature_pos <- c(causal_feature_pos,data[[i]]$subset_true+p*(i-1))
    }
    
    Data$Y <- Data$Y + rnorm( N+N_test , 0 ,sigma[idx])
    
    Data$Y <- rbinom(n = N+N_test, size = 1 ,prob = round(exp(Data$Y )/(1+exp(Data$Y )),3))
    
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

Gen_Data_Poisson_Simulation <- function(idx , N, P , k , causal_feature_pos,causal_features_coef,Num_simu, Methods,Metrics){
  
  for(j in 1:Num_simu){
    
    sigma <- c(0, 0.03)

    p = P/2
    
    causal_feature_pos <- c()
    
    data<-list()
    
    N_test <- floor(0.2*N)
    
    for( i in 1:(P/p)){
      
      data[[i]]<- Gen_Data(n = N+N_test , p = p, pos_truecoef = 100+seq(1,11,by=2), family = "gaussian"
                           ,effect_truecoef = c(1,-1,1,-1,1,-1),correlation = "AR",rho = 0.7+0.1*i ,sigma=0)
    }
    
    Data<-list(Y = rep(0,N+N_test),X= NULL)
    
    for(i in 1:(P/p)){
      Data$Y <- Data$Y + data[[i]]$Y
      Data$X <- cbind(Data$X,data[[i]]$X)
      causal_feature_pos <- c(causal_feature_pos,data[[i]]$subset_true+p*(i-1))
    }
    
    Data$Y <- Data$Y +  rnorm( N+N_test , 0 ,sigma[idx])
    
    Data$Y <- rpois(n = N+N_test ,lambda  = exp(Data$Y ))
    
    Data_test_X <- Data$X[(N+1):(N+N_test),]
    Data_test_y <- Data$Y[(N+1):(N+N_test)]
    
    
    Data$X <- Data$X[1:N,]
    Data$Y <- Data$Y[1:N]
    
    for(i in 1:length(Methods)){
      
      time1 <- proc.time()
      
      test_Data <- Data
      
      id <- Run_method(Methods[[i]], k ,poisson() ,Data)
      
      model <- glm(Y~., data = data.frame(X = Data$X[,id], Y = Data$Y),family = poisson() )
      
      newdata <- data.frame( X = Data_test_X[,id], Y = Data_test_y)
      
      lambda <- predict(model,newdata,type = "response")
      
      Table_test_error[i,j] <-  poisson_log_likelihood_loss(Data_test_y,lambda)/N_test
      
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
Results_G <- lapply(simulation_times, function(x) Gen_Data_Gaussian_Simulation(x,600, 8000, 20, causal_feature_pos, 2*c(1,-1,1,-1,1,-1), Num_simu, Methods, Metrics))

Results_B <- lapply(simulation_times, function(x) Gen_Data_Binomial_Simulation(x,900, 4000, 20, causal_feature_pos,  4*c(1,-1,1,-1,1,-1), Num_simu, Methods, Metrics))

Results_P <- lapply(simulation_times, function(x) Gen_Data_Poisson_Simulation(x,650, 8000, 40, causal_feature_pos,  c(1,-1,1,-1,1,-1), Num_simu, Methods, Metrics))

# Optional: Name the list elements based on simulation times for clarity
names(Results) <- paste("Simulation", simulation_times)
print(Results)
save(Results, file = "Table3.RData")

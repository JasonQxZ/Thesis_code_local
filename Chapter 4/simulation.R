library(SIS)
library(glmnet)
library(abess)
NumIter <- 30
num_casual_featrues = 25
num_methods = 6
###########Linear Regression ####################################

Lasso_avg_length<-rep(1 : NumIter)
PRR<-matrix(0, nrow = num_methods , ncol = NumIter)
SSR<-matrix(0, nrow = num_methods , ncol = NumIter)
TIME<- matrix(0, nrow = num_methods , ncol = NumIter)
time<-rep(0, num_methods )  
k = 50
correlation = "AR"
for(j in 1 : NumIter){
  
  N = 1600
  
  train_size = 1200
  
  P = 8000
  
  p = 1600
  
  data<-list()
  
  k = 50
  
  for( i in 1:(P/p)){
    
    data[[i]]<- Gen_Data(n = N, p = p, pos_truecoef = 1:5, family = "binomial"
                         ,effect_truecoef = i/5*c(4,-5,3,-5,4),correlation = correlation  )
    
  }
  
  Data<-list(Y = rep(0,N),X= NULL,coef_true= NULL,subset_true= NULL)
  
  for(i in 1:5){
    Data$Y <- Data$Y + data[[i]]$Y
    Data$X <- cbind(Data$X,data[[i]]$X)
    Data$coef_true<- c(Data$coef_true,data[[i]]$coef_true)
    Data$subset_true <- c(Data$subset_true,data[[i]]$subset_true+p*(i-1))
    
  }
  pi <- exp(Data$Y) / (1 + exp(Data$Y))
  
  Data$Y  <- rbinom(N, size = 1, prob = pi)
  
  
  ##models fitting
  #SIS
  sis_a <- proc.time()
  SIS <- SIS(x = Data$X,y = Data$Y,family = "binomial",iter = F,nsis = k)
  sis_b <- proc.time()
  time[1] <- (sis_b-sis_a)[3]
  #ISIS
  isis_a <- proc.time()
  SMLE_bess_fit <- Sbess(Data$Y,Data$X, k = k,family = binomial(),lasso_initial = F)
  isis_b <- proc.time()
  time[2] <- (isis_b-isis_a)[3]
  #Lasso
  Lasso_a <- proc.time()
  Lasso_fit <- glmnet(x = Data$X, y = Data$Y, family = "binomial", pmax = k)
  Lasso_b <- proc.time()
  time[3] <- (Lasso_b-Lasso_a)[3]
  #SMLE
  SMLE_a <- proc.time()
  SMLE_package_lasso <- SMLE(X = Data$X,Y = Data$Y,family = "binomial",categorical = FALSE, k = k,intercept  = F)
  SMLE_b <- proc.time()
  time[4] <- (SMLE_b-SMLE_a)[3]
  #SMLE_fast
  fast_a <- proc.time()
  SMLE_package_zero <- SMLE(X = Data$X,Y = Data$Y,family = "binomial",categorical = FALSE,k = k,coef_initial = rep(0,2000),intercept  = F)
  fast_b <- proc.time()
  time[5] <- (fast_b - fast_a)[3]
  abess_a <- proc.time()
  abess_fit <- abess(x = Data$X,y = Data$Y,family = "binomial", support.size = k)
  abess_b <- proc.time()
  time[6]<-(abess_b-abess_a)[3]
  
  Lasso_index <- (1 : dim(Data$X)[2])[Lasso_fit$beta[,dim(Lasso_fit$beta)[2]] != 0]
  Lasso_avg_length[j] <- length(Lasso_index)
  ind <- extract(abess_fit)$support.vars
  ind <- as.numeric(regmatches(ind, gregexpr("[[:digit:]]+", ind)))
  index_set <- list(SIS$ix0, SMLE_bess_fit$IHT_res_out$ID_None_s, Lasso_index, SMLE_package_lasso$ID_retained,
                    SMLE_package_zero$ID_retained,ind)
  
  PRR[,j] <- unlist(lapply(1 : num_methods ,function(i){
    
    sum(Data$subset_true %in% index_set[[i]])
    
  }))
  SSR[,j] <- unlist(lapply(1:num_methods ,function(i){
    
    sum(Data$subset_true %in% index_set[[i]]) == num_casual_featrues
    
  }))
  TIME[,j] <- time
}


# Order: SIS / ISIS / glmnet / SMLE / SMLE_fast 

# Screening Accuracy
PRR_mean <- apply(PRR, 1, mean)/num_casual_featrues
SSR_mean <- apply(SSR, 1, mean)
# Computational Time
Time_mean <- apply(TIME, 1, mean)

result_df <- data.frame(PRR = PRR_mean, SSR = SSR_mean, Time = Time_mean)
row.names(result_df) <- c("SIS","SMLE_Bess","lasso","SMLE_lasso","SMLE_zero","abess")
print(result_df)

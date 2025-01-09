library(SMLE)
source("~/Library/CloudStorage/OneDrive-UniversityofOttawa/Research/Streaming/tool.R")
N= 300

train_size = 200

P= 1500

data<-list()


for( i in 1:5){
  
  data[[i]]<- Gen_Data(n = 300, p = 300, pos_truecoef = 1:5, family = "gaussian"
                   ,effect_truecoef = c(5,2,-3,-6,4),correlation = "AR" , rho=0.8)
  
}

Data<-list(Y = rep(0,N),X= NULL,coef_true= NULL,subset_true= NULL)

for(i in 1:5){
  Data$Y <- Data$Y + data[[i]]$Y
  Data$X <- cbind(Data$X,data[[i]]$X)
  Data$coef_true<- c(Data$coef_true,data[[i]]$coef_true)
  Data$subset_true <- c(Data$subset_true,data[[i]]$subset_true+300*(i-1))
  
}

s = 20

k = 30

yt <- matrix(0, nrow = N - train_size ,ncol =  P)

ytt <- matrix(0, nrow = train_size ,ncol =  P)

error_art <- rep(0,P/s)

error_atr <- rep(0,P/s)

error_Train_Streaming <- rep(0,P/s)

error_str <- rep(0,P/s)

error_srt <- rep(0,P/s)

error_test_Streaming<- rep(0,P/s)

error_real<- rep(0,P/s)

error_realt<- rep(0,P/s)

modelsizea<- rep(0,P/s)

error_size_rate_smle<- rep(0,P/s)

modelsizes<- rep(0,P/s)

beta_s <- rep(0, k + s)

beta_t <- rep(0,k)

ttAI <- rep(0,P/s)

ttSMLE <-rep(0,P/s)

ttstreaming <- rep(0,P/s)



Streaming_data <- function(data){
  
  x <- data$X[1:train_size,]
  
  x_test <- data$X[(train_size+1):N,]
  
  y <- data$Y[1:train_size]
  
  y_test <- data$Y[(train_size+1):N]
  
  df_train<- data.frame(x= x,y =y)
  
  df_test<- data.frame(x=x_test,y=y_test)
  
  # AI 
  
  wealth = 0.5
  
  delta_alpha = 0.5
  
  n <- dim(x)[1]
  
  p <- dim(x)[2]
  
  model <- {}
  
  model_index<-{}
  
  pre_model <- {}
  
  # SMLE 
  
  cand <- NULL
  
  Streaming_cand <- NULL
  
  timeAI <-0
  
  timeSMLE<-0
  
  r_index <- 1:p 
  
  for( i in 1:p ){
    
    j <- sample(r_index,1)
    
    r_index <- r_index[-which(r_index %in% j)]
    
    xi<- matrix(x[,j],ncol =1)
    
    #AI
    t1<-proc.time()
    
    alpha<- wealth/(2*i)
    
    pval<- get_p_val(xi,model,y)
    
    if( pval < alpha){
      
      model <- cbind(model,xi)
      
      model_index <- c(model_index,j)
      
      wealth <- wealth - alpha + delta_alpha
      
    }else{
      
      wealth <- wealth - alpha
      
    }
    timeAI<-timeAI + (proc.time()-t1)[3]
    
    #SMLE
    
    t2<-proc.time()
    
    if( (i-k) %% s == 0 && i >= k){
      
      fit <- SMLE( X = cand , Y = y, k = k)
      
      cand <- cand[,fit$ID_retained]
      
    }else{
      
      cand<- cbind(cand,x[,j])
      
    }
    
    timeSMLE<- timeSMLE + (proc.time()-t2)[3]
    
    #streaming MLE
    
    Streaming_cand <- cbind(Streaming_cand,x[,j])
  
    if( (i-k) %% s == 0 & i >= k){
      
      if(i == k ){
        
        #initial model with size k
        
        X_i <- data.frame(x = Streaming_cand, y = y)
        
        fit<- glm(y~., data=X_i )

        beta_t <- matrix(fit$coefficients[-1],ncol=1)
        
        Residuals <- y - Streaming_cand%*%beta_t
        
        }else{
          
        u = 0.01/sqrt(svd(t(Streaming_cand)%*% Streaming_cand)$d[1])
        
        beta_t_temp <- matrix(c(beta_t,rep(0,s)),ncol=1) + u*(c(crossprod(Streaming_cand[,(1:k)],Residuals),crossprod(Streaming_cand[,(k+1):(k+s)],y)))
        
        beta_t <- Hard( beta_t_temp , k)
        
        retained_index <-which(beta_t_temp==beta_t)
        
        Residuals <- y - Streaming_cand%*%beta_t
        
        beta_t <- beta_t[retained_index]
        
        Streaming_cand <- Streaming_cand[,retained_index]
        
        }
      
    }
    
    
    
    #error 
    
    if( (i-k) %% s == 0 & i > k  & length(model_index)>1){
      
      causal <- data$subset_true < i
      
      yt[,i] <- x_test[,data$subset_true[causal]]%*%data$coef_true[causal]
      
      ytt[,i] <- x[,data$subset_true[causal]]%*%data$coef_true[causal]
      
      A_fit <- glm(y~. , data =  df_train[,model_index])
      
      y_ahat <- predict(A_fit, newdata = df_test)
      
      y_ahat_tr <- predict(A_fit)
      
      error_atr[floor(i/s)] <-sum(abs(y_ahat_tr-y))
      
      error_art[floor(i/s)]<- sum(abs(y_ahat - y_test))
      
      fit_streaming <- glm(y~. , data =  data.frame(x= Streaming_cand,y=y))
      
      error_Train_Streaming[floor(i/s)] <- sum(abs(predict(fit_streaming) - y))
      
      error_test_Streaming[floor(i/s)]<- sum(abs(predict(fit_streaming, newdata = df_test)-y_test))
      
      #PSR_a[i/s] <- sum(data$subset_true %in% model_index) / length(data$subset_true)
      
      smle_model <- find_index(cand,x)
      
      S_fit <-glm(y~. ,data =  cbind(df_train[,smle_model],y))
      
      y_shat <- predict(S_fit , newdata = df_test)
      
      y_shat_tr <- predict(S_fit)
      
      error_str[floor(i/s)] <-sum(abs(y_shat_tr-y))
      
      error_srt[floor(i/s)]<- sum(abs(y_shat - y_test))
      
      #PSR_s[i/s]<- sum(data$subset_true %in% smle_model) / length(data$subset_true)
      
      error_real[floor(i/s)] <- sum(abs(yt[,i] - y_test))
      
      error_realt[floor(i/s)] <- sum(abs(ytt[,i] - y))
      
      modelsizea[floor(i/s)] <- sum(abs(y_ahat - y_test))*length(model_index)
      
      modelsizes[floor(i/s)] <- sum(abs(y_shat - y_test))*length(smle_model)
      
      error_size_rate_smle[floor(i/s)] <- sum(abs(predict(fit_streaming, newdata = df_test)-y_test))*length(smle_model)
      
      ttAI[floor(i/s)] <- timeAI
      
      ttSMLE[floor(i/s)] <- timeSMLE
      
    }
  }
  
  list(error_art, error_srt, error_real, error_size_rate_smle,
       
       modelsizea, modelsizes, ttAI, ttSMLE,error_atr,error_str,error_realt,
       
       error_Train_Streaming, error_test_Streaming
       )
  
}


result<- Streaming_data(Data)

layout(matrix(c(1,2,3,4),2,2))

plot(result[[1]],lty=1, type = 'l', col= 'red', 
     
     ylim = c(0, max(result[[1]],result[[2]],result[[13]])),
     
     xlab = "t", ylab = "Error")

lines(result[[2]],lty=2, col= 'blue')

lines(result[[3]],lty=3, col= 'black')

lines(result[[13]],lty=4, col= 'orange')

title("Prediction Error")

plot(result[[9]],lty=1, type = 'l', col= 'red', 
     
     ylim = c(0, max(result[[9]],result[[10]],result[[12]])),
     
     xlab = "t", ylab = "Error")

lines(result[[10]],lty=2, col= 'blue')

lines(result[[11]],lty=3, col= 'black')

lines(result[[12]],lty=4, col= 'orange')

title("Train Error")

#legend("topright", legend = c("Marginal", "SMLE"), lty =1:2,col = c("red", "blue"))

plot(result[[5]],col="red",xlab = "t", ylab = "ratio",
     
     ylim = c(0,max(result[[5]],result[[6]])) )

points(result[[6]],col="blue")

points(result[[4]],col="orange")
title("Prediction Error/Model Size")

plot(result[[8]],type = 'l', col= 'blue', xlab = "t", ylab = "Time")

lines(result[[7]],lty=2, col= 'red')

title("Cumulative Time")
















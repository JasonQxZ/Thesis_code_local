library(SMLE)

nt <- 10

p= 2000
  
Rep = 10

PSR_SMLE<- c()

PSR_AI<-c()

s=5

for( j in 1:Rep){
  
  d1 <- Gen_Data(n=500, p = p,num_truecoef = nt ,correlation = "AR" , rho =0.8)
  
  cand <- NULL
  for(i in 1:p){
    
    if( i %% s == 0 && i >= 20){
      
      fit <- SMLE( X = cand , Y = d1$Y ,k = 15)
      
      cand <- cand[,fit$ID_retained]
      
    }else{
      
      cand<- cbind(cand,d1$X[,i])
      
    }
  }
  
  
  find_index<-function(A,B){
    
    # find index of A as a submatrix of B
    
    index<-apply(A,2,FUN=function(y){(1:dim(B)[2])[apply(B,2,FUN = function(x){identical(x,y)})]})
    
    index
    
  }
  
  PSR_SMLE[j] <- sum(d1$subset_true %in% find_index(cand,d1$X)) / nt
  
  model_index<-alpha_investing(d1$X,d1$Y)
  
  PSR_AI[j] <- sum(d1$subset_true %in% model_index) / nt
  
}

mean(PSR_SMLE)



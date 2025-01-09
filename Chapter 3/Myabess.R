library(SMLE)
Data<-Gen_Data()
X<-Data$X
y<-as.matrix(Data$Y, nrow = dim(X)[1])

s=10

#ABESS 

Myabess<-function(X, y, s ,tau_s = NULL, k_max = NULL){
  
  # X@matrix-dim-n,p
  # y@matrix-dim-n,1 or @list-length-n.
  # s@integer 
  # s is support size and k_max is the maximum splicing size, and k_max < x.
  # tau@numeric 
  
  n <- dim(X)[1]
  
  p <- dim(X)[2]
  
  if(is.null(tau_s)){tau_s <- 0.01*s*log(p)*log(log(n))/n}
  
  if(is.null(k_max)){k_max <- s}
  
  ## Initialization  ----------
  
  # Initialize A_0
  
  xjTxj <- diag(crossprod(X,X)) 
  
  xjTy<- crossprod(X,y)
  
  marginal_effects <- xjTy/xjTxj
  
  A_index <- order(abs(marginal_effects),decreasing = TRUE)[1:s]
  
  I_index <- sub_off(1:p, A_index) 
  
  # Initialize beta_0 and d_0
  
  beta_0 <- matrix(rep(0,p),ncol = 1) 
  
  beta_0_A <- marginal_effects[A_index]
  
  beta_0[A_index] <- beta_0_A
  
  d_0 <- matrix(rep(0,p),ncol = 1) 
  
  d_0_I <- (crossprod(X,y-X%*%beta_0)/n)[I_index]
  
  d_0[I_index] <- d_0_I
  
  ## Splicing loop  ----------
  
  num_iter <- 0
  
  splicing_object_0<-list(beta_0 = beta_0, d = d_0, A_index = A_index, I_index=I_index)
  
  while( num_iter < 100 ){
    
    splicing_object_1 <- Splicing(X, y, tau_s, k_max, splicing_object_0)
    
    if(identical(splicing_object_0,splicing_object_1)){break}
    
    num_iter <- num_iter + 1
    
    splicing_object_0 <- splicing_object_1
    
  }
  print(num_iter)
  return(splicing_object_1)
  
}


T <- Myabess(X,y,s)


Splicing<-function(X, y, tau_s, k_max, splicing_object){
  
  beta_0 <-splicing_object$beta_0
  # beta@ matrix-dim-p,1
  d <- splicing_object$d
  # d@matrix-dim-p,1
  A_index <- splicing_object$A_index
  # A_index@list-length-s
  I_index <- splicing_object$I_index
  # I_index@list-length-s
  # k_max@integer
  # tau_s@numeric
  
  n <- dim(X)[1]; p <- dim(X)[2]
  
  L_0 <- norm(y-X %*% beta_0,"2")/(2*n)
  
  L <- L_0
  
  xjTxj <- diag(crossprod(X,X)) 
  
  xjTy<- crossprod(X,y)
  
  ## Backward sacrifice
  
  epsilon_j <- (xjTxj*beta_0^2)[A_index]
  
  ## Forward sacrifice
  
  d_j <- (xjTy - xjTxj*beta_0)[I_index]
  
  eta_j <- d_j^2/(xjTxj[I_index])
  
  
  for( k in 1:k_max){
    
    A_k_index <- A_index[order(epsilon_j,decreasing = FALSE)[1:k]]
    
    I_k_index <- I_index[order(eta_j,decreasing = TRUE)[1:k]]
    
    A_k_index_t <- c(sub_off(A_index,A_k_index),I_k_index)
    
    I_k_index_t <- c(sub_off(I_index,I_k_index),A_k_index)
    
    beta_k <- matrix(rep(0,p),ncol = 1) 
    
    beta_k_A <- (xjTy/xjTxj)[A_k_index_t]
    
    beta_k[A_k_index_t] <- beta_k_A
    
    d_k <- matrix(rep(0,p),ncol = 1) 
    
    d_k_I <- (crossprod(X,y-X%*%beta_k)/n)[I_k_index_t]
    
    d_k[I_k_index_t] <- d_k_I
    
    L_k <- norm(y-X%*%beta_k,"2")/(2*n)
    
    if(L_k < L){
      
      beta_1 <- beta_k
      
      d_1 <- d_k
      
      A_index_1 <- A_k_index_t
      
      I_index_1 <- I_k_index_t
      
      L <- L_k
      
    }
    
  }
  # End For
  
  if(L_0 - L < tau_s){
    
    beta_1 <- beta_0
    
    d_1 <- d_k
    
    A_index_1 <- A_index
    
    I_index_1 <- I_index
    
  }
  
  return(list(beta_0 = beta_1, d = d_1, A_index = A_index_1, I_index=I_index_1))
  
}
sub_off<-function(all,sub){
  
  #----------------------------------------------------------------------------#
  # Take elements off a list                                                   #
  #----------------------------------------------------------------------------#
  
  return(all[!(all %in% sub)])
}

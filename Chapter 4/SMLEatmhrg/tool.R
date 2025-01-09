logit_loss <- function(actual_labels, predicted_probabilities) {
  -mean(actual_labels * log(predicted_probabilities) + 
          (1 - actual_labels) * log(1 - predicted_probabilities))
}

ebic <- function(Beta_s = NULL ,Y, X, family , s = NULL,  gamma_ = 0.5){
  
  if(is.null(s)){
    
    p <- dim(X)[2]
    
    k <- sum(Beta_s != 0)
    
    n <- length(Y)
    
    ebic <- -2*lh( Y, X, Beta_s, family = family ) + k * log(n) +2* gamma_ * log(p)
    
    return(ebic)
    
  }else{
    
    p <- dim(X)[2]
    
    Beta_s <- coefficients(glm.fit(x = X[,s],y = Y, family = family))
    
    k <- length(s)
    
    n <- length(Y)
    
    ebic <- -2*lh( Y, X, Beta_s, family = family ) + k * log(n) +2* gamma_ * log(p)
    
    return(ebic)
    
  }
 
}

IHT_w_ucheck <-function(Beta_s, k , Y, X, family,U, U_rate, tol, max_iter,p){

  ###############################################################
    
  # Input : \beta_0 , k ,X , Y and some parameters of convergence 
  
  # Ouput : \beta_s of sparsity k
  
  ###############################################################

  # Set up U and v_s
  
  if( sum(Beta_s != 0) == 0 ){
    
    U_0 <- U/(sqrt(p))
    
    v_s <-  crossprod(X , Y)
    
  }else{
    
    ID_None0 <- which(Beta_s != 0) 
    
    coef_None0 <- as.matrix( Beta_s[ ID_None0 ] , ncol=1) 
    
    Xs_0 <- as.matrix( X[, ID_None0] )
    
    theta_0  <- Xs_0 %*% coef_None0
    
    U_0 <- U/(sqrt(p)*norm(Xs_0,"i")^2)
    
    theta_0 <- make.link(family$link)$linkinv(theta_0)
    
    v_s <- crossprod(X , Y - theta_0)
    
  }

  
  # Check if the sparsity of beta_0 > k
  
  if( sum(Beta_s != 0) > k ){ Beta_s <- Hard( t = Beta_s, k = k)}
  
  i <- 1
  
  LH <- rep(0,max_iter)
  
  FD <- NULL
  
  repeat{
  
    count<- 0
    
    u_s <- U_0
    
    repeat{
      
      Beta_t <- Beta_s + u_s * v_s
      
      # Hard threshold
      
      Beta_t <- Hard(t = Beta_t , k= k)
      
      ucheck<- Ucheck(Y = Y, X = X, beta1 = Beta_s, beta2 = Beta_t, family = family)
      
      if (ucheck >= 0){
        
        break
        
      }else{
        
        u_s <- U_rate * u_s
        
        count <- count + 1
        
      }
      
    }
    
    sindex <- which(Beta_s!= 0)
    
    #  Retained features at step t
    
    tindex <- which(Beta_t!= 0)
    
    # Difference of retained features between two steps
    
    fs<-sum(!(tindex %in% sindex))
    
    FD<-c(FD,fs)
    
    # Difference of retained features between two steps
    
    LH[i] <- lh( Y, X, Beta_t, family = family )
    
    valid_LH_diff<- 0.01 *(LH[2] - LH[1])
    
    if(i>1){
      
      MSE<- sqrt(sum(( Beta_s-Beta_t )^2))
      
      if( MSE < tol ){break}
      else if(i>10){if(sum(tail(FD,10))==0){break}}
      
    } 
    
    Beta_s <- as.vector(Beta_t)
    
    ID_None_s <- sort((1:p)[Beta_s!= 0])
    
    coef_None_s <- as.matrix(Beta_s[ID_None_s],ncol=1)
    
    Xs_s <- X[, ID_None_s]
    
    theta_s  <- make.link(family$link)$linkinv(Xs_s %*% coef_None_s)
    
    v_s <- crossprod(X, Y - theta_s)
    
    u_s  <- U_0
    
    LH_last <- lh( Y, X, Beta_s, family = family )
    
    i <- i + 1
    
  }
  
  IHT_res <- list(rho= v_s, ID_None_s = ID_None_s, Beta_s = coef_None_s,LH = LH_last)
  
  IHT_res
  
}

Ucheck <-function(Y,X, beta1, beta2, family){
  
  #----------------------------------------------------------------------------#
  # This is the ucheck function                                                #
  #----------------------------------------------------------------------------#
  # This function compare the likelihood between beta1 and beta2 (the          #
  # coefficient at the current step and the coefficient for the next step.     #
  #----------------------------------------------------------------------------#
  # Parameters:  
  # Y : The response vector.
  # X : Feature matrix.      
  # beta1 : Coefficient of the current step. 
  # beta2 : Coefficient of the next step. 
  # family: Model assumption.                                                  #
  #----------------------------------------------------------------------------#
  
  ll_1 <- lh(Y=Y,X=X,beta=beta1,family=family)
  
  ll_2 <- lh(Y=Y,X=X,beta=beta2,family=family)
  
  return(ll_2 - ll_1)
  
}


Hard<-function(t=c(0,0), k=-1, lam=1)
{
  
  y<-t
  t<-abs(t)
  
  if(k > 0)
  { lam<-sort(t, decreasing=TRUE)[k] }
  
  y[t<lam]<-0
  return(y)
}




sub_off<-function(all,sub){
  S<-all[!(all %in% sub)]
  return(S)
}


lh<-function(Y, X, beta, family){
  
  #----------------------------------------------------------------------------#
  # This function calculate the log-likelihood                                 #
  #----------------------------------------------------------------------------#
  # This function call the aic() from family object to calculate the likelihood#
  # given the response Y, the data matrix X the coefficient beta and the model #
  # assumption family.                                                         #
  #----------------------------------------------------------------------------#
  # Parameters:  
  # Y : The response vector.
  # X : Feature matrix.      
  # beta : Coefficient input.
  # family: Model assumption.                                                  #
  #----------------------------------------------------------------------------#
  
  #Calculating log likelihood from family object in stats package
  
  linkinv <- family$linkinv
  aic <- family$aic
  dev.resids <- family$dev.resids
  
  y <- Y
  n = length(y)
  ind <- which(beta != 0)
  k <- sum(beta!= 0)
  eta <- X[,ind]%*%as.matrix(beta[ind],nrow = length(ind))
  mu <- linkinv(eta)
  dev <- sum(dev.resids(y, mu, 1))
  return((-1/2)*aic(y, n, mu, 1, dev) + k)
}

Standardize <- function(X) {
  name<-dimnames(X)
  center = colMeans(X)
  X.c = sweep(X, 2, center)
  unit.var = apply(X.c, 2, sd)
  if(sum(unit.var==0)!=0){
    zero_var = (1:dim(X)[2])[unit.var==0]
    X.c<-X.c[,-zero_var]
    unit.var<-unit.var[-zero_var]
  }
  val = sweep(X.c, 2, unit.var, "/")
  Xs=as.matrix(val,dimnames=name)
  return(list(Xs=Xs,X_mean = center, X_sd= unit.var))
}
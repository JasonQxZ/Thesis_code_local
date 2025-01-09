require("DescTools")
library(abess)
library(SIS)
library(parallel)
library(fastglm)

Iter_SMLE <- function(Memory, Gt, y, k,family){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  num_feature_used <- Memory$num_feature_used
  
  Xt<- cbind(St,Gt)
  
  Xt_index  <- c(St_index,(1:dim(Gt)[2])+num_feature_used)
  
  t_start = proc.time()[3]
  
  fit<- SMLE(X = Xt,  Y = y , k = k,family = family)
  
  retained_index <- fit$ID_retained
  
  retained_set <- Xt[,retained_index]
  
  coef <- fit$coef_retained
  
  residual <- y - predict(fit)
  
  St_index <- Xt_index[retained_index]

  time = proc.time()[3] -t_start
  
  return(list(time = time, St_index = St_index, St = retained_set, coefficients =  coef, residual = residual,num_feature_used = num_feature_used + dim(Gt)[2]))
  
}

Iter_SIS <- function(Memory, Gt, y, k,family){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  num_feature_used <- Memory$num_feature_used
  
  Xt<- cbind(St,Gt)
  
  Xt_index  <- c(St_index,(1:dim(Gt)[2])+num_feature_used)
  
  t_start = proc.time()[3]
  
  fit<- SIS(x = Xt,  y = y ,nsis = k,iter = FALSE,family = family$family)
  
  retained_index <- fit$sis.ix0
  
  retained_set <- Xt[,retained_index]
  
  St_index <- Xt_index[retained_index]
  
  time = proc.time()[3] -t_start
  
  return(list(time = time, St_index = St_index, St = retained_set,num_feature_used = num_feature_used + dim(Gt)[2]))
  
}

Iter_abess <- function(Memory, Gt, y, k,family){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  num_feature_used <- Memory$num_feature_used

  Xt<- cbind(St,Gt)
  
  colnames(Xt) <- paste0('x',1:(dim(Xt)[2]))
  
  Xt_index  <- c(St_index,(1:dim(Gt)[2]) + num_feature_used)
  
  t_start = proc.time()[3]
  
  if(length(Xt_index) < k){
    
    return(list(time = proc.time()[3] -t_start, St_index = Xt_index , St = Xt, num_feature_used = num_feature_used + dim(Gt)[2]))
    
  }
  
  abess_fit <- abess(x = Xt,  y = y, support.size = k,family = family$family)
  
  ind_name <- extract(abess_fit)$support.vars
  
  retained_index <- as.numeric(regmatches(ind_name, gregexpr("[[:digit:]]+", ind_name)))
  
  retained_set <- Xt[,retained_index]
  
  St_index <- Xt_index[retained_index]
  
  return(list(time = proc.time()[3] -t_start, St_index = St_index, St = retained_set,num_feature_used = num_feature_used + dim(Gt)[2]))
  
}

Group_alpha_investing <- function(Memory, Gt, y){
  
  St <- Memory$St
  
  wealth <-  Memory$wealth
  
  alpha <- Memory$alpha
  
  St_index <- Memory$St_index
  
  num_feature_used <- Memory$num_feature_used
  
  t_start = proc.time()[3]
  
  for( i in 1:dim(Gt)[2] ){
    
    xi<- Gt[,i]
    
    alpha<- wealth/(2*(i+num_feature_used))
    
    # get p value from git package.
    
    # H0 : beta_i = 0 vs H1 : beta_i != 0
    
    pval<- get_p_val(xi,St,y)
    
    if( pval < alpha){
      
      St <- cbind(St,xi)
      
      wealth <- wealth - alpha + 0.5
      
      St_index <- c(St_index, i + num_feature_used)
      
    }else{
      
      wealth <- wealth - alpha
      
    }
    
  }
  time = proc.time()[3] -t_start
  
  return(list(time =time,St_index = St_index, St = St , wealth = wealth, alpha = alpha, num_feature_used = num_feature_used + dim(Gt)[2]))
}

Streaming_SMLE_0 <- function(Memory, Gt, y,k){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  s <- dim(Gt)[2]
  
  num_feature_used <- Memory$num_feature_used
  
  residual <- Memory$residual
  
  beta_t <- Memory$beta_t # k dimensional list (numeric)
  
  ###########################################

  Xt<- cbind(St,Gt)
  
  Xt_index  <- c(St_index,(1:s)+num_feature_used)
  
  t_start = proc.time()[3]
  
  ###########################################
  # All features retained
  
  if(is.null(beta_t)){
    
    fit <- glm(y~.,data = data.frame(x = Xt, y = y))
    
    beta_t <- fit$coefficients[-1]
    
    residual <- fit$residuals
    
    St <- Xt
    
    St_index <- Xt_index
    
    time = proc.time()[3] - t_start
    
    return(list(time = time, St_index = St_index, St = St, beta_t =  beta_t, residual = residual,num_feature_used = num_feature_used + s))
    
  }
  
  ###########################################
  
  # U -serach 

  UR_inf = norm(crossprod(Xt[,1:k],residual),'i')
  U_inf_R_inf = norm(Xt[,1:k],'i')*max(abs(residual))
  u = UR_inf/U_inf_R_inf/500

  ###########################################
  
  beta_t_temp <- matrix(c(beta_t,rep(0,s)),ncol = 1) # k+s vector (n by 1 matrix)
  
  beta_t_tilde <- beta_t_temp + u*(matrix(c(crossprod(Xt[,1:k],residual),crossprod(Xt[,(k+1):dim(Xt)[2]],y)),ncol =1))
  
  beta_t <- Hard( beta_t_tilde , k) # k vector
  ###########################################
  
  retained_index <-which(beta_t_tilde %in% beta_t)

  retained_set <- Xt[,retained_index]

  ###########################################
  
  residual <- y - retained_set%*% matrix( beta_t[retained_index], ncol =1 )
  
  St_index <- Xt_index[retained_index]
  
  time = proc.time()[3] -t_start
   
  return(list(time = time, St_index = St_index, St = retained_set, beta_t =  beta_t[retained_index], residual = residual,num_feature_used = num_feature_used + s))
  
}

Streaming_SMLE_1 <- function(Memory, Gt, y,k){
  
  St <- Memory$St
  
  St_index <- Memory$St_index

  s <- dim(Gt)[2]
  
  num_feature_used <- Memory$num_feature_used
  
  residual <- Memory$residual
  
  beta_t <- Memory$beta_t # k dimensional list (numeric)
  
  ###########################################
  
  Xt<- cbind(St,Gt)
  
  Xt_index  <- c(St_index,(1:s)+num_feature_used)
  
  t_start = proc.time()[3]
  
  ###########################################
  # All features retained
  
  if(is.null(beta_t)){
    
    fit <- glm(y~.,data = data.frame(x = Xt, y = y))
    
    beta_t <- fit$coefficients[-1]
    
    residual <- fit$residuals
    
    St <- Xt
    
    St_index <- Xt_index
    
    time = proc.time()[3] - t_start
    
    return(list(time = time, St_index = St_index, St = St, beta_t =  beta_t, residual = residual,num_feature_used = num_feature_used + s))
    
  }
  
  ###########################################
  
  # U -serach 
  
  ###########################################
  
  beta_t_temp <- matrix(c(beta_t,rep(0,s)),ncol = 1) # k+s vector (n by 1 matrix)
  
  beta_update <- matrix(c(crossprod(Xt[,1:k],residual),crossprod(Xt[,(k+1):dim(Xt)[2]],y)),ncol =1)
  
  u  <- min(abs(beta_t))/max(abs(beta_update))
  
  beta_t_tilde <- beta_t_temp + u*beta_update
  
  beta_t <- Hard( beta_t_tilde , k) # k vector
  
  ###########################################
  
  retained_index <-which(beta_t_tilde %in% beta_t)
  
  retained_set <- Xt[,retained_index]
  
  ###########################################
  
  residual <- y - retained_set%*% matrix( beta_t[retained_index], ncol =1 )
  
  St_index <- Xt_index[retained_index]
  
  time = proc.time()[3] -t_start
  
  return(list(time = time, St_index = St_index, St = retained_set, beta_t =  beta_t[retained_index], residual = residual,num_feature_used = num_feature_used + s))
  
}

Streaming_SMLE_fortran <- function(Memory, Gt, y, k,family, l = 20){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  s <- dim(Gt)[2]
  
  num_feature_used <- Memory$num_feature_used
  
  residual <- Memory$residual
  
  beta_t <- Memory$beta_t # k dimensional list (numeric)
  
  ###########################################
  
  Xt<- cbind(St,Gt)
  
  Xt_index  <- c(St_index,(1:s)+num_feature_used)
  
  t_start = proc.time()[3]
  
  ###########################################
  # All features retained
  
  if(is.null(beta_t)){
    
    fit <- glm(y~.,data = data.frame(x = Xt, y = y))
    
    beta_t <- fit$coefficients[-1]
    
    residual <- fit$residuals
    
    St <- Xt
    
    St_index <- Xt_index
    
    time = proc.time()[3] - t_start
    
    return(list(time = time, St_index = St_index, St = St, beta_t =  beta_t, 
                
                residual = residual, num_feature_used = num_feature_used + s,
                
                ll = logLik(glm(y~.,data.frame(X = St,y = y),family = gaussian()))))
    
  }
  
  beta_t_temp <- matrix(c(beta_t,rep(0,s)),ncol = 1) # k+s vector (n by 1 matrix)
  
  goodness_of_fit <- crossprod(St,residual)
  
  new_margin_effect <- crossprod(Gt,y)
  
  beta_update <- matrix(c(goodness_of_fit,new_margin_effect),ncol =1)
  
  #ll_ <- unlist(lapply(1:l,FUN = u_select, beta_t, beta_t_temp,new_margin_effect, beta_update, Xt, y))
  
  ll_ <- .Fortran
  
  ll<- c( Memory$ll , ll_)
  
  best_model <- which( ll %in% max(ll) )
  
  best_model <- best_model[order(best_model)]
  
  if( 1 %in% best_model ){
    
    retained_index <- 1:k
    
    retained_set <- St
    
    residual <- Memory$residual
    
    St_index <- St_index
    
  }else{
    
    ii <- best_model[1]-1
    
    u <- sort(abs(beta_t),decreasing = FALSE)[ii]/sort(abs(new_margin_effect),decreasing = TRUE)[ii]
    
    if( u == 0){    
      
      retained_index <- 1:k
      
      retained_set <- St
      
      residual <- Memory$residual
      
      St_index <- St_index
      
    }else{
      
      beta_t_tilde <- beta_t_temp + u*beta_update
      
      beta_t <- Hard( beta_t_tilde , k) # k vector
      
      retained_index <- which(beta_t_tilde %in% beta_t)
      
      retained_set <- Xt[,retained_index]
      
      residual <- y - retained_set %*% matrix( beta_t[retained_index], ncol =1 )
      
      St_index <- Xt_index[retained_index]
      
    }
    
  }
  time = proc.time()[3] -t_start
  
  return(list(time = time, St_index = St_index, St = retained_set, 
              beta_t =  beta_t[retained_index], residual = residual,
              num_feature_used = num_feature_used + s,
              ll = logLik(glm(y~.,data.frame(x = retained_set,y = y),family = family))))
  
}

Streaming_SIS <- function(Memory, Gt, y,k){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  s <- dim(Gt)[2]
  
  num_feature_used <- Memory$num_feature_used
  
  correlation_St <- Memory$correlation_St

  ###########################################
  
  Xt<- cbind(St,Gt)
  
  Xt_index  <- c(St_index,(1:s)+num_feature_used)
  
  t_start = proc.time()[3]
  
  ###########################################


  if(is.null(correlation_St)){
    
    correlation_St <- crossprod(Xt, y)
    
    time = proc.time()[3] -t_start
    
    return(list(time = time, St_index = Xt_index, St = Xt, correlation_St =  correlation_St, num_feature_used = num_feature_used + s))
    
  }
  
  correlation_Gt <- crossprod(Gt, y)
  
  retained_index <- order(abs(c(correlation_St,correlation_Gt)),decreasing = TRUE)[1:k]
  
  retained_set <- Xt[,retained_index]
  
  St_index <- Xt_index[retained_index]
  
  correlation_St <- crossprod(retained_set, y)
  
  time = proc.time()[3] -t_start
  
  if( length(retained_index)!= k) {stop()}
  
  return(list(time = time, St_index = St_index, St = retained_set, correlation_St =  correlation_St, num_feature_used = num_feature_used + s))
  
}

OFS_fisher <- function(Memory, Gt, y , alpha = 0.05 ){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  num_feature_used <- Memory$num_feature_used
  
  t_start = proc.time()[3]
  
  for( i in 1:dim(Gt)[2] ){
    
    n <- dim(Gt)[1]
    
    xi <- matrix(Gt[,i],ncol = 1 )
    
    r <- cor(xi, y)
  
    t <- (0.5*log( (1+r)/(1-r) ))*(sqrt(n - 3))
    
    # get p value
    
    # H0 : beta_i = 0 vs H1 : beta_i != 0
    
    cutoff  = qnorm(1-alpha/2)
    
    if( abs(t) >= cutoff ){
      
      # reject H0 : z = 0
      
      # xi related to y, the feature passed the related test and be added to the model
      
      St <- cbind(St, xi)
      
      St_index<- c(St_index, i+num_feature_used)
      
      pp <- length(St_index)
      
      # redundant test when more than three features are included.
      
      screened_index <- c()
      
      if(pp >= 3){
        
        for( j in 1:length(St_index)){
          
          temp_cond <- St[,-j]
          
          xj <- matrix(St[,j],ncol =1)

          drop_off <- cond_dep_fisher_z(xj, y, temp_cond, alpha)
          
          if(drop_off){ screened_index <- c(screened_index,j) }
          
        }
        
        if(! is.null(screened_index)){

          St_index <- St_index[-screened_index]
          
          St <- St[,-screened_index]
          
        }

      }
      
    }
    
  }
  
  time = proc.time()[3] -t_start

  return(list(time =time, St_index = St_index, St = St , num_feature_used = num_feature_used + dim(Gt)[2]))
  
}

Saola_mi <- function(Memory, Gt, y , threshold = 0 ){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  num_feature_used <- Memory$num_feature_used
  
  t_start = proc.time()[3]
  
  for( i in 1:dim(Gt)[2] ){
    
    n <- dim(Gt)[1]
    
    xi <- matrix(Gt[,i], ncol = 1 )
    
    dep_i <- symmetrical_uncertainty(xi,y)
    
    screened_index <- c()
    
    if( dep_i > threshold ){
      
      St <- cbind(St, xi)
      
      St_index<- c(St_index, i+num_feature_used)
      
      pp <- length(St_index)

      # redundant test when more than three features are included.
      
      if(pp >= 3){
        
        for( j in 1:(length(St_index)-1)){
          
          xj <- matrix( St[,j], ncol =1 )

          dep_ij <- symmetrical_uncertainty(xi,xj)
          
          dep_j <- symmetrical_uncertainty(xj,y)
          
          if( dep_j > dep_i && dep_ij >= dep_i){
            
            # the feature is redundant
            
            St_index <- St_index[- length(St_index)]
            
            St <- St[,-length(St_index)]
            
            break
            
          }
          
          if(dep_i >  dep_j && dep_ij >= dep_j ){
            
            # there are redundant features in active set
            
            screened_index <- c(screened_index,j)
            
          }
          
        }        
        
        if(! is.null(screened_index)){

          St_index <- St_index[-screened_index]
          
          St <- St[,-screened_index]
          
        }


      }
      
    }
    
  }
  
  time = proc.time()[3] -t_start

  return(list(time =time, St_index = St_index, St = St , num_feature_used = num_feature_used + dim(Gt)[2]))
  
}

Saola_z <- function(Memory, Gt, y ){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  num_feature_used <- Memory$num_feature_used
  
  t_start = proc.time()[3]
  
  for( i in 1:dim(Gt)[2] ){
    
    n <- dim(Gt)[1]
    
    xi <- matrix(Gt[,i], ncol = 1 )
    
    ri <- cor(xi,y)
    
    St <- cbind(St, xi)
      
    St_index<- c(St_index, i+num_feature_used)
      
    pp <- length(St_index)
      
    # redundant test when more than three features are included.
      
    if(pp >= 3){
      
      St1 <- St
      
      St1_index <- St_index
        
      for( j in 1:(length(St_index)-1)){
          
        xj <- matrix( St[,j], ncol = 1 )
          
        rj <- cor(xj,y)
          
        rij <- cor(xi,xj)
          
        if( abs(rj) >= abs(ri) && abs(rij) > abs(ri)){
            
          # the feature is redundant
            
          St1_index <- St_index[ - length(St_index) ]
            
          St1 <- St[,-length(St_index)]
            
          break
            
          }
          
          if(abs(ri) >  abs(rj) && abs(rij) > abs(rj) ){
            
            # there are redundant features in active set
            
            id <- apply(St1,2,function(x){identical(x,St[,j])})
            
            if(sum(id)){
              
              St1 <- St1[,which(!id)]
              
              St1_index <- St1_index[which(!id)]
              
            }
            
            
          }
          
        }        
        
          
          St <- St1
          
          St_index <- St1_index
          
        }
        
        
    }
      
    
    
  
  time = proc.time()[3] -t_start
  
  return(list(time =time, St_index = St_index, St = St , num_feature_used = num_feature_used + dim(Gt)[2]))
  
}

BANS <- function(Memory, Gt, y, k,family, l = 5){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  s <- dim(Gt)[2]
  
  num_feature_used <- Memory$num_feature_used
  
  residual <- Memory$residual
  
  beta_t <- Memory$beta_t # k dimensional list (numeric)
  
  ###########################################
  
  Xt <- cbind(St,Gt)
  
  if( s != 1){Xt_index  <- c(St_index,(1:s)+num_feature_used)}else{Xt_index  <- c(St_index,1+num_feature_used)}
  
  ###########################################
  # All features retained
  
  if(dim(Xt)[2] < k){
    
    # Limits not reached, keep all features 
    
    fit <- fastglm(x = Xt, y = y,family = family)
    
    beta_t <- fit$coefficients
    
    residual <- fit$residuals
    
    St <- Xt
    
    St_index <- Xt_index
    
    return(list(time = 0, St_index = St_index, St = St, beta_t =  beta_t, 
                
                residual = residual, num_feature_used = num_feature_used + s,
                
                ll = -Inf))
    
  }else if(is.null(beta_t)){
    
    # Limits reached when first group features came
    
    fit <- SMLE(X = Xt , Y = y, family = family , k = k)
    
    beta_t <- fit$coef_retained
    
    residual <- y - predict(fit)
    
    St <- Xt[,fit$ID_retained]
    
    St_index <- Xt_index[fit$ID_retained]
    
    return(list(time = 0, St_index = St_index, St = St, beta_t =  beta_t, 
                
                residual = residual, num_feature_used = num_feature_used + s,
                
                ll = logLik(fit)))
  }
  
  t_start = proc.time()[3]
  
  beta_t_temp <- matrix(c(beta_t,rep(0,s)),ncol = 1) # k+s vector (n by 1 matrix)
  
  goodness_of_fit <- crossprod(St,residual)
  
  new_margin_effect <- crossprod(Gt,(y-family$linkinv(rep(0,length(y))))) 
  
  beta_update <- matrix(c( goodness_of_fit, new_margin_effect ),ncol = 1)
  
  tryCatch(dim(beta_update) != dim(beta_t_temp))
  
  ll_ <- unlist(lapply(1:l,FUN = u_select, beta_t, beta_t_temp,new_margin_effect, beta_update, Xt, y,family))
  
  ll <- c( Memory$ll , ll_)
  
  best_model <- sort(which( ll %in% max(ll) ),decreasing = FALSE)
  
  if( 1 %in% best_model ){
    
    # Keep the original model 
    
    beta_t <- beta_t[1:k]
    
    retained_set <- Xt[,1:k]
    
    residual <- Memory$residual
    
    St_index <- St_index
    
  }else{
    
    ii <- best_model[1]-1

    if(ii > s ){ u <- sort(abs(beta_t),decreasing = FALSE)[ii]/max(abs(new_margin_effect))}else{
      
      u <- sort(abs(beta_t),decreasing = FALSE)[ii]/sort(abs(new_margin_effect),decreasing = TRUE)[ii]
      
    }
    
    if( u == 0){    
      
      beta_t <- beta_t[1:k]
      
      retained_set <- Xt[,1:k]
      
      residual <- Memory$residual
      
      St_index <- St_index
      
    }else{
      
      beta_t_tilde <- beta_t_temp + u * beta_update
      
      beta_t1 <- Hard( beta_t_tilde , k) # k vector
      
      beta_t <- beta_t1[ beta_t1 != 0 ]
      
      retained_index <- which(beta_t_tilde %in% beta_t)
      
      retained_set <- Xt[,retained_index]
      
      residual <- y - family$linkinv(retained_set %*% matrix( beta_t, ncol =1 ))
      
      St_index <- Xt_index[retained_index]
      
    }
    
  }
  time = proc.time()[3] -t_start
  
  return(list(time = time, St_index = St_index, St = retained_set, 
              beta_t =  beta_t, residual = residual,
              num_feature_used = num_feature_used + s,
              ll = logLik(fastglm(x = retained_set,y = y,family =family))
  ))
  
}

OSFS_FI <- function(Memory, Gt, y,gamma =0.3){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  num_feature_used <- Memory$num_feature_used
  
  t_start = proc.time()[3]
  
  for( i in 1:dim(Gt)[2] ){
    
    n <- dim(Gt)[1]
    
    xi <- matrix(Gt[,i], ncol = 1 )
    
    ri <- cor(xi,y)
    
    int_val <- dep_interaction(xi, y , St, alpha =0.05 , ri)
    
    if(int_val > gamma){
      
      St <- cbind(St, xi)
      
      St_index<- c(St_index, i+num_feature_used)
      
    }else if(int_val > 0){
      
      St <- cbind(St, xi)
    
      St_index<- c(St_index, i+num_feature_used)
    
      pp <- length(St_index)
    
    # redundant test when more than three features are included.
    
      if(pp >= 3){
        
        St1 <- St
        
        St1_index <- St_index
        
        for( j in 1:(length(St_index)-1)){
          
          xj <- matrix( St[,j], ncol = 1 )
          
          rj <- cor(xj,y)
          
          rij <- cor(xi,xj)
          
          if( abs(rj) >= abs(ri) && abs(rij) > abs(ri)){
            
            # the feature is redundant
            
            St1_index <- St_index[ - length(St_index) ]
            
            St1 <- St[,-length(St_index)]
            
            break
            
          }
          
          if(abs(ri) >  abs(rj) && abs(rij) > abs(rj) ){
            
            # there are redundant features in active set
            
            id <- apply(St1,2,function(x){identical(x,St[,j])})
            
            if(sum(id)){
              
              St1 <- St1[,which(!id)]
              
              St1_index <- St1_index[which(!id)]
              
            }
            
            
          }
          
        }        
        
        
        St <- St1
        
        if(is.null(dim(St))){St <- matrix(St,ncol=1)}
        
        St_index <- St1_index
        
      }
    
    }
  }
  
  time = proc.time()[3] -t_start
  
  return(list(time =time, St_index = St_index, St = St , num_feature_used = num_feature_used + dim(Gt)[2]))
  
}

OSFS_A3M <- function(Memory, Gt, y){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  num_feature_used <- Memory$num_feature_used
  
  if(is.null(Memory$dep_Mean)){dep_Mean <- 0}else{dep_Mean <- Memory$dep_Mean } 
     
  if(is.null(Memory$dep_Set)){dep_Set <- 0}else{dep_Set <- Memory$dep_Set }
  
  if(is.null(Memory$dep_Array)){dep_Array <- list()}else{dep_Array <- Memory$dep_Array } 

  Y <- as.matrix(y,ncol=1)
  
  t_start = proc.time()[3]
  
  for( i in 1:dim(Gt)[2] ){
    
    xi <- matrix(Gt[,i], ncol = 1 )
    
    dep_single <- dep_an(xi,Y)
    
    dep_Array <- unlist(c(dep_Array,dep_single))
    
    if(dep_single > dep_Set){
      
      St <- cbind(St, xi)
      
      St_index<- c(St_index, i+num_feature_used)

      dep_New <- dep_an(St,Y)
      
      if(dep_New > dep_Set){
        
        dep_Set <- dep_New
        
        dep_Mean <- sum(dep_Array[St_index])/length(St_index)
        
      } else if(dep_New == dep_Set){
        
        index_del <- non_signf(St,Y)
        
        St <- St[,-index_del]
        
        St_index <- St_index[-index_del]
        
        dep_Mean <- sum(dep_Array[St_index])/length(St_index)
        
      } else {
        
        St <- St[,-length(St_index)]
        
        St_index <- St_index[-length(St_index)]
      }
      
    }

  }
  
  time = proc.time()[3] -t_start
  
  return(list(time =time, St_index = St_index, St = St , num_feature_used = num_feature_used + dim(Gt)[2],
         dep_Set =dep_Set, dep_Mean = dep_Mean))
  
}

OSFS_Density <- function(Memory, Gt, y){
  
  St <- Memory$St
  
  St_index <- Memory$St_index
  
  num_feature_used <- Memory$num_feature_used
  
  if(is.null(Memory$dep_Mean)){dep_Mean <- 0}else{dep_Mean <- Memory$dep_Mean } 
  
  if(is.null(Memory$dep_Set)){dep_Set <- 0}else{dep_Set <- Memory$dep_Set }
  
  if(is.null(Memory$dep_Array)){dep_Array <- list()}else{dep_Array <- Memory$dep_Array } 
  
  Y <- as.matrix(y,ncol=1)
  
  t_start = proc.time()[3]
  
  for( i in 1:dim(Gt)[2] ){
    
    xi <- matrix(Gt[,i], ncol = 1 )
    
    dep_single <- dep_an2(xi,Y)
    
    dep_Array <- unlist(c(dep_Array,dep_single))
    
    if(dep_single > dep_Mean){
      
      St <- cbind(St, xi)
      
      St_index<- c(St_index, i+num_feature_used)
      
      dep_New <- dep_an2(St,Y)
      
      if(dep_New > dep_Set){
        
        dep_Set <- dep_New
        
        dep_Mean <- sum(dep_Array[St_index])/length(St_index)
        
        next
      } 
      
      if( abs(dep_New - dep_Set) / dep_Set <= 0.05){
        
        St_temp <- St
        
        St_index_temp <- St_index
        
        index_del <- non_signf2(St,Y)
        
        St <- St_temp[,-index_del]
        
        St_index <- St_index[-index_del]
        
        dep_Del <- dep_an2(St,Y)
        
        if (dep_Del >= dep_Set) {
          
          dep_Set <- dep_Del
          
          dep_Mean <- sum(dep_Array[St_index])/length(St_index)
          
        } else {
          
         St <- St_temp
         
         St_index <- St_index_temp

        }
        
        next
        
      }
      
      St <- St[,-length(St_index)]
        
      St_index <- St_index[-length(St_index)]
      
    }
    
  }
  
  time = proc.time()[3] -t_start
  
  return(list(time =time, St_index = St_index, St = St , num_feature_used = num_feature_used + dim(Gt)[2],
              dep_Set =dep_Set, dep_Mean = dep_Mean,dep_Array =dep_Array))
  
}
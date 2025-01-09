library( infotheo )
library(glmnet)
logit_loss <- function(actual_labels, predicted_probabilities) {
  -mean(actual_labels * log(predicted_probabilities) + 
          (1 - actual_labels) * log(1 - predicted_probabilities))
}


IHT_w_ucheck <-function(Beta_s, U_0 , v_s, k , Y, X, family, U_rate, tol, max_iter,p){
  
  i <- 0
  
  LH <- rep(0,max_iter)
  
  FD <- NULL
  
  repeat{
    
    i <- i+1
    
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
      
      
      if( MSE < tol || i >= max_iter  || (LH[i]-LH[i-1])< valid_LH_diff){break}
      
    }
    
    Beta_s <- Beta_t
    
    ID_None_s <- sort((1:p)[Beta_s!= 0])
    
    coef_None_s <- as.matrix(Beta_s[ID_None_s],ncol=1)
    
    Xs_s <- X[, ID_None_s]
    
    theta_s  <- make.link(family$link)$linkinv(Xs_s %*% coef_None_s)
    
    v_s <- crossprod(X, Y - theta_s)
    
    u_s  <- U_0
    
    LH_last <- LH[i]
    
    i <- i + 1
    
  }
  
  IHT_res <- list(rho= v_s, ID_None_s = ID_None_s, Beta_s = Beta_s,LH = LH_last)
  
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
  eta <- X[,ind]%*%as.matrix(beta[ind],nrow = length(ind))
  mu <- linkinv(eta)
  dev <- sum(dev.resids(y, mu, 1))
  return((-1/2)*aic(y, n, mu, 1, dev))
}

find_index<-function(A,B){
  
  # find index of A as a submatrix of B
  
  index<-apply(A,2,FUN=function(y){(1:dim(B)[2])[apply(B,2,FUN = function(x){identical(x,y)})]})
  
  index
  
}

alpha_investing <- function(x, y){
  
  #  configure parameters 
  
  ## 0.5 recommended by the paper
  wealth = 0.5
  
  delta_alpha = 0.5
  
  n <- dim(x)[1]
  
  p <- dim(x)[2]
  
  model <- {}
  
  model_index<-{}
  
  for( i in 1:p ){
    
    xi<- x[,i]
    
    alpha<- wealth/(2*i)
    
    # get p value from git package.
    
    # H0 : beta_i = 0 vs H1 : beta_i != 0
    
    pval<- get_p_val(xi,model,y)
    
    if( pval < alpha){
      
      model <- cbind(model,xi)
      
      model_index <- c(model_index,i)
      
      wealth <- wealth - alpha+ delta_alpha
      
      
    }else{
      
      wealth <- wealth - alpha
      
    }
  }
  model_index
}

Group_alpha_investing_test <- function(x,y){
  
  wealth = 0.5
  
  delta_alpha = 0.5
  
  n <- dim(x)[1]
  
  p <- dim(x)[2]
  
  model <- {}
  
  model_index<-{}
  
  s = 20
  
  for( j in 1 : (p/s)){
    
    Gt <- x[,((j-1)*s+1):(j*s)]
    
    num_feature_used <- (j-1)*s
    
    for( i in 1:dim(Gt)[2] ){
      
      xi<- Gt[,i]
      
      alpha<- wealth/(2*(i+num_feature_used))
      
      # get p value from git package.
      
      # H0 : beta_i = 0 vs H1 : beta_i != 0
      
      pval<- get_p_val(xi,model,y)
      
      if( pval < alpha){
        
        model <- cbind(model,xi)
        
        wealth <- wealth - alpha + 0.5
        
        model_index <- c(model_index,i+num_feature_used)
        
      }else{
        
        wealth <- wealth - alpha
        
      }
      
    }
    
  
  }

  model_index
}




get_p_val<-function(x, model, y){
  
  model<- cbind(model,x)
  
  fit <- glm.fit(model,y)
  
  class(fit) <- "glm"
  
  if(sum(is.na(fit$coefficients))){pval = 1}else{pval <- coefficients(summary(fit))[dim(model)[2],4]}
  
  pval
}

Hard<-function(t,k,lam =1)
{
  #----------------------------------------------------------------------------#
  # Hard threshold function                                                    #
  #----------------------------------------------------------------------------#
  # This function keep the k elements with the larges absolute value and       #
  # truncated the rest to zero.                                                #
  #----------------------------------------------------------------------------#
  # Parameters:                                                                #
  # t : A list of values.                                                      #
  # k:  number of values retained.                                             #
  #----------------------------------------------------------------------------#
  
  
  y<-as.vector(t)
  t<-abs(y)
  
  lam<-sort(t, decreasing=TRUE)[k]
  
  y[t<lam]<-0

  return(y)
}

# install.packages("infotheo")

saola<-function(x,y,threshold){
  
  n <- dim(x)[1]
  
  p <- dim(x)[2]
  
  model <- {}
  
  model_index<-{}
  
  for( i in 1:p ){
    
    xi <- x[,i]
    
    dep_score <- symmetrical_uncertainty(xi,y)
    
    if(dep_score > threshold){
      
      # Feature related to response  
      
      model_index <- c(model_index , i)
      
      model <- cbind(model, xi)
      
      pp <- length(model_index)
      
      if( pp > 3 ){
        
        remove_index <- NULL
        
        for ( j in 1: pp){
          
          dep_score_ij<- symmetrical_uncertainty(xi,model[,j])
          
          if( dep_score_ij > threshold){
            
            # check related features in candidate set.
            
            if(dep_score_ij > dep_score){
              
              #remove the last feature in the candidate set.
              
              model<- model[,-dim(model)[2]]
              
              model_index <- model_index[-length(model_index)]
              
              break
              
            }else{
              
              # remove j th feature in the candidate set.
              
              remove_index <-c( remove_index, j) 
            }
            
          }
          
        }
        if(!is.null(remove_index)){
          
          model_index<- model_index[-remove_index]
          
          model <- model[,-remove_index]
          
        }
        
      } 
      
    }
    
  }
  
  
  model_index 
  
}

symmetrical_uncertainty<- function(x,y){
  
  x <- discretize(x,"equalwidth")
  
  hx <- entropy(x)
  
  y <- discretize(y,"equalwidth")
  
  hy <- entropy(y)
  
  mutual_info <- mutinformation(x,y)
  
  score <- (2 * mutual_info)/(hx + hy)  
  
  score
}



cond_dep_fisher_z <- function(x,y,condition,alpha = 0.05){
  
  cond <- c()
  
  drop_off <- 0
  
  for( i in 1:dim(condition)[2]){
    
    cond <- cbind(cond,condition[,i])
    
    data <- cbind(x,y,cond)
    
    n <- dim(data)[1]
    
    p <- dim(cond)[2] + 2
    
    Cor_matrix <- cor(data)
    
    S <- Cor_matrix[1:2,1:2] - Cor_matrix[1:2,3:p] %*% solve(Cor_matrix[3:p,3:p]) %*% Cor_matrix[3:p,1:2]
    
    r <- S[1,2] /sqrt(S[1,1]*S[2,2])
    
    z <- (0.5*log( (1+r)/(1-r) ))*(sqrt(n - p- 3))
    
    cutoff  = qnorm(1-alpha/2)
    
    if( abs(z) < cutoff ){
      
      drop_off = 1
      
      return(drop_off)
      
    }
      
  }
  
  return(drop_off)

}
cond_dep_fisher_z_r <- function(x,y,condition,alpha = 0.05){
  
  cond <- c()
  
  drop_off <- 0
  
  for( i in 1:dim(condition)[2]){
    
    cond <- cbind(cond,condition[,i])
    
    data <- cbind(x,y,cond)
    
    n <- dim(data)[1]
    
    p <- dim(cond)[2] + 2
    
    Cor_matrix <- cor(data)
    
    S <- Cor_matrix[1:2,1:2] - Cor_matrix[1:2,3:p] %*% solve(Cor_matrix[3:p,3:p]) %*% Cor_matrix[3:p,1:2]
    
    r <- S[1,2] /sqrt(S[1,1]*S[2,2])
    
  }
  
  return(abs(r))
  
}

dep_interaction <- function(xi, y , St, alpha,ri){
  
  p <- dim(St)[2]
  
  if(is.null(St)){return(1)}
  
  else{  
    
    if(dim(St)[2]==0){return(1)}
    
    val <- apply(St,2,FUN =cond_dep_fisher_z_r,y=y,condition = xi,alpha = alpha )
  
  return(mean(val)-ri)
    
    }
  

}



u_select <- function(i, beta_t, beta_t_temp,new_margin_effect, beta_update, Xt, y, family = gaussian()){
  
  if( i > length(new_margin_effect)){
    
    u <- sort(abs(beta_t),decreasing = FALSE)[i]/min(abs(new_margin_effect))
    
  }else{
    
    u <- sort(abs(beta_t),decreasing = FALSE)[i]/sort(abs(new_margin_effect),decreasing = TRUE)[i]
    
    }
  
  
  
  beta_t_tilde_i <- beta_t_temp + u * beta_update
  
  beta_t1 <- Hard( beta_t_tilde_i , k) # k vector
  
  beta_t_i <- beta_t1[ beta_t1 != 0 ] 
  
  retained_index <- which( beta_t_tilde_i %in% beta_t_i )

  fit <- fastglm(x = Xt[,retained_index], y = y, family = family)
  
  return(  logLik( fit ))
  
}

standardize <- function(x){
  
  y <- apply(x,2,FUN = function(xi){ (xi-mean(xi))/sqrt(var(xi)) })
  
  matrix(y,nrow = dim(x)[1])
}


card <- function(sets,Y,label,N){
  D <- sort(sets,index.return=TRUE)
  
  min_d <- D$x[2]
  max_d <- D$x[N]
  mean_d <- 1.5*(max_d-min_d)/(N-2)
  
  cNum <- 0
  cTotal <- 1
  ind2 <- D$ix[2]
  if(Y[ind2,1] == label){
    cNum <- cNum + 1
  }
  
  for(j in 3:N){
    if(D$x[j] - D$x[j-1] < mean_d){
      ind <- D$ix[j]
      cTotal <- cTotal + 1
      if(Y[ind,1] == label){
        cNum <- cNum + 1
      }
    } else {
      break
    }
  }
  
  c <- cNum/cTotal
  return(c)
}

dep_an <- function(data, Y){
  n <- dim(data)[1]
  card_ND <- 0
  DArray <- as.matrix(stats::dist(data, method = "euclidean"))
  dep <- mean(unlist(lapply(1:n,FUN = function(i){card(DArray[,i], Y, Y[i,1], n)})))
  return(dep)
}

sig <- function(St,Y,j){
  
  d_B <- dep_an(St, Y)
  St<-St[,-j]
  d_F <- dep_an(St, Y)

  return((d_B - d_F) / d_B)
}

non_signf <- function(St,Y){
  
  Num <- dim(St)[2]-1
  screening_index <- rep(0,Num)
  A <- sample(Num)
  
  for(i in 1:Num){
    rnd <- A[i]
    if(sig(St,Y,rnd) == 0){
      screening_index[rnd] <- 1
    }
  }
  index_del <- (1:Num)[screening_index == 1]
  return(index_del)
}


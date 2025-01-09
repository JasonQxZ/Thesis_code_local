library( infotheo )
library(glmnet)

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
    
    u <- sort(abs(beta_t),decreasing = FALSE)[i]/max(abs(new_margin_effect))
    
  }else{
    
    u <- sort(abs(beta_t),decreasing = FALSE)[i]/sort(abs(new_margin_effect),decreasing = TRUE)[i]
    
    }
  
  beta_t_tilde_i <- beta_t_temp + u * beta_update
  
  beta_t1 <- Hard( beta_t_tilde_i , k) # k vector
  
  beta_t_i <- beta_t1[ beta_t1 != 0 ] 
  
  retained_index <- which( beta_t_tilde_i %in% beta_t_i )

  fit <- fastglm(x= Xt[,retained_index],y = y,family =family)
  
  return(  logLik( fit ))
  
}

standardize <- function(x){
  
  y <- apply(x,2,FUN = function(xi){ (xi-mean(xi))/sqrt(var(xi)) })
  
  matrix(y,nrow = dim(x)[1])
}

card_density <- function(sets, Y, label, N) {
  sets<- as.matrix(sets,nrow = Y)
  D <- sort(sets, index.return = TRUE)
  cTotal <- 0
  cNum <- 0
  last_density <- 0
  if(is.na(D$x[2])){
    return(0)
  }
  
  if (D$x[2] >= 0) {
    last_density <- D$x[2]
  }
  for (j in 2:N) {
    density <- D$x[j] / (j - 1)
    ind <- D$ix[j]
    if (density == 0 || density <= last_density) {
      last_density <- density
      cTotal <- cTotal + 1
      if (Y[ind, 1] == label) {
        cNum <- cNum + 1
      }
    } else {
      break
    }
  }
  c <- cNum / cTotal
  return(c)
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
  Y <- as.matrix(Y, ncol =1)
  n <- nrow(Y)
  data <- as.matrix(data,nrow = n)
  card_U <- n
  card_ND <- 0
  D <- stats::dist(data, method = "euclidean")
  DArray <- as.matrix(D)
  for(i in 1:n){
    d <- DArray[,i]
    class <- Y[i,1]
    card_ND <- card_ND + card(d, Y, class, n)
  }
  dep <- card_ND / card_U
  return(dep)
}

dep_an2 <- function(data, Y){
  
  Y <- as.matrix(Y, ncol =1)
  n <- nrow(Y)
  data <- as.matrix(data,nrow = n)
  card_U <- n
  card_ND <- 0
  D <- stats::dist(data, method = "euclidean")
  DArray <- as.matrix(D)
  for(i in 1:n){
    d <- as.matrix(DArray[,i],ncol=1)
    class <- Y[i,1]
    card_ND <- card_ND + card_density(d, Y, class, n)
  }
  dep <- card_ND / card_U
  return(dep)
}


sig <- function(St,Y,j){
  
  d_B <- dep_an(St, Y)
  St<-St[,-j]
  d_F <- dep_an(St, Y)

  return((d_B - d_F) / d_B)
}

sig2 <- function(St,Y,j){
  
  d_B <- dep_an2(St, Y)
  St<-St[,-j]
  d_F <- dep_an2(St, Y)
  
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
non_signf2 <- function(St,Y){
  Num <- dim(St)[2]-1
  screening_index <- rep(0,Num)
  A <- sample(Num)
  for(i in 1:Num){
    rnd <- A[i]
    if(sig2(St,Y,rnd) <= 0){
      screening_index[rnd] <- 1
    }
  }
  index_del <- (1:Num)[screening_index == 1]
  return(index_del)
}


library(MASS)

Set_A_h <- function( h, s, X , Y,family){
  
  # Input : \beta_tilde , h , s_tilde 
  
  # Ouput : A set of top h features in magnitude of \beta_tilde in s_tilde
  
  beta <- coefficients(glm.fit(x = X[,s],y = Y, family = family))
  
  out <- sort(s[order(abs(beta),decreasing = T)[1:h]])
  
  out
  
}

Set_B_h <- function( h , s , X , Y ,family){
  
  ##--------------------------------------------------------------------------##
  
  # Input : \beta_tilde , h , s_tilde
  
  # Output : A set of top h features in magnitude of  \rho = x_j ^t (Y - b'( X_s\beta_s)) in s_tilde compliment
  
  ##--------------------------------------------------------------------------##
  
  s_compliment <- sub_off( 1:dim(X)[2],s )
  
  beta <- coefficients(glm.fit(x = X[,s],y = Y, family = family))
  
  theta <- make.link(family$link)$linkinv(X[,s]%*%beta)
  
  rho <- crossprod(X[,s_compliment],(Y-theta))
  
  out <- sort(s_compliment[order(abs(rho),decreasing = T)[1:h]])
  
  out
  
}


Set_C_h <-  function(h , s, X , Y,family){
  
  ##--------------------------------------------------------------------------##
  
  # Input : \beta_tilde , h , s_tilde
  
  # Output : A set of top h features in magnitude of  \rho = x_j ^t (Y - b'( X_s\beta_s)) in s_tilde compliment
  
  ##--------------------------------------------------------------------------##

  s_compliment <- sub_off( 1:dim(X)[2],s )
  
  X_s <- X[,s]
  
  beta <- coefficients(glm.fit(x = X_s,y = Y, family = family))
  
  theta <- make.link(family$link)$linkinv(X_s%*%beta)

  X_sTX_s <- cor(X_s)
  
  PI_s <- X_s %*% (chol2inv(chol(X_sTX_s))) %*% t(X_s)
    
  rho <- crossprod(X[,s_compliment],(diag(I(length(Y)))-PI_s)%*%(Y-theta))
  
  out <- sort(s_compliment[order(abs(rho),decreasing = T)[1:h]])
  
  out
  
}

#-----------------       TEST       -----------------#
# X <- matrix(rnorm(100000),ncol=1000,nrow = 100)
# 
# beta_star <- c(2,3,4,-2,6,rep(0,995))
# 
# Y <- X %*% beta_star + rnorm(100,0,0.5)
# 
# s <- c(1,2,3,6,7)
# 
# beta_s <- coefficients(glm.fit(x = X[,s], y = Y))
# 
# beta <- rep(0,1000)
# 
# beta[s] <- beta_s
# 
# 
# 
# print(Set_A_h( 3 , s, X, Y,family =gaussian()))
# 
# print(Set_B_h(beta, 10, s, X, Y, family = gaussian()))
# 
# print(Set_C_h(beta, 10, s, X, Y, family = gaussian()))
#----------------------------------------------------#


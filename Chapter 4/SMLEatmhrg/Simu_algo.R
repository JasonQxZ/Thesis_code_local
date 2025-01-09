Run_method <- function(METD , k , family , Data ){
  
  switch(METD, 
         
         SIS = run_SIS(x = Data$X,y = Data$Y, k = k , family = family),
         
         Lasso = run_Lasso(x = Data$X, y = Data$Y,  k = k, family = family),
         
         Sbess_S1 = run_Sbess( x = Data$X,y = Data$Y, k = k, family = family , Strategy = 1 ),
         
         Sbess_S2 = run_Sbess( x = Data$X,y = Data$Y, k = k, family = family , Strategy = 2 ),
         
         Sbess_S3 = run_Sbess( x = Data$X,y = Data$Y, k = k, family = family , Strategy = 3 ),
         
         Sbess = run_Sbess( x = Data$X,y = Data$Y, k = k, family = family , Strategy = 4 ),
         
         abess = run_abess( x = Data$X,y = Data$Y, k = k, family = family),
         
         SMLE_zero = run_SMLE_zero( x = Data$X,y = Data$Y, k = k, family = family),
  
         SMLE = run_SMLE_Lasso( x = Data$X,y = Data$Y, k = k, family = family),
        
         Holp = run_Holp( x = Data$X,y = Data$Y, k = k, family = family),
  
         Tilting = run_tilting(x = Data$X,y = Data$Y, k = k, family = family))
}


run_SIS <-  function(x , y,family , k ){
  
  SIS(x , y , family = family$family ,iter = F,nsis =k)$ix0
  
}

run_Lasso <-  function(x , y,family , k ){
  
  Lasso_fit <- glmnet(x = x, y = y, family = family$family, pmax = k)
  
  Lasso_index <- (1 : dim(x)[2])[ Lasso_fit$beta[, dim(Lasso_fit$beta)[2]] != 0 ]
  
  Lasso_index 
  
}

run_abess <-  function(x , y,family , k ){
  
  abess_fit <- abess(x = x, y = y, family = family$family , support.size = k)
  
  ind <- extract(abess_fit)$support.vars
  
  ind <- as.numeric(regmatches(ind, gregexpr("[[:digit:]]+", ind)))
  
  ind
  
}

run_SMLE_zero <- function(x , y, family , k ){
  
  ini <- rep(0,dim(x)[2])
  
  SMLE_fit <- SMLE(X =x,Y = y,family = family, categorical = FALSE, k = k, fast = T, coef_initial = ini,intercept  = F)

  
  SMLE_fit$ID_retained 
}

run_SMLE_Lasso <- function(x , y, family , k ){

  SMLE_fit <- SMLE(X =x,Y = y,family = family, categorical = FALSE, k = k, fast = T,intercept  = F)
  
  SMLE_fit$ID_retained 
}

run_Sbess <- function(x, y , k, family ,Strategy, ...){
  
  Sbess_fit <- Sbess(X = x,Y = y,family = family,k = k, Strategy = Strategy, ...)
  
  Sbess_fit$IHT_res_out$ID_None_s

}

run_Holp <- function(x , y ,k, family ){
  
  Holp_fit <- screening(x , y ,method = "holp",num.select = k, family=family$family)
  
  Holp_fit$screen
  
}

run_tilting<- function(x , y ,k, family ){
  
  tilting_fit <- tilting(x , y ,max.size = k, op = 2)
  
  tilting_fit$active
  
}



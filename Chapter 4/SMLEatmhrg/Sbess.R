library(glmnet)

Sbess<-function(Y, X, k, family=gaussian(), d = 2*k, m = 6,
                    max_iter = 500, tol = 10^(-3), U = 1,
                    U_rate= 0.5, lasso_initial = F, Strategy = c(1,2,3,4)
                ){
  
  
  # Initialize 
  
  LH <- 0
  
  FD <- NULL
  
  number_of_Ucheck <- rep(0,max_iter)
  
  n<-dim(X)[1]; p<-dim(X)[2]
  
  llh_out <- c()
  
  Beta_out <- c()
  
  ID_out <- c()
  
  if(lasso_initial){ 
    
    fit_pre <- glmnet(x = X, y = Y,family = family)
    
    Beta0 <- Hard( t = c( fit_pre$beta[,dim(fit_pre$beta)[2]] ), k = k)
    
    ID_None0 <- which(Beta0 != 0) 
    
    coef_None0 <- as.matrix( Beta0[ ID_None0 ] , ncol=1) 
    
    Xs_0 <- as.matrix( X[, ID_None0] )
    
    theta_0  <- Xs_0 %*% coef_None0
    
    U_0 <- U/(sqrt(p)*norm(Xs_0, "i")^2)
    
    theta_0 <- make.link(family$link)$linkinv(theta_0)
    
    v_s <- crossprod(X , Y - theta_0)
    
    
  }else{
    
    Beta0 <- rep( 0,dim(X)[2] )
    
  }
  
  Beta_s <- Beta0
  
  IHT_res <- IHT_w_ucheck(Beta_s, k , Y, X, family,U, U_rate, tol, max_iter, p)

  IHT_res_out <- IHT_res
  
  L <- IHT_res$LH 

  rho <- IHT_res$rho  # p*1
  
  ID_out <- cbind(ID_out, IHT_res$ID_None_s)
  
  Beta_out <- cbind(Beta_out, IHT_res$Beta_s)
  
  llh_out <- c(llh_out, IHT_res$LH)
  
  end_count <- 0
  
  repeat{
    
    # Splicing 
    if(Strategy == 1){
      
      Splicing_Set <- Set_B_h(h = d - k , s = IHT_res$ID_None_s, X = X , Y = Y , family = family)
      
      Candidate_Set <- sort(c(IHT_res$ID_None_s,Splicing_Set))
      
    }else if(Strategy == 2){
      
      Splicing_Set <- Set_C_h(h = d - k , s = IHT_res$ID_None_s, X = X , Y = Y , family = family)
      
      Candidate_Set <- sort(c(IHT_res$ID_None_s,Splicing_Set))
      
    }else if(Strategy == 3){
      
      A_m_S <- Set_A_h(h = m , s = IHT_res$ID_None_s, X = X , Y = Y , family = family)
      
      Splicing_Set <- Set_B_h(h = d- m , s = A_m_S, X = X , Y = Y , family = family)
      
      Candidate_Set <- sort(c(A_m_S,Splicing_Set))
      
    }else if(Strategy == 4){
      
      A_m_S <- Set_A_h(h = m , s = IHT_res$ID_None_s, X = X , Y = Y , family = family)
      
      Splicing_Set <- Set_C_h( h = d- m , s = A_m_S, X = X , Y = Y , family = family)
      
      Candidate_Set <- sort(c(A_m_S,Splicing_Set))
      
    }

    Beta_new <- rep(0,dim(X)[2])
    
    Beta_new[Candidate_Set] <- coefficients(glm.fit(x = X[,Candidate_Set],y = Y, family = family))
    
    IHT_res <- IHT_w_ucheck(Beta_new, k , Y, X, family, U, U_rate, tol, max_iter, p)
    
    ID_out <- cbind(ID_out, IHT_res$ID_None_s)
    
    Beta_out <- cbind(Beta_out, IHT_res$Beta_s)
    
    llh_out <- c(llh_out, IHT_res$LH)
    
    if(IHT_res$LH > L ){ 
      
      IHT_res_out <- IHT_res
      
      L <- IHT_res$LH 
      
      end_count <- 0
      
    }
    
    end_count <- end_count+1
    
    if(end_count == 3) { break }
      
  }
  
  fit<-list(llh_out = llh_out,ID_out = ID_out,Beta_out = Beta_out, IHT_res_out = IHT_res_out )
  
}
   
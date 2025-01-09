aSbess<-function(Y, X, family=gaussian(), m = 3,
                max_iter = 500, tol = 10^(-3), U = 1,
                U_rate= 0.5
){
  
  
  n<-dim(X)[1]; p<-dim(X)[2]
  
  k = n^(1/3)*log(p)
  
  # Initialize 

  LH <- 0
  
  FD <- NULL
  
  number_of_Ucheck <- rep(0,max_iter)

  ebic_out <- c()
  
  ID_out <- c()

  Beta_s <-  rep( 0,dim(X)[2] )
  
  IHT_res <- IHT_w_ucheck(Beta_s, k , Y, X, family,U, U_rate, tol, max_iter, p)
  
  E <- ebic(IHT_res$Beta_s, Y, X, family)
  
  IHT_res_out <- IHT_res
  
  rho <- IHT_res$rho  # p*1
  
  ID_out <- cbind(ID_out, IHT_res$ID_None_s)

  ebic_out <- c(ebic_out, E)
  
  end_count <- 0
  
  repeat{
    
    if( k+2 < p ){d <- k+2}else{d <- p}
    
    repeat{
      
      A_n <- Set_A_h(h = m - 1 , s = IHT_res$ID_None_s, X = X , Y = Y , family = family)
      
      A_m <- Set_A_h(h = m , s = IHT_res$ID_None_s, X = X , Y = Y , family = family)
    
      A_t <- Set_A_h(h = m + 1 , s = IHT_res$ID_None_s, X = X , Y = Y , family = family)
      
      ebic_n <- ebic(IHT_res$Beta_s, Y, X, family, s = A_n)
      
      ebic_m <- ebic(IHT_res$Beta_s, Y, X, family, s = A_m)
      
      ebic_t <- ebic(IHT_res$Beta_s, Y, X, family, s = A_t)
      
      if(ebic_m > ebic_t){
        
        m = m + 1  
      
      }else if(ebic_n < ebic_m & m > 2){
        
        m = m - 1
        
      }else{
        
        A_m_S <- A_m
        
        break
      
      }
       
    }
    
    Splicing_Set <- Set_C_h( h = d - m , s = A_m_S, X = X , Y = Y , family = family)
      
    Candidate_Set <- sort(c(A_m_S,Splicing_Set))
    
    Beta_new <- rep(0,dim(X)[2])
    
    Beta_new[Candidate_Set] <- coefficients(glm.fit(x = X[,Candidate_Set],y = Y, family = family))
    
    IHT_res <- IHT_w_ucheck(Beta_new, k , Y, X, family, U, U_rate, tol, max_iter, p)
    
    E_temp <- ebic(IHT_res$Beta_s, Y, X, family)
    
    ID_out <- cbind(ID_out, c(IHT_res$ID_None_s,rep(0,dim(ID_out)[1]-length(IHT_res$ID_None_s))))

    ebic_out <- c(ebic_out, E_temp)

    if(E_temp < E ){ 
      
      IHT_res_out <- IHT_res
      
      E <- E_temp
      
      end_count <- 0
      
    }
    
    end_count <- end_count+1
    
    if(end_count == 3) { break }
    
    k <- floor((m+k)/2)
  
    
  }
  
  fit<-list(ebic_out = ebic_out,ID_out = ID_out, IHT_res_out = IHT_res_out )
  
}

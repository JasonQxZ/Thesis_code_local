# This is the source code for "The Sparse MLE for Ultrahigh-Dimensional Feature Screening" 
# Journal of the American Statistical Association, 109(507), 1257–1269.
# Authors: Chen Xu and Jiahua Chen

#####################################################################################

# Hard thresholding function--------------------------------

Hard<-function(t=c(0,0), k=-1, lam=1)
{  y<-t
   t<-abs(t)

   if(k > 0)
   { lam<-sort(t, decreasing=TRUE)[k] }

   y[t<lam]<-0
   return(y)
}




# A: p*n,  a: n*1-------------------------------

Fmt<-function(A, a, rr=5)
{
  pp<-dim(A)[1];  r<-pp/rr

  b<- matrix(0, nrow=pp, ncol=1)

  for(i in 1:rr )
  {
   b[((i-1) * r +1) : (i*r)] <- A[ ((i-1) * r +1) : (i*r), ] %*% a
  }

return(b)
}



# log-likelihood function-------------------------

lh<-function(Y, X, beta, sig=1, tag="N")
{
  X<-as.matrix(X)
  dx<-dim(X);  n<-dx[1];  p<-dx[2]

  ind_0<-(1:p)[beta != 0]

  bc_0<- beta[ind_0]
  bc_0<-as.matrix(bc_0, ncol=1)

  Xs_0<-X[, ind_0]

  R_0<- Xs_0 %*% bc_0

  if(tag=="N")
  { R_1 <- R_0^2 /2    }

  if(tag=="B")
  { R_1 <- log(1 + exp(R_0) ) }

  if(tag=="P")
  { R_1 <- exp(R_0) }

  ll<- sum( (Y*R_0 -  R_1)/sig^2 )

  if(tag=="P")
  ll<- sum( dpois(Y, R_1, log = TRUE) )

  return(ll)

}



# function for u check------------------------------

Uh<-function(A, uh, b0, b1, tag="N")
 {
   if(tag != "P")
  {
   if(tag=="N")
   { rho<-1 }

   if(tag=="B")
   { rho<- 0.25 }

    bb<-b1-b0
    in_b<- (1:length(bb))[bb != 0]

    bt<-0

    if(length(in_b) != 0)
    {
     bt<- bb[in_b]
     bt<- as.matrix(bt, ncol=1)
     }

    tm1 <- sum(bt^2) / uh
    tm2 <- sum((A[,in_b]%*%bt)^2) * rho

    Bt <- exp(tm1 - tm2)
   }

   if(tag=="P")
    {
     theta_0<-0; theta_1<-0

     in_b0<- (1:length(b0))[b0 != 0]
     in_b1<- (1:length(b1))[b1 != 0]

      if( length(in_b0) !=0)
     {
      bt0<- b0[in_b0];  bt0<- as.matrix(bt0, ncol=1)
      theta_0<- A[,in_b0] %*% bt0
      }

      if( length(in_b1) !=0)
      {
       bt1<- b1[in_b1];  bt1<- as.matrix(bt1, ncol=1)
       theta_1<- A[,in_b1] %*% bt1
      }

      theta_t<- pmax(theta_0, theta_1)

      bb<- b1-b0
      in_b<- (1:length(bb))[bb != 0]

      bt<-0

      if( length(in_b) != 0)
      { bt<- bb[in_b];  bt<- as.matrix(bt, ncol=1) }

      tm1 <-  sum(bt^2) / uh
      tm2 <-  sum( exp(theta_t) * (theta_0 - theta_1)^2 )

      Bt <- exp(tm1 - tm2)
    }

    return(Bt)
 }


###################################
## X: n*p;  Y: n*1---------------
###################################

AHT<-function(Y, X, beta0, k=20, rr=5, T=500, U=1, er=10^(-3), tag="N")
{

 # initializing------------------------------

 dX<-dim(X);  n<-dX[1]; p<-dX[2]
 X_t<-t(X);

 ind_0<- (1:p)[beta0 != 0]

 # Case if beta0==0--------------------------

 if(length(ind_0) == 0 )
 {
   R_0<-matrix(0, ncol=1, nrow=n)

   if(tag=="B")
   { R_0<- exp(R_0) / (1+exp(R_0) ) }

    if(tag=="P")
   { R_0<- exp(R_0) }


   V_0<-Fmt(A= X_t, a=Y - R_0, rr=rr)
   beta1<- Hard(t=V_0, k=k)

   ind_1<- (1:p)[beta1 != 0]
   Xs_1<-X[, ind_1]
   tau_1 <- svd(Xs_1)$d[1]

   ############################
    uu<-  U/ tau_1^2
   #############################
  }

 # Case if beta0!=0-----------------------

 if(length(ind_0) != 0 )
 {
   bc_0<- beta0[ind_0]
   bc_0<- as.matrix(bc_0, ncol=1)

   Xs_0<-X[, ind_0]

   R_0<- Xs_0 %*% bc_0

   if(tag=="B")
   { R_0<- exp(R_0) / (1+exp(R_0) ) }

   if(tag=="P")
   { R_0<- exp(R_0) }

    V_0<-Fmt(A= X_t, a=Y - R_0, rr=rr)

   tau_0 <- svd(Xs_0)$d[1]

   ############################
    uu<-  U/ tau_0^2
   
  
   #############################
 }


# iteration start---------------------------
   u_c <- uu
   i<-1

repeat{

   repeat
      {
        gamma0<- beta0 + uu * V_0
        beta1<- Hard(t=gamma0, k=k)

        ########## u-check #################
         ucheck<- Uh(A=X, uh=uu, b0=beta0, b1=beta1, tag=tag)

         if( ucheck >= 1 )
         {break}else{ uu <- 0.5 * uu}

        ####################################

      }

      ######## convergence check ##############

        MSE<- sqrt(sum((beta1-beta0)^2))
        if( (i>T) ||  (MSE < er) ){break}

      #########################################

      beta0<-beta1
      ind_0<- (1:p)[beta0 != 0]

      bc_0<- beta0[ind_0]
      bc_0<-as.matrix(bc_0, ncol=1)
      Xs_0<-X[, ind_0]

      R_0<- Xs_0 %*% bc_0

      if(tag=="B")
      { R_0<- exp(R_0) / (1+exp(R_0) ) }

      if(tag=="P")
      { R_0<- exp(R_0) }

      V_0<-Fmt(A= X_t, a=Y - R_0, rr=rr)

      uu<-u_c
      i<-i+1

  }

 index<- (1:p)[beta1 != 0 ]

list(index=index, B=bc_0, step=i)
}


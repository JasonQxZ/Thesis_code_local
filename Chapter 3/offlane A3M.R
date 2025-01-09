
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
  n <- nrow(data)
  card_U <- nrow(Y)
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

sig <- function(X,Y,mode,f_i){
  
  d_B <- dep_an(X[,mode == 1], Y)
  mode[f_i] <- 0
  d_F <- dep_an(X[,mode == 1], Y)
  
  s <- (d_B - d_F) / d_B
  return(s)
}

non_signf <- function(X,Y,mode,i){
  
  B <- rep(0,length(mode))
  R <- mode
  Ts <- mode
  Ts[i] <- 0
  
  indexs <- which(Ts == 1)
  Num <- length(indexs)
  A <- sample(Num)
  
  for(i in 1:Num){
    rnd <- A[i]
    ind <- indexs[rnd]
    if(sig(X,Y,R,ind) == 0){
      B[ind] <- 1
      R[ind] <- 0
    }
  }
  index_del <- B == 1
  return(index_del)
}

OFS_A3M <- function(X, Y){
  Y <- as.matrix(Y,ncol=1)
  start <- Sys.time()
  p <- dim(X)[2]
  
  mode <- rep(0,p)
  dep_Mean <- 0
  dep_Set <- 0
  depArray <- rep(0,p)
  
  for(i in 1:p){
    col <- as.matrix(X[,i],ncol=1)
    dep_single <- dep_an(col,Y)
    depArray[i] <- dep_single
    if(dep_single > dep_Mean){
      mode[i] <- 1
      cols <- as.matrix(X[,mode == 1],ncol= sum(mode==1))
      dep_New <- dep_an(cols,Y)
      
      if(dep_New > dep_Set){
        dep_Set <- dep_New
        dep_Mean <- dep_Set/sum(mode)
      } else if(dep_New == dep_Set){
        index_del <- non_signf(X,Y,mode,i)
        mode[index_del] <- 0
        dep_Mean <- dep_Set/sum(mode)
      } else {
        mode[i] <- 0
      }
    }
  }
  
  selectedFeatures <- which(mode == 1)
  
  end <- Sys.time()
  time <- end - start
  
  return(list(selectedFeatures = selectedFeatures, time = time))
}
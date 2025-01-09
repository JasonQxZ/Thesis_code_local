card_density <- function(sets, Y, label, N) {
  D <- sort(sets, index.return = TRUE)
  cTotal <- 0
  cNum <- 0
  last_density <- 0
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

dep_an <- function(data, Y) {
  n <- nrow(data)
  card_U <- nrow(Y)
  card_ND <- 0
  D <- stats::dist(data, method = "euclidean")
  DArray <- as.matrix(D)
  for (i in 1:n) {
    d <- DArray[, i]
    class <- Y[i, 1]
    card_ND <- card_ND + card_density(d, Y, class, n)
  }
  dep <- card_ND / card_U
  return(dep)
}

sig <- function(X, Y, mode, f_i) {
  d_B <- dep_an(X[,mode == 1, drop = FALSE], Y)
  mode[f_i] <- 0
  d_F <- dep_an(X[,mode == 1, drop = FALSE], Y)
  s <- (d_B - d_F)
  return(s)
}

non_signf <- function(X, Y, mode, i) {
  B <- numeric(length(mode))
  R <- mode
  T <- mode
  T[i] <- 0
  indexs <- which(T == 1)
  Num <- length(indexs)
  A <- sample(Num)
  for (i in 1:Num) {
    rnd <- A[i]
    ind <- indexs[rnd]
    if (sig(X, Y, R, ind) <= 0) {
      B[ind] <- 1
      R[ind] <- 0
    } else {
      R[ind] <- 1
    }
  }
  index_del <- which(B == 1)
  return(index_del)
}

OFS_Density <- function(X, Y) {
  p <- ncol(X)
  mode <- numeric(p)
  dep_Mean <- 0
  dep_Set <- 0
  depArray <- numeric(p)
  start <- Sys.time()
  for (i in 1:p) {
    col <- X[, i]
    dep_single <- dep_an(matrix(col, ncol=1), Y)
    depArray[i] <- dep_single
    if (dep_single >= dep_Mean) {
      mode[i] <- 1
      cols <- X[,mode == 1, drop = FALSE]
      dep_New <- dep_an(cols, Y)
      if (dep_New > dep_Set) {
        dep_Set <- dep_New
        dep_Mean <- sum(depArray[mode == 1]) / sum(mode)
        next
      }
      if (abs(dep_New - dep_Set) / dep_Set <= 0.05) {
        index_del <- non_signf(X, Y, mode, i)
        mode[index_del] <- 0
        dep_Del <- dep_an(X[,mode == 1, drop = FALSE], Y)
        if (dep_Del >= dep_Set) {
          dep_Set <- dep_Del
          dep_Mean <- sum(depArray[mode == 1]) / sum(mode)
        } else {
          mode[index_del] <- 1
        }
        next
      }
      mode[i] <- 0
    }
  }
  selectedFeatures <- which(mode == 1)
  time <- Sys.time() - start
  return(list(selectedFeatures = selectedFeatures, time = time))
}

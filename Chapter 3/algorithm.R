library(methods)
source("./Streaming_methods.R")
source("./tool.R")
library(SMLE)
library(SIS)
library(abess)
###############################################################################
## Data Class Definition 
###############################################################################

setClass("Streaming_Data", slots= c( name = "character", 
                                     X = "matrix",
                                     y = "numeric",
                                     causal_index = "numeric",
                                     train_index = "numeric"),
         prototype = list(name = NA_character_,
                          X = matrix(),
                          y = NA_real_,
                          causal_index = NA_real_,
                          train_index = NA_real_))
## Data Methods Definition

setGeneric("setData_", function(x,...) standardGeneric("setData_"))
setMethod("setData_", "Streaming_Data", function(x,train_index,...){
  x@train_index <- train_index
  x
})

#Data<-SMLE::Gen_Data()
#raw_data<-new("Streaming_Data" , X =Data$X, y =Data$Y , causal_index = Data$subset_true)
#processed_data<- setData_(raw_data,train_index = 1:60)

###############################################################################
## Result Class Definition
###############################################################################




###############################################################################
## Algorithm Class Definition
###############################################################################
setClass("Algorithm", slots = c(Processed_Data = "Streaming_Data",
                               Methods = "list",
                               Memory = "list",
                               Result = "list",
                               shuffle_order = "integer"),
         prototype = list(Processed_Data = new("Streaming_Data"),
                          Methods = list(),
                          Memory = list(),
                          Result = list(),
                          shuffle_order = NA_integer_))

setMethod("initialize", "Algorithm", function(.Object,...){
  .Object <- callNextMethod()
  
  # initialize the result matrix by the number of methods.
  .Object@Result <-vector(mode = "list", length = length(.Object@Methods))
  .Object@Memory <-vector(mode = "list", length = length(.Object@Methods))
  
  validObject(.Object)
  .Object
})

#a1 <- new("Algorithm", Methods = list("alpha_investing","Iter_SMLE"), Processed_Data = processed_data)

setGeneric("run",function(x,...){standardGeneric("run")})
setMethod("run", "Algorithm", function(x, shuffle = FALSE, s =20, k =25, family = gaussian(), ...){
  
############################################################################### 
  
  #Initialize
  n <- dim(x@Processed_Data@X)[1]
  
  p <- dim(x@Processed_Data@X)[2]
  
  size_Gt <- s
  
###############################################################################
  
  train_index <- x@Processed_Data@train_index
   
  if(is.na(train_index[1])){
    
    train_x <- x@Processed_Data@X
    
    train_y <- x@Processed_Data@y
    
    test_x <- NULL
    
    test_y <- NULL
    
  }else{
    
    train_x <- x@Processed_Data@X[train_index,]
    
    train_y <- x@Processed_Data@y[train_index]
    
    test_index <- (1:n)[!(1:n) %in% train_index]
    
    test_x <- x@Processed_Data@X[test_index,]
    
    test_y <- x@Processed_Data@y[test_index]
     
  }
  
###############################################################################
  
  if(shuffle){
  
    shuffle_order <- sample(1:p,p)
    
    train_x <- train_x[,shuffle_order]
    
    test_x <- test_x[,shuffle_order]
    
  }else{
    
    shuffle_order <- 1:p 
    
  }
###############################################################################
# Initialize Momory 
################################################################################
  
  for(j in 1:length(x@Methods)){
    
    # Public
    
    x@Memory[[j]]$time <- 0 
    
    x@Memory[[j]]$St <- NULL 
    
    x@Memory[[j]]$num_feature_used <- 0
    
    # For SMLE
    
    x@Memory[[j]]$coefficients <- NULL
    
    x@Memory[[j]]$residual<- NULL
    
    # For alpha-investing
    
    x@Memory[[j]]$wealth <- 0.5
    
    x@Memory[[j]]$alpha <- 0.5
    
  }
  
  feature_remain <- p
  
  train_remain <- train_x
  
  while( feature_remain >= 0){
    
    if(feature_remain > size_Gt){Gt <- matrix(train_remain[,1:size_Gt],ncol=s)}else{Gt <- matrix(train_remain,ncol =feature_remain )}
    
    for(j in 1:length(x@Methods)){
      
      METD <- x@Methods[[j]]
      
      # Update Active set
      
      x@Memory[[j]]<- switch(METD , 
             Iter_SMLE = Iter_SMLE(x@Memory[[j]],Gt,train_y, k,family),
             Iter_SIS = Iter_SIS(x@Memory[[j]],Gt,train_y, k,family),
             Iter_abess = Iter_abess(x@Memory[[j]],Gt,train_y, k,family),
             alpha_investing = Group_alpha_investing(x@Memory[[j]], Gt,train_y),
             BANS = BANS(x@Memory[[j]], Gt, train_y, k, family),
             BANS_l2 = BANS(x@Memory[[j]], Gt, train_y, k,family,l=2),
             BANS_l5 = BANS(x@Memory[[j]], Gt, train_y, k,family,l=5),
             BANS_l10 = BANS(x@Memory[[j]], Gt, train_y, k,family,l=10),
             BANS_l15 = BANS(x@Memory[[j]], Gt, train_y, k,family,l=15),
             BANS_l20 = BANS(x@Memory[[j]], Gt, train_y, k,family,l=20),
             Streaming_SIS = Streaming_SIS(x@Memory[[j]], Gt, train_y, k),
             OFS_fisher = OFS_fisher(x@Memory[[j]], Gt, train_y, alpha = 0.05),
             Saola_mi = Saola_mi(x@Memory[[j]], Gt, train_y,0),
             Saola_z = Saola_z(x@Memory[[j]], Gt, train_y),
             OSFS_FI = OSFS_FI(x@Memory[[j]], Gt, train_y, 0.25),
             OSFS_A3M = OSFS_A3M(x@Memory[[j]], Gt, train_y),
             OSFS_Density = OSFS_Density(x@Memory[[j]], Gt, train_y)
             
      )
      
      x@Result[[j]] <- c(x@Result[[j]],list( time = x@Memory[[j]]$time, subset_index = x@Memory[[j]]$St_index))
    
    }
    
    feature_remain <- feature_remain - size_Gt
    
    if(feature_remain <= 0) {break}
    
    train_remain <- train_remain[,-(1:size_Gt)]
  
  }
  x@shuffle_order <- shuffle_order
  x
  
})
#R <- run(a1)


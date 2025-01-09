
######################################################################################
# This script provides R code that produces the numerical results shown in the paper.#
######################################################################################
# Last update: July 9 2024
# Require R > 4.0.0

# Packages needed to run this code:
# ----------------------------------------------------------------------------------------
## Package glmnet is needed to construct Figure 2 and Table 4.
## Package SIS is needed to construct Figures 7-9 and Table 4.
## Package qqman is needed to construct Figures 6-8. 
## Package SMLE is the main package to implement joint feature screening in most examples.
## Package abess is needed to construct Table 4
## Remark: using lower-version packages may lead to results different from the paper.

library("glmnet")                 # version 4.1.7
library("SIS")                    # version 0.8.8
library("qqman")                  # version 0.1.9
library("SMLE")                   # version 2.1.1
library("abess")                  # version 0.4.7
library("VariableScreening")      # version 0.2.1

###########################
# Code for Section 3      #
###########################

# Section 3.2---Demo code---Simulating ultrahigh-dimensional GLM data---------------------

set.seed(1)
Data_ctg <- Gen_Data(n = 200, p = 1000, family = "gaussian", pos_ctgidx = c(1, 2, 3), 
                     level_ctgidx = c(3, 4, 5))
head(Data_ctg$X)[, 1:5]


# Section 3.2---Demo code---Joint feature screening--------------------------------

fit <- SMLE(Y = Data_ctg$Y, X = Data_ctg$X, k = 15, family = "gaussian", keyset = c(1, 4, 5), 
            categorical = TRUE, group = TRUE)
fit

# Section 3.2---Demo code---Post-screening Selection-------------------------------

fit_s <- smle_select(fit, criterion = "ebic", gamma_seq = seq(0,1,0.2), vote = TRUE)
fit_s


# Section 3.2---Figure 2-----------------
# Figure 2 is constructed based on NumIter=100 repetitions and it 
# will have a long running time. Set NumIter=3 for checking purposes.

NumIter <- 100

# "OSMLE" is the original IHT implementation used in Xu and Chen (2014).
source("AHT_Xu2014.r") 

OSMLE <- function(X,Y,family,k){  
  fit_pre <- glmnet(x = X,y = Y,family = family) 
  Beta0 <- c(fit_pre$beta[,dim(fit_pre$beta)[2]])
  tag <- switch(family, "gaussian" = "N" ,"binomial" = "B","poisson" = "P")
  AHT(Y = Y, X = X, beta0 = Beta0, k = k, U = 1, rr = 1, T = 100, tag = tag ,er = 10^(-3))
}

# _g and _b stands for the running time in linear model and logistic model.
SMLE_time_g <- matrix(0, nrow = 9 , ncol = NumIter)
IHT_time_g <- matrix(0, nrow = 9 , ncol = NumIter)

SMLE_time_b <- matrix(0, nrow = 9 , ncol = NumIter)
IHT_time_b <- matrix(0, nrow = 9 , ncol = NumIter)

p <- seq(from = 4000,to = 20000, length.out = 9)
n <- 1000
k <- 10

for(j in 1 : NumIter){
  
  for(i in 1 : 9){
    
    family = "gaussian" 
    
    DD <- Gen_Data(n, p[i], family = family)
    Y <- DD$Y
    X <- DD$X
    
    iT <- proc.time()
    OSMLE(X = X, Y = Y, k = k,family = family)
    IHT_time_g[i, j] <- (proc.time() - iT)[3]
    
    
    sT <- proc.time()
    SMLE(X = X, Y = Y, k = k, family = family, intercept = FALSE, standardize = FALSE,
         categorical = FALSE, max_iter = 100, tol = 10^(-3))
    SMLE_time_g[i,j] <- (proc.time()-sT)[3]
    
    
  }  
  
}

# linear model comparison

IHT_average_time <- apply(IHT_time_g,1, mean)
plot(x = p,y = IHT_average_time,type = "b", pch = 18, col = "blue", cex.main=1.5,cex.lab=1.5,
     xlab="The number of features", ylab="Elapsed Time(s)",main = "Linear model",
     ylim = c(0,1.2*max(IHT_average_time)))
lines(x = p, y = apply(SMLE_time_g, 1, mean), pch = 19, col = "red", type = "b", lty = 2)
legend("topleft",legend = c("Original Code", "SMLE()"),
       pch = 18 : 19,lty = 1 : 2,  col = c("blue","red"))


# logistic model comparison


for(j in 1 : NumIter){
  
  for(i in 1 : 9){

    DD <- Gen_Data(n = n, p = p[i], num_truecoef = 10,family = "binomial",
                   effect_truecoef = rep(1.5, 10), rho = 0.3)
    Y <- DD$Y; X <- DD$X
    
    iT <- proc.time()
    OSMLE(X = X,Y = Y, k = k,family = "binomial")
    IHT_time_b[i, j]<-(proc.time()-iT)[3]
    sT<-proc.time()
    
    SMLE(X = X, Y = Y, k = k, family = "binomial" , intercept = F, standardize = FALSE,
         categorical = FALSE, max_iter = 100,tol = 10^(-3), U = 4)
    SMLE_time_b[i,j] <- (proc.time()-sT)[3]
    
  }  
  
}

IHT_average_time <- apply(IHT_time_b, 1, mean)
plot(x = p,y = IHT_average_time,type = "b", pch = 18, col = "blue", 
     xlab = "The number of features", ylab = "Elapsed Time(s)",main = "Logistic model",
     ylim = c(0,1.2 * max(IHT_average_time)), cex.main=1.5,cex.lab=1.5)
lines(x = p,y = apply(SMLE_time_b, 1, mean), pch = 19, col = "red", type = "b", lty = 2)
legend("topleft",legend = c("Original Code", "SMLE()"),
        pch = 18 : 19,lty = 1 : 2,  col = c("blue", "red"))

# Section 3.2---Figures 3 & 4 & 5 --------

set.seed(1)

Data_sim <- Gen_Data(n = 200, p = 1000, correlation = "AR", rho = 0.9, sigma = 1,
                     family = "gaussian", pos_truecoef = c(1,3,5,7,9),
                     effect_truecoef = (0.8)*c(2,-3,-3,2,-4))

fit_path <- SMLE(Y = Data_sim$Y, X = Data_sim$X , k = 10 , coef_initial = rep(0,1000))
plot(fit_path,cex = 0.8,cex.axis = 1.5, cex.lab = 1.5)

# The following comment is a reminder that the object fit_s is defined previously.
# fit_s <- smle_select(fit, criterion = "ebic", gamma_seq = seq(0,1,0.2), vote = TRUE)
plot(fit_s,cex.axis = 1.5, cex.lab = 1.5)



###########################
# Code for Section 4      #
###########################

# Section 4.1---Demo Code for SMLE-screening---

set.seed(1)

Data_eg<-Gen_Data(n = 400, p = 1000, family = "binomial",  correlation="AR", rho = 0.9,
                  pos_truecoef = c(1,3,5,7,9),effect_truecoef = c(2,3,-3,3,-4))

print(Data_eg)

fit1<-SMLE(Y = Data_eg$Y, X = Data_eg$X,family = "binomial", k = 10, coef_initial = rep(0,1000))
summary(fit1)

coef(summary(glm(Data_eg$Y ~ Data_eg$X[ , 2], family = "binomial")))


fit1_s <- smle_select(fit1, criterion = "ebic")
fit1_s

summary(fit1_s)

#
# Section 4.1---Demo Code for SMLE-screening---

set.seed(1)

# Function to run SMLE with different initial values
run_AR_smle_simulation <- function(initialization_method ) {
  # Generate example data
  Data_eg <- Gen_Data(n = 400, p = 1000, family = "binomial", correlation = "AR", rho = 0.9,
                      pos_truecoef = c(1, 3, 5, 7, 9), effect_truecoef =  c(2, 3, -3, 3, -4))
  
  # Set initial coefficients based on the method
  if (initialization_method == "default") {
    coef_initial <- NULL # Assuming the default method will handle this
  } else if (initialization_method == "zero") {
    coef_initial <- rep(0, 1000)
  }
  
  # Fit the SMLE model
  fit <- SMLE(Y = Data_eg$Y, X = Data_eg$X, family = "binomial", k = 10, coef_initial = coef_initial)
  
  return(fit)
}

# Section 4.2---Demo Code for SMLE-screening---

run_ID_smle_simulation <- function(initialization_method ) {
  # Generate example data
  Data_eg <- Gen_Data(n = 300, p = 1000, family = "binomial", correlation = "ID", rho = 0.9,
                      pos_truecoef = c(1, 3, 5, 7, 9), effect_truecoef =  c(2, 3, -3, 3, -4))
  
  # Set initial coefficients based on the method
  if (initialization_method == "default") {
    coef_initial <- NULL # Assuming the default method will handle this
  } else if (initialization_method == "zero") {
    coef_initial <- rep(0, 1000)
  }
  
  # Fit the SMLE model
  fit <- SMLE(Y = Data_eg$Y, X = Data_eg$X, family = "binomial", k = 10, coef_initial = coef_initial)
  
  return(fit)
}

# Function to summarize the results of the simulation
summarize_results <- function(results_list) {
  # Example: summarizing the coefficients
  SSR <- sapply(results_list, function(fit) floor(sum(Data_eg$subset_true %in% fit$ID_retained)/5))
  PSR <- sapply(results_list, function(fit) sum(Data_eg$subset_true %in% fit$ID_retained)/5)
  return(list(SSR = SSR, PSR = PSR))
}

# Number of simulations
n_simulations <- 100

# Run simulations for different initialization methods
results_default <- replicate(n_simulations, run_AR_smle_simulation("default"), simplify = FALSE)
results_zero <- replicate(n_simulations, run_AR_smle_simulation("zero"), simplify = FALSE)

# Summarize results
summary_default <- summarize_results(results_default)
summary_zero <- summarize_results(results_zero)

# Combine results into a data frame for plotting
results <- data.frame(
  Initialization = rep(c("Default", "Zero"), each = 100),
  SSR = c(summary_default$SSR, summary_zero$SSR),
  PSR = c(summary_default$PSR, summary_zero$PSR)
)

#Calculate mean SSR and PSR for each Initialization method
mean_results_AR <- aggregate(cbind(SSR, PSR) ~ Initialization, data = results, FUN = mean)

# Create a bar plot
barplot_heights <- as.matrix(mean_results_AR[, 2:3])
barplot_names_AR <- mean_results_AR$Initialization

# Plot the barplot
barplot(
  t(barplot_heights),
  beside = TRUE,
  names.arg = barplot_names_AR,
  ylab = "Rate",
  cex.main = 1.5,  # Increase main title size
  cex.lab = 1.5,   # Increase label size
  ylim = c(0, 1),  # Set y-axis range from 0 to 1
  legend.text = c("SSR", "PSR"),
  args.legend = list(x = "topright", cex = 1.5)  # Increase legend text size
)
results_default <- replicate(n_simulations, run_ID_smle_simulation("default"), simplify = FALSE)
results_zero <- replicate(n_simulations, run_ID_smle_simulation("zero"), simplify = FALSE)

# Summarize results
summary_default <- summarize_results(results_default)
summary_zero <- summarize_results(results_zero)

# Combine results into a data frame for plotting
results <- data.frame(
  Initialization = rep(c("Default", "Zero"), each = 100),
  SSR = c(summary_default$SSR, summary_zero$SSR),
  PSR = c(summary_default$PSR, summary_zero$PSR)
)

#Calculate mean SSR and PSR for each Initialization method
mean_results_ID <- aggregate(cbind(SSR, PSR) ~ Initialization, data = results, FUN = mean)

# Create a bar plot
barplot_heights_ID <- as.matrix(mean_results_ID[, 2:3])
barplot_names_ID <- mean_results_ID$Initialization

# Plot the barplot
# Section 4.2---Figures 6 --------
barplot(
  t(barplot_heights_ID),
  beside = TRUE,
  names.arg = barplot_names_ID,
  ylab = "Rate",
  cex.main = 1.5,   # Increase main title size
  cex.lab = 1.5,    # Increase label size
  ylim = c(0, 1),   # Set y-axis range from 0 to 1
  legend.text = c("SSR", "PSR"),
  args.legend = list(x = "topright", cex = 1.5)  # Increase legend text size
)

#Section 4.2 --- Table 4 ---------------------------

# Table 4 is summarized based on NumIter=500 repetitions, which takes a long time
# to run. Set NumIter=5 for quick checking purposes.
# SIS/ISIS is implemented by Package "SIS".
# Lasso is implemented by package "glmnet".

# We have suppressed warnings from GLMNET so that the output file is not too large.
# The warning message from GLMNET package is due to the fact that Lasso may  
# not always be able to produce an estimate with sparsity exactly matching the
# pre-specified screening size k. When it happens, a Lasso estimate with the 
# largest sparsity not exceeding k is returned. This is expected and the warning
# message is safe to ignore.

# We have suppressed routine iteration output and warnings from ISIS being 
# printed to screen since the resulting output file would be much too large. 
# The warnings relate to ISIS fitting a saturated model periodically with Poisson
# data. However, SIS, glmnet, SMLE and abess do not give warnings with this data.

# The abess function modifies the seed of the random number generator in a 
# way that affects the ability to reproduce results exactly. For this reason, we 
# randomly create a vector of seeds with length equal to the number of iterations
# and we re-seed the random number generator each iteration. 

NumIter <- 500
num_methods = 6


###########Linear Regression ####################################
###########Linear Regression ####################################
set.seed(1)
seed.vec <- round(sample(c(-1,1), size=NumIter, replace=TRUE)*runif(NumIter,min=0,max=1e7))

Lasso_avg_length<-rep(1 : NumIter)
PRR<-matrix(0, nrow = num_methods , ncol = NumIter)
SSR<-matrix(0, nrow = num_methods , ncol = NumIter)
TIME<- matrix(0, nrow = num_methods , ncol = NumIter)
time<-rep(0, num_methods )  
k = 20

for(j in 1 : NumIter){
  
  set.seed(seed.vec[j])
  Data<- Gen_Data(n = 100, p = 2000, pos_truecoef = c(1,2,3,4), family = "gaussian"
                  , effect_truecoef = rep(2.5, 4), correlation = "CS", rho = 0.3)
  ##models fitting
  #SIS
  sis_a <- proc.time()
  SIS <- SIS(x = Data$X, y = Data$Y, family = "gaussian", iter = F, nsis = k)
  sis_b <- proc.time()
  time[1] <- (sis_b - sis_a)[3]
  
  #ISIS
  isis_a <- proc.time()
  suppressWarnings(capture.output(ISIS <- SIS(x = Data$X, y = Data$Y, family = "gaussian", iter = T, nsis = k)))
  isis_b <- proc.time()
  time[2] <- (isis_b - isis_a)[3]
  
  #Lasso
  Lasso_a <- proc.time()
  suppressWarnings(Lasso_fit <- glmnet(x = Data$X, y = Data$Y, family = "gaussian", pmax = k))
  Lasso_b <- proc.time()
  time[3] <- (Lasso_b - Lasso_a)[3]
  
  #abess
  abess_a <- proc.time()
  abess_fit <- abess(x = Data$X, y = Data$Y, family = "gaussian", support.size = k)
  abess_b <- proc.time()
  time[4]<-(abess_b - abess_a)[3]
  
  #VariableScreening
  Vscreening_a <- proc.time()
  Vscreening_fit <- screenIID(X= Data$X, Y = Data$Y)
  Vscreening_b <- proc.time()
  time[5]<-(Vscreening_b - Vscreening_a)[3]
  
  
  #SMLE
  SMLE_a <- proc.time()
  SMLE_fit <- SMLE(X = Data$X, Y = Data$Y, family = "gaussian", standardize = F, categorical = FALSE, k = k, fast = F)
  SMLE_b <- proc.time()
  time[6] <- (SMLE_b - SMLE_a)[3]
  
  Lasso_index <- (1 : dim(Data$X)[2])[Lasso_fit$beta[, dim(Lasso_fit$beta)[2]] != 0]
  Lasso_avg_length[j] <- length(Lasso_index)
  ind <- extract(abess_fit)$support.vars
  ind <- as.numeric(regmatches(ind, gregexpr("[[:digit:]]+", ind)))
  
  rank_id <- order(Vscreening_fit$rank)[1:k]
  index_set <- list(SIS$ix0, ISIS$ix0, Lasso_index, ind,rank_id, SMLE_fit$ID_retained)
  
  PRR[,j] <- unlist(lapply(1 : num_methods ,function(i){
    
    sum(Data$subset_true %in% index_set[[i]])
    
  }))
  SSR[,j] <- unlist(lapply(1:num_methods ,function(i){
    
    sum(Data$subset_true %in% index_set[[i]]) == 4
    
  }))
  TIME[,j] <- time
}


# Order: SIS / ISIS / glmnet / abess / SMLE / SMLE_fast 

# Screening Accuracy
PRR_linear <- apply(PRR, 1, mean)/4
SSR_linear <- apply(SSR, 1, mean)
# Computational Time
Time_linear <- apply(TIME, 1, mean)


###########Logistic Regression ####################################
set.seed(1)
seed.vec <- round(sample(c(-1,1), size=NumIter, replace=TRUE)*runif(NumIter,min=0,max=1e7))

Lasso_avg_length <- rep(1:NumIter)
PRR <- matrix(0, nrow = num_methods , ncol = NumIter)
SSR <- matrix(0, nrow = num_methods , ncol = NumIter)
TIME <- matrix(0, nrow = num_methods , ncol = NumIter)
time <- rep(0, num_methods )  
k=30
for(j in 1:NumIter){
  
  set.seed(seed.vec[j])
  Data<- Gen_Data(n = 600, p = 2000, pos_truecoef = c(1,2,3,4), family = "binomial",
                  effect_truecoef = rep(1.1, 4),correlation = "CS", rho = 0.3)
  ##models fitting
  #SIS
  sis_a <- proc.time()
  SIS <- SIS(x = Data$X, y = Data$Y, family = "binomial", iter = F, nsis = k)
  sis_b <- proc.time()
  time[1] <- (sis_b - sis_a)[3]
  
  #ISIS
  isis_a <- proc.time()
  suppressWarnings(capture.output(ISIS <- SIS(x = Data$X, y = Data$Y, family = "binomial", iter = T, nsis = k)))
  isis_b <- proc.time()
  time[2] <- (isis_b - isis_a)[3]
  
  #Lasso
  Lasso_a <- proc.time()
  suppressWarnings(Lasso_fit <- glmnet(x = Data$X, y = Data$Y, family = "binomial", pmax = k))
  Lasso_b <- proc.time()
  time[3] <- (Lasso_b - Lasso_a)[3]
  
  #abess
  abess_a <- proc.time()
  abess_fit <- abess(x = Data$X, y = Data$Y, family = "binomial", support.size = k)
  abess_b <- proc.time()
  time[4]<-(abess_b - abess_a)[3]
  
  #VariableScreening
  Vscreening_a <- proc.time()
  Vscreening_fit <- screenIID(X= Data$X, Y = Data$Y)
  Vscreening_b <- proc.time()
  time[5]<-(Vscreening_b - Vscreening_a)[3]
  
  #SMLE
  SMLE_a <- proc.time()
  SMLE_fit <- SMLE(X = Data$X, Y = Data$Y, family = "binomial", standardize = F, categorical = FALSE, k = k, fast = F)
  SMLE_b <- proc.time()
  time[6] <- (SMLE_b - SMLE_a)[3]
  
  Lasso_index <- (1 : dim(Data$X)[2])[Lasso_fit$beta[ , dim(Lasso_fit$beta)[2]] != 0]
  Lasso_avg_length[j] <- length(Lasso_index)
  ind <- extract(abess_fit)$support.vars
  ind <- as.numeric(regmatches(ind, gregexpr("[[:digit:]]+", ind)))
  index_set <- list(SIS$ix0, ISIS$ix0, Lasso_index, ind,rank_id, SMLE_fit$ID_retained)
  
  PRR[,j] <- unlist(lapply(1 : num_methods, function(i){
    
    sum(Data$subset_true %in% index_set[[i]])
    
  }))
  SSR[,j] <- unlist(lapply(1:num_methods, function(i){
    
    sum(Data$subset_true %in% index_set[[i]]) == 4
    
  }))
  TIME[,j] <- time
}


# Order: SIS / ISIS / glmnet / abess / SMLE / SMLE_fast 

# Screening Accuracy
PRR_binomial<-apply(PRR, 1, mean)/4
SSR_binomial<-apply(SSR, 1, mean)
# Computational Time
Time_binomial<-apply(TIME,1, mean)



###########Poisson Regression ####################################
set.seed(1)
seed.vec <- round(sample(c(-1,1), size=NumIter, replace=TRUE)*runif(NumIter,min=0,max=1e7))

Lasso_avg_length <- rep(1 : NumIter)
PRR<-matrix(0, nrow = num_methods, ncol = NumIter)
SSR<-matrix(0, nrow = num_methods, ncol = NumIter)
TIME<- matrix(0, nrow = num_methods, ncol = NumIter)
time<-rep(0, num_methods)  
k=10
for(j in 1:NumIter){
  
  set.seed(seed.vec[j])
  Data<- Gen_Data(n = 250, p = 2000, pos_truecoef = c(1,2,3,4), family = "poisson",
                  effect_truecoef = rep(0.7,4), correlation = "CS", rho = 0.3)
  ##models fitting
  #SIS
  sis_a <- proc.time()
  SIS <- SIS(x = Data$X, y = Data$Y, family = "poisson", iter = F, nsis = k)
  sis_b <- proc.time()
  time[1] <- (sis_b - sis_a)[3]
  
  #ISIS
  isis_a <- proc.time()
  suppressWarnings(capture.output(ISIS <- SIS(x = Data$X, y = Data$Y, family = "poisson", iter = T, nsis = k)))
  isis_b <- proc.time()
  time[2] <- (isis_b - isis_a)[3]
  
  #Lasso
  Lasso_a <- proc.time()
  suppressWarnings(Lasso_fit <- glmnet(x = Data$X, y = Data$Y, family = "poisson", pmax = k))
  Lasso_b <- proc.time()
  time[3] <- (Lasso_b - Lasso_a)[3]
  
  #abess
  abess_a <- proc.time()
  abess_fit <- abess(x = Data$X, y = Data$Y, family = "poisson", support.size = k)
  abess_b <- proc.time()
  time[4]<-(abess_b - abess_a)[3]
  
  #VariableScreening
  Vscreening_a <- proc.time()
  Vscreening_fit <- screenIID(X= Data$X, Y = Data$Y)
  Vscreening_b <- proc.time()
  time[5]<-(Vscreening_b - Vscreening_a)[3]
  
  #SMLE
  SMLE_a <- proc.time()
  SMLE_fit <- SMLE(X = Data$X, Y = Data$Y, family = "poisson", standardize = F, categorical = FALSE, k = k, fast = F)
  SMLE_b <- proc.time()
  time[6] <- (SMLE_b - SMLE_a)[3]
  
  Lasso_index <- (1 : dim(Data$X)[2])[Lasso_fit$beta[,dim(Lasso_fit$beta)[2]] != 0]
  Lasso_avg_length[j] <- length(Lasso_index)
  ind <- extract(abess_fit)$support.vars
  ind <- as.numeric(regmatches(ind, gregexpr("[[:digit:]]+", ind)))
  index_set <- list(SIS$ix0, ISIS$ix0, Lasso_index, ind,rank_id, SMLE_fit$ID_retained)
  
  PRR[,j] <- unlist(lapply(1 : num_methods, function(i){
    
    sum(Data$subset_true %in% index_set[[i]])
    
  }))
  SSR[,j] <- unlist(lapply(1:num_methods, function(i){
    
    sum(Data$subset_true %in% index_set[[i]]) == 4
    
  }))
  TIME[,j] <- time
}


# Order: SIS / ISIS / glmnet / abess / SMLE 

# Screening Accuracy
PRR_poisson <- apply(PRR, 1, mean)/4
SSR_poisson <- apply(SSR, 1, mean)
# Computational Time
Time_poisson <- apply(TIME,1, mean)


#Summary
# Order: SIS / ISIS / glmnet / abess / SMLE 
PRR_linear
PRR_poisson
PRR_binomial

SSR_linear
SSR_poisson
SSR_binomial

Time_linear
Time_poisson
Time_binomial


#Section 4.3 --- Application to high-dimensional genetic data-----------

data("synSNP")
Y_SNP <- synSNP[, 1]
X_SNP <- as.matrix(synSNP[, -1])


#Casual SNPs' positions and effects generating the data was saved in "answers.txt'.
#We use this to point casual SNPs and estimate the screening results.
answers = readLines(paste0(path.package("SMLE"), "/answers.txt"))

level1 = strsplit(answers[1], split = " ")[[1]]
effects1 = c(answers[2], answers[2])
level2 = strsplit(answers[3], split = " ")[[1]]
effects2 = strsplit(answers[4], split = " ")[[1]]
level3 = strsplit(answers[5], split = " ")[[1]]
effects3 = strsplit(answers[6], split = " ")[[1]]
level4 = strsplit(answers[7], split = " ")[[1]]
effects4 = strsplit(answers[8], split = " ")[[1]]

snps_chr = list(level1, level2, level3, level4)


#Section 4.3 ---  figure 7 ---------
# The marginal p-values of each SNPs is pre-calculated and saved in "pvals".
data("pvals")

# Code corresponding to the calculation of p-values
pvals$P <- unlist(lapply(2:10032,function(i){summary(lm(Y_SNP~synSNP[,i],data=synSNP))$coef[2,4]}))

# The first three columns of "pvals" are position information on chromosome, which are 
# provided in the original data set.
qqman::manhattan(pvals,ylim = c(0,8),genomewideline = FALSE)
pvals$logp <- -log10(pvals$P)
pvals$index = rep.int(seq_along(unique(pvals$CHR)), times = tapply(pvals$SNP,pvals$CHR,length))
pvals$pos = pvals$BP
lastbase = 0
for (i in unique(pvals$index)){
  if (i == 1) {
    pvals[pvals$index == i, ]$pos = pvals[pvals$index == i, ]$BP
  } else {
    lastbase = lastbase +max(pvals[pvals$index == (i-1),"BP"])   
    pvals[pvals$index == i,"BP"] = pvals[pvals$index == i,"BP"]-min(pvals[pvals$index == i,"BP"]) +1
    pvals[pvals$index == i, "pos"] = pvals[pvals$index == i,"BP"] + lastbase    
  }
  
}
colors<-c("red","orange","blue","green3")

for(i in 1:4){pvals.highlight = pvals[which(pvals$SNP %in% snps_chr[[i]]), ]
with(pvals.highlight, points(pos, logp, col=colors[[i]], pch = 20))}

#Section 4.3 ---  figure 8 ---------
SMLE_fit <- SMLE(Y = Y_SNP, X = X_SNP, family = "gaussian", k = 60, fast = F)
SIS_fit <- SIS(y = Y_SNP, x = X_SNP, family = "gaussian", nsis = 60, iter = F)

suppressWarnings(Lasso_fit <- glmnet(x = X_SNP,y = Y_SNP, family = "gaussian" , pmax = 60 ))
Lasso_index <- (1 : dim(X_SNP)[2])[Lasso_fit$beta[,dim(Lasso_fit$beta)[2]] != 0]
par(mar = c(5, 6, 4, 2) + 0.1)
# SMLE screening results are showed in red line.
qqman::manhattan(pvals,ylim = c(0,8),genomewideline = FALSE,cex.lab=2,cex.main = 2,main = "SMLE()")
pvals.screened = pvals[SMLE_fit$ID_retained, ]
with(pvals.screened, abline(v = pos, col = "red", pch = 20)) 
for(i in 1 : 4){pvals.highlight = pvals[which(pvals$SNP %in% snps_chr[[i]]), ]
with(pvals.highlight, points(pos, logp, col = colors[[i]], pch = 20))}

# SIS screening results are showed in red line.
qqman::manhattan(pvals,ylim = c(0,8),genomewideline = FALSE, cex.lab=2,cex.main = 2,main = "SIS()")
pvals.screened = pvals[SIS_fit$ix0, ]
with(pvals.screened, abline(v = pos, col = "red", pch = 20)) 
for(i in 1 : 4){pvals.highlight = pvals[which(pvals$SNP %in% snps_chr[[i]]), ]
with(pvals.highlight, points(pos, logp, col = colors[[i]], pch = 20))}

# Lasso screening results are showed in red line.
qqman::manhattan(pvals,ylim = c(0,8),genomewideline = FALSE, cex.lab=2,cex.main = 2,main = "glmnet()")
pvals.screened = pvals[Lasso_index, ]
with(pvals.screened, abline(v = pos, col = "red", pch = 20)) 
for(i in 1 : 4){pvals.highlight = pvals[which(pvals$SNP %in% snps_chr[[i]]), ]
with(pvals.highlight, points(pos, logp, col = colors[[i]], pch = 20))}


# Section 4.3 --- Screening with varying k----

k = c(40,60,80,100)
MD <- list()
for(i in 1 : length(k)){
  SMLE_fit <- SMLE(Y = Y_SNP, X = X_SNP, k = k[i], family = "gaussian", standardize = F,
                 fast = F)
  SIS_fit <- SIS(x = X_SNP, y = Y_SNP, family = "gaussian", nsis = k[i], iter = F,
               standardize = F)$ix0
  suppressWarnings(Lasso_fit <- glmnet(x = X_SNP,y = Y_SNP, family = "gaussian" , pmax = k[i]))
  
  Lasso_index <- (1 : dim(X_SNP)[2])[Lasso_fit$beta[,dim(Lasso_fit$beta)[2]] != 0]
  
  screened <- c(list(colnames(X_SNP)[SIS_fit]),list(colnames(X_SNP)[SMLE_fit$ID_retained]),
                
                list(colnames(X_SNP)[Lasso_index]))
  #Causal Snps' positions
  snps_postion <- match(unlist(snps_chr),pvals$SNP)
  #Causal Snps' effects
  effects = as.numeric(c(effects1,effects2,effects3,effects4))
  results = matrix(nrow = length(levels), ncol = 3)
  distances = data.frame(snps = snps_postion, effects = effects, results)
  colnames(distances)[3 : 5] = c("SIS","SMLE","glmnet")
  for (l in 1:nrow(distances)){
    for (j in 1:3){
      subset = match(screened[[j]],pvals$SNP)
      distances[l,j+2] = min(abs(distances[l,1] - subset))
    }
  }
  #distances$w_min_SIS = abs(distances$effects) * distances$SIS
  #distances$w_min_SMLE = abs(distances$effects) * distances$SMLE
  md <- colMeans(distances[,3 : 5])
  MD <- rbind(MD,md)
  
  # Section 4.3 ---- figure  ------
  
  if (i == 3) {
    # Set y-axis limits for this specific case
    ylimits <- c(-100, 2000)
  
    distances_long <- data.frame(
      EffectSize = abs(distances[, 2]),
      Distance = c(distances[, 3], distances[, 4], distances[, 5]),
      Method = factor(rep(c("SIS", "SMLE", "glmnet"), each = nrow(distances)))
    )
    
    md <- c(md[1], md[2], md[3])  
    labels <- paste("Averaged MRD =", round(md, 2))

    plots <- lapply(unique(distances_long$Method), function(method) {
 
      data_plot <- subset(distances_long, Method == method)
      label_text <- labels[which(unique(distances_long$Method) == method)]
      ggplot(data_plot, aes(x = EffectSize, y = Distance)) +
        geom_point() +
        labs(x = "Effect size", y = "Distance", title = paste0(method, "()")) +
        ylim(ylimits) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),  # 增大标题
          axis.title = element_text(size = 20),  # 增大轴标签
          axis.text = element_text(size = 16),  # 增大轴刻度
          plot.margin = margin(20, 20, 20, 20)  # 调整边距
        ) +
        annotate(
          "text", x = max(data_plot$EffectSize) * 0.75, 
          y = ylimits[2] * 0.75, label = label_text, size = 6  # 增大文本注释
        )
    })
  }
}

row.names(MD) <- paste0("k=", k)
colnames(MD) <- c("SIS_MRD","SMLE_MRD","GLMNET_MRD")
MD_df <- as.data.frame(MD)

# Add k values as a separate column by extracting from row names
MD_df$Screening_Size <- as.numeric(sub("k=", "", rownames(MD_df)))

# Reorder columns to put Screening_Size first
MD_df <- MD_df[, c("Screening_Size", "SIS_MRD", "SMLE_MRD", "GLMNET_MRD")]
MD_df$Screening_Size <- as.numeric(MD_df$Screening_Size)
MD_df$SIS_MRD <- as.numeric(MD_df$SIS_MRD)
MD_df$SMLE_MRD <- as.numeric(MD_df$SMLE_MRD)
MD_df$GLMNET_MRD <- as.numeric(MD_df$GLMNET_MRD)
MD_long$Method <- factor(MD_long$Method, levels = c("SIS_MRD", "SMLE_MRD", "GLMNET_MRD"), 
                         labels = c("SIS()", "SMLE()", "glmnet()"))

# Plot using ggplot2
ggplot(MD_long, aes(x = Screening_Size, y = MRD, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("SIS()" = "blue", "SMLE()" = "red", "glmnet()" = "orange")) +
  labs(
    x = "Screening size (k)",
    y = "Averaged MRD",
    color = "Method"
  ) +
  theme_minimal(base_size = 20) +  # Increase base size to make text larger
  theme(
    legend.position = c(0.95, 0.95),  # Top right position
    legend.justification = c("right", "top"),  # Align the legend to the top right
    legend.box.just = "right",  # Box alignment to the right
    plot.title = element_text(size = 24, face = "bold"),  # Larger title
    axis.title = element_text(size = 20),  # Larger axis titles
    axis.text = element_text(size = 18),  # Larger axis labels
    legend.text = element_text(size = 18),  # Larger legend text
    legend.title = element_text(size = 20)  # Larger legend title
  )


library(stringr)
library(SMLE)
library(SIS)
library(abess)

Root<- setwd("C:\\Users\\Gamer PC\\Desktop\\SMLE_bess (1)")

source(paste0(Root,"/tool.R"))
source(paste0(Root,"/Simu_algo.R"))
source(paste0(Root,"/Sets_define.R"))
source(paste0(Root,"/Sbess.R"))
source(paste0(Root,"/aSBess.R"))

Metrics <- list("time" ,"model_size","test_error")

Methods <- list(   "lasso_cv")

Num_methods <- length(Methods)

Metadata <- read.csv("C:\\Users\\Gamer PC\\Dropbox\\My PC (DESKTOP-1G2FU5J)\\Downloads\\Final_Metabolites_IOP_IMP_No_NA.csv")

X <- as.matrix(Metadata[,2:1072],ncol= 1071)
Y <- Metadata[,1073]

Num_simu <- 100

for( i in Metrics ){assign(paste0("Table_",i),matrix(0,nrow=Num_methods,ncol = Num_simu))}

family = gaussian()

for(j in 1:Num_simu){
  
  Data <- list()
  
  train_index <- sample(1:8229,floor(0.8*8229),replace  = FALSE)
  
  test_index <- (1:8229)[!(1:8229 %in% train_index)]
  
  test_size <- length(test_index)
  
  Data$X <- X[train_index,]
  
  testX <- X[test_index,]
  
  Data$Y <- Y[train_index]
  
  testY <- Y[test_index]
   
  for(i in 1:length(Methods)){
    
    time1 <- proc.time()
    
    id <- Run_method(Methods[[i]], k = NULL ,family ,Data)
    
    time2 <- proc.time()
    
    data <- data.frame(Y = Data$Y, X = Data$X[,id])
    
    model <- glm(formula = Y~ .,family = family, data = data)
    
    new_model <- data.frame(Y = testY, X = testX[,id])
    
    mu_hat <- predict(model , newdata = new_model, type = "response")
    
    Table_test_error[i,j] <- crossprod(testY- mu_hat)/length(testY)
    
    Table_time[i,j] <- (time2- time1)[3]
    
    Table_model_size[i,j] <- length(id)
    
    }
  
}

Result_name <- str_subset(ls(),"Table_")

Result <- lapply(Result_name,function(i){
  round(rowSums(get(i))/Num_simu,2)
  
})
names(Result) <- str_replace(Result_name,"Table_", "")

df.Result <- data.frame(Result,row.names = Methods)
print(df.Result)
write.table(df.Result, file = "result_select_log", sep = "\t")
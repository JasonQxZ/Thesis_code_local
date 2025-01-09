
# install.packages(c("MASS","stringr"))

library(stringr)
library(SMLE)
library(SIS)
library(abess)

Root<- setwd("~/SMLE_bess/SMLEatmhrg/")

source(paste0(Root,"/tool.R"))
source(paste0(Root,"/Simu_algo.R"))
source(paste0(Root,"/Sets_define.R"))
source(paste0(Root,"/Sbess.R"))

Metrics <- list("time" ,"model_size","test_error")


Metadata <- read.csv("/home/shared/CLSA/vitCmeta/IOP_Metabolites_Data/Final_Metabolites_IOP_IMP.csv")
cleaned_data <- na.omit(Metadata)[,-1]
#cleaned_data <- data.frame(replicate(1092, rnorm(184)))
#colnames(cleaned_data) <- c(paste0("Std_" ,1:1091),"IOP_adjusted")

for( i in Metrics ){assign(paste0("Table_",i),matrix(0,nrow=Num_methods,ncol = Num_simu))}

family = gaussian()

Num_simu <- 10

for(j in 1:Num_simu){
  
  Data <- list()
  
  train_index <- sample(1:184,160,replace  = FALSE)
  
  test_index <- (1:184)[!(1:184 %in% train_index)]
  
  test_size <- length(test_index)
  
  Data$X <- as.matrix(cleaned_data[train_index,-1092],ncol = 1091)
  
  colnames(Data$X) <- NULL
  
  testX <- as.matrix(cleaned_data[test_index,-1092],ncol = 1091)
  
  colnames(testX) <- NULL
  
  Data$Y <- cleaned_data$IOP_adjusted[train_index]
  
  testY <- cleaned_data$IOP_adjusted[test_index]
   
  for(i in 1:length(Methods)){
    
    time1 <- proc.time()
    
    id <- Run_method(Methods[[i]], k ,family ,Data)
    
    time2 <- proc.time()
    
    data <- data.frame(Y = Data$Y, X = Data$X[,id])
    
    model <- glm(formula = Y~.,family = family, data = data)
    
    new_model <- data.frame(Y= testY, X= testX)
    
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
write.table(df.Result, file = "result_log", sep = "\t")


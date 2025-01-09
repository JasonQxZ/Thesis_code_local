# install.packages(c("MASS","stringr"))

library(stringr)
library(SMLE)
library(SIS)
library(abess)

Root <- "~/Library/CloudStorage/OneDrive-UniversityofOttawa/SMLE_bess/SMLE_bess"

source(paste0(Root,"/tool.R"))
source(paste0(Root,"/Simu_algo.R"))
source(paste0(Root,"/Sets_define.R"))
source(paste0(Root,"/Sbess.R"))

Methods <- list( "SIS", "Lasso","abess","SMLE_Lasso","Sbess")

#Methods <- list("SMLE_zero","SMLE_Lasso","Sbess_S3","Sbess_S4")

Metrics <- list("time", "psr" ,"model_size", "ssr")

Num_methods <- length(Methods)

Num_simu <- 100

for( i in Metrics ){assign(paste0("Table_",i),matrix(0,nrow=Num_methods,ncol = Num_simu))}

N = 300

P = 2000

k = 20

family = poisson()

causal_feature_pos <- c(101,103,105,107,109,111)

causal_features_coef = 0.2*c(2,-2,3,-3,-3,3)


for(j in 1:Num_simu){
  
  Data<- Gen_Data(n = N, p = P, family = family$family,  correlation = "AR", rho = 0.8,
                  pos_truecoef = causal_feature_pos, effect_truecoef = causal_features_coef)
  
  for(i in 1:length(Methods)){
    
    time1 <- proc.time()
    
    id <- Run_method(Methods[[i]], k ,family ,Data)
    
    Table_time[i,j] <- (proc.time()-time1)[3]
    
    Table_psr[i,j] <- sum(causal_feature_pos %in% id)/length(causal_feature_pos)
    
    Table_model_size[i,j] <- length(id)
    
    Table_ssr[i,j] <- floor(sum(causal_feature_pos %in% id)/length(causal_feature_pos))
    
  }
  
}

Result_name <- str_subset(ls(),"Table_")

Result <- lapply(Result_name,function(i){
  round(rowSums(get(i))/Num_simu,2)
  
})
names(Result) <- str_replace(Result_name,"Table_", "")

df.Result <- data.frame(Result,row.names = Methods)
write.table(df.Result, file = "model_comparison.txt", sep = "\t")


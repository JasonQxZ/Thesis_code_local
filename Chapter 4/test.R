

set.seed(123)

N = 250

P = 2000

k = 20

family = poisson()

causal_feature_pos <- seq(1,11,2)

causal_features_coef = 0.2*c(2,-2,3,-3,-3,3)

Data<- Gen_Data(n = N, p = P, family = family$family,  correlation = "AR", rho = 0.8,
                pos_truecoef = causal_feature_pos, effect_truecoef = causal_features_coef)

fit_1 <- Sbess(X = Data$X,Y = Data$Y,family = family,k = k, Strategy = 1)

print(fit_1$ID_out)

fit_2 <- Sbess(X = Data$X,Y = Data$Y,family = family,k = k, Strategy = 4)

print(fit_2$ID_out)

fit_smle <- SMLE(X = Data$X,Y = Data$Y,family = family,k = k,fast=T,coef_initial =rep(0,P),intercept = F)

print(fit_smle$ID_retained)

fit_smle_lasso <- SMLE(X = Data$X,Y = Data$Y,family = family,intercept = F,k = k,fast=T)

print(fit_smle_lasso$ID_retained)
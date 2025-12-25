
setwd("C:/Users/taiwo/OneDrive - University of Cincinnati/HDSSPG")
rm(list = ls())


#######

source("LOGITPGPTmass.R")
source("SIMDATA.R")
source("Eval.R")

## Download and install the skinnybasad package from the supplementary materials on the paper,
## "Skinny Gibbs: A Consistent and Scalable Gibbs Sampler for Model Selection"
## link: https://www.tandfonline.com/doi/suppl/10.1080/01621459.2018.1482754?scroll=top

library(skinnybasad)
library(ncvreg)
library(glmnet)
library(MASS)
library(truncnorm)
library(BayesLogit)


## Cases: 1= indept covariates;  2= ar1 correlation; 3= unstructured correlation
## settings;  1= weak signal + only; 2= strong signal + only 
## settings;  3= weak signal +  and -; 4= strong signal + and -
## rho: the parameter for ar1 under case 2
## rho1- correlation between active covariates :
## rho2- correlation between inactive covariates :
## rho3- correlation between active and inactive covariates :
## p = total covariates:  n = total observations
## density  =  % of active covariates


##############


hdsspg.mat   = c()
fix.mat   = c()
hsgpg.mat = c()
hegpg.mat = c()
st.mat    = c()
et.mat    = c()
lasso.mat = c()
mcp.mat   = c()
scad.mat  = c()
nsim=25


               
n = n               #c(100, 200, 300)
p = p               #c(100, 200, 500)
density = density   #c(1, 2,5, 5, 10, 25, 50, 75,100) 

for(i in 1:nsim){
  
  df0 =simdata(n=n,p=p, density=density, case=1,  set=2,  rho=0.5,  rho1=0.5, rho2=0.5, rho3=0.1)
  
  K=10
  
  choicep <-function(x){return(x - K + qnorm(0.90)*sqrt((x/df0$p)*(1- x/df0$p)))}
  cp = uniroot(choicep, c(1,K))$root
  pr = cp/df0$p
  
  hdsspg.model =HDSS_LPG(Y =df0$E, X=df0$X, a=1, b=1, c=2, s = 1, nmc=5000,method = "HDSSPG")
  
  fix.model = FIX_LPG(Y=df0$E, X=df0$X, a0 =5, delta0=0.01, nmc=5000, method= "FIXPG")
  
  hsgpg.model = HSPG_LPG(Y=df0$E, X=df0$X, c=2, s=1, rq = pr, nmc=5000, method= "HSGPG")
  
  hegpg.model = HEGPG_LPG(Y=df0$E, X=df0$X, c=2, s=1, rq = pr, nmc=5000, method= "HEGPG")
  
  model.lasso <- glmnet::cv.glmnet(x= as.matrix(df0$X), y = unlist(df0$E), family = "binomial", 
                                   alpha = 1, standardize = F)
 
   mcp.model   = ncvreg(X= as.matrix(df0$Xunscaled), y = unlist(df0$E), family = "binomial",
                        penalty = "MCP", gamma = 3,  nlambda=50,  warn = F)

   scad.model  = ncvreg(X= as.matrix(df0$Xunscaled), y = unlist(df0$E), family = "binomial",
                        penalty = "SCAD", gamma = 3.7, nlambda=50,  warn = F)
  
  gammaZ0 = rep(0, df0$p); ind = sample(1:df0$p, K);gammaZ0[ind] =1
  n0 = nrow(df0$X)
  st.model = skinnybasad(X= as.matrix(df0$X), E=df0$E, pr=pr, B0=rep(0,df0$p), Z0=gammaZ0,
                         nsplit=10,a0=0.01,b0=1, modif=1, nburn=2000, niter=3000,
                         printitrsep=1000,maxsize=max(K,sqrt(n0)))

  et.model  =skinnybasad(X=as.matrix(df0$X), E=df0$E, pr=pr,B0=rep(0,df0$p), Z0=gammaZ0,
                        nsplit=10,a0=0.01,b0=1,modif=0, nburn=2000, niter=4000,
                        printitrsep=2000, maxsize=max(K,sqrt(n0)))
  
  ###################################################
  p.act = ceiling((df0$density/100)*df0$p)
  l_act = c(rep(1,p.act),rep(0,df0$p-p.act))

  hdsspglvarhat = ifelse(hdsspg.model$lvarZ>0.5, 1, 0)
  fixlvarhat = ifelse(fix.model$lvarZ>0.5, 1, 0)
  hsgpglvarhat = ifelse(hsgpg.model$lvarZ>0.5, 1, 0)
  hegpglvarhat = ifelse(hegpg.model$lvarZ>0.5, 1, 0)
 stlvarhat = ifelse(st.model$marZ>0.5, 1, 0)
 etlvarhat = ifelse(et.model$marZ>0.5, 1, 0)
  lassolvarhat = ifelse(coef(model.lasso, s = "lambda.1se")[-1]!= 0, 1, 0)
  mcplvarhat   =  ifelse(mcp.model$beta[, mcp.model$convex.min][-1]!= 0, 1, 0)
  scadlvarhat  =  ifelse(scad.model$beta[, scad.model$convex.min][-1]!= 0, 1, 0)
  
  #ifelse(hdsspglvarhat ==1, hdsspg.model$betahat, 0)
  
  eval.hdsspg = evaluation(actual = l_act, predicted = hdsspglvarhat, beta_est= ifelse(hdsspglvarhat ==1, hdsspg.model$betahat, 0), 
                        beta =df0$B, X=as.matrix(df0$X), Y = df0$E, method = "logit")
  
  eval.fix = evaluation(actual = l_act, predicted = fixlvarhat, beta_est = ifelse(fixlvarhat ==1, fix.model$betahat, 0), 
                        beta =df0$B,X=as.matrix(df0$X), Y = df0$E,method = "logit")
  
  eval.hsgpg = evaluation(actual = l_act, predicted=hsgpglvarhat, beta_est=ifelse(hsgpglvarhat ==1, hsgpg.model$betahat, 0), 
                          beta =df0$B, X=as.matrix(df0$X), Y = df0$E, method = "logit")
  
  eval.hegpg = evaluation(actual = l_act, predicted=hegpglvarhat, beta_est=ifelse(hegpglvarhat ==1, hegpg.model$betahat, 0),  
                          beta =df0$B,X=as.matrix(df0$X), Y = df0$E,  method = "logit")
  eval.st = evaluation(actual = l_act, predicted=stlvarhat, beta_est=ifelse(stlvarhat ==1, hsgpg.model$betahat, 0),
                          beta =df0$B, X=as.matrix(df0$X), Y = df0$E, method = "logit")

  eval.et = evaluation(actual = l_act, predicted=etlvarhat, beta_est=ifelse(etlvarhat ==1, hegpg.model$betahat, 0),
                          beta =df0$B,X=as.matrix(df0$X), Y = df0$E,  method = "logit")
  
  eval.lasso = evaluation(actual = l_act, predicted=lassolvarhat, beta_est=coef(model.lasso, s = "lambda.1se")[-1], 
                        beta =df0$B, X=as.matrix(df0$X), Y = df0$E, method = "logit")
  
  eval.mcp = evaluation(actual = l_act, predicted=mcplvarhat, beta_est=mcp.model$beta[, mcp.model$convex.min][-1], 
                       beta =df0$B, X=as.matrix(df0$X), Y = df0$E, method = "logit")
  
  eval.scad = evaluation(actual = l_act, predicted=scadlvarhat, beta_est=scad.model$beta[, scad.model$convex.min][-1],  
                         beta =df0$B,X=as.matrix(df0$X), Y = df0$E,  method = "logit")
  
  
  
  
  hdsspg.mat   = rbind(hdsspg.mat, unlist(eval.hdsspg))
  fix.mat   = rbind(fix.mat, unlist(eval.fix))
  hsgpg.mat = rbind(hsgpg.mat, unlist(eval.hsgpg))
  hegpg.mat = rbind(hegpg.mat, unlist(eval.hegpg))
  st.mat   = rbind(st.mat, unlist(eval.st))
  et.mat   = rbind(et.mat, unlist(eval.et))
  lasso.mat   = rbind(lasso.mat, unlist(eval.lasso))
  mcp.mat   = rbind(mcp.mat, unlist(eval.mcp))
  scad.mat = rbind(scad.mat, unlist(eval.scad))
  
}

results = data.frame(
  HDSSPG = apply(hdsspg.mat, 2, mean),
  FIXPG = apply(fix.mat, 2, mean),
  HSGPG = apply(hsgpg.mat, 2, mean),
  HEGPG = apply(hegpg.mat, 2, mean),
  SGT  = apply(st.mat, 2, mean),
  EGT  = apply(et.mat, 2, mean),
  LASSO   = apply(lasso.mat, 2, mean),
  MCP   = apply(mcp.mat, 2, mean),
  SCAD  = apply(scad.mat, 2, mean)                                                                                                    
)


row.names(results) = c("SEN", "SPE", "MCC", "TP", "FP", "TN", "FN", "rMSE", "MSPE", 
                      "SEN_y", "SPE_y", "ACC_y","TP_y", "FN_y", "FP_y", "TN_y")


#results

write.csv(results,  paste0("lasso_and_others_result/p=",df0$p,"_n=",df0$n,"_density=",df0$density,"_case=",df0$case,"_set=",df0$set,"K=",K,".csv"))


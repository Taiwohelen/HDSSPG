library(mvnfast)
library(BayesLogit)
### Auxilliaryfunction for extracting the dimension
dim_func=function(X, Y,paralist, method){
  n <- nrow(X); p <- ncol(X); Yminus05 = Y-0.5;
  
  if(method == "HDSSPG"| method == "HSGPG" | method == "HEGPG"){
    return(list("n"=n,"p"=p , "Yminus05"= Yminus05))
  }
  
  if(method == "FIXPG"){
    delta0 =paralist$delta0;  a0 = paralist$a0;  q = a0/p;
    lambdasq = delta0*(p^(2+delta0))/n; const = (log(q) - log((1-q))) 
    return(list("n"=n,"p"=p , "Yminus05"= Yminus05, "X"=X,"const" = const, 
                "lambdasq" = lambdasq))
  }
}


### Initial values Function 
library(invgamma)
initial_func =function(Y,X,paralist, method){
  mod_int <- glmnet::glmnet(X, Y, alpha = 1,  family = "binomial")
  beta = mod_int$beta[ , which.min(mod_int$lambda)] ## initial value for beta
  w = BayesLogit::rpg.sp(length(Y), 1, X%*%beta) ### initial for PG latent
  #Yminus05 = Y-0.5
  res <- (Y-0.5)/w  - (X%*%beta);
  e <- sum(beta!=0); rq <-  1/100;
  y <- NA # initial for 1/lambda^2
  if(method == "FIXPG") return(list("beta"=beta,"res"=res,"w"=w))
  if(method == "HDSSPG") return(list("beta"=beta,"res"=res,"e"=e,"w"=w,"rq"=rq,"y"=y))
  if(method == "HSGPG" | method == "HEGPG"){
    gammaZ =  rep(0, ncol(X)); ind  = sample(1:ncol(X), ncol(X)/10)
     gammaZ[ind] =1
  return(list("beta"=beta,"w"=w, "y"= rinvgamma(1, shape=2, rate = 1),
                                     "lvar" = gammaZ)) 
     }
}


HDSS_step=function(X,paralist,auxlist,iterpara){
  n=auxlist$n;p=auxlist$p; Yminus05= auxlist$Yminus05
  a=paralist$a;b=paralist$b;c=paralist$c;s=paralist$s;
  beta=iterpara$beta; res=iterpara$res; e=iterpara$e; rq=iterpara$rq; y=iterpara$y;
  lvar = iterpara$lvar
  w=iterpara$w; #Yminus05 =iterpara$Yminus05
  ### UPDATE rq
  rq=rbeta(1,a+e,p-e+b);
  ### UPDATE y
  y <- rgamma(1, shape = c+(e/2), rate = (s + sum(beta*beta*lvar)/2 ));
  
  ### UPDATE beta and lj (sequentially)
  diagXtWX = diag(t(X)%*%diag(w)%*%X)
  lvar= numeric()
  for(j in 1:p){
    C1 <- diagXtWX[j] + y 
    res <- res +  beta[j]*X[,j] 
    C2 <- sum(res*w*X[,j])
    logcj <- log(rq) - log((1-rq)) +log(sqrt(y/C1)) + (C2^2/(2*C1)); 
    if(is.infinite(exp(logcj))){
      beta[j] = rnorm(1,C2/C1,sqrt(1/C1))
      lj =1
    }else{
      lj = rbinom(1,1,exp(logcj)/(1+exp(logcj)) ) ## update lj
      ### update betaj
      if(lj == 0 ) beta[j] = 0
      if(lj ==1 ) beta[j] = rnorm(1,C2/C1,sqrt(1/C1))
    }
    lvar[j] = lj
    res <- res -  beta[j]*X[,j] 
  }
  e <- sum(beta!=0); ### updated no. of active betas
  w <- BayesLogit::rpg(n,1, X%*%beta)### updated pg latent variables
  res = (Yminus05/w) - (X%*%beta)
  return(list("beta"=beta,"res"=res,"e"=e, "w" = w , "rq"=rq,"y"=y, "lvar" = lvar));
}

###List of functions  for method specification
HDSSLOGITPG=list(init_fun=initial_func,aux_fun=dim_func,step_fun=HDSS_step)


Gib_samp<- function(listfunc,Y,X, paralist, MAX_STEPS, 
                                 hard_burnin=2000, nmc=5000, method){
    ### INITIAL-VALUE-FUNC 
  auxlist=listfunc$aux_fun(X, Y,paralist, method); iterpara=listfunc$init_fun(Y,X,paralist ,method);
  ### Resutls storage
  betamat <- c(); q <- c(); lambdasq <- c();  Wmat <- c();lvarmat<-c()
  
  ## First-step burnin and first sample
  for(i in 1:(nmc+hard_burnin)){
    iterpara=listfunc$step_fun(X,paralist,auxlist,iterpara);
    betamat=cbind(betamat,iterpara$beta);
    q=c(q,iterpara$rq);
    lambdasq=c(lambdasq,1/(iterpara$y));
    Wmat = cbind(Wmat,iterpara$w);
    lvarmat = cbind(lvarmat,iterpara$lvar);
  }
  
  lvarmat = lvarmat[,-(1:(hard_burnin))]
  betamat = betamat[,-(1:(hard_burnin))]
  q = q[-(1:(hard_burnin))]
  lambdasq = lambdasq[-(1:(hard_burnin))]
  
  return(list("betamat" = betamat, "q" = q, "lambdasq" = lambdasq,"lvarmat"=lvarmat));
}


## function for running hdsspg method 
HDSS_LPG <- function(Y, X, a=1, b=1, c=2, s = 1, nmc=5000,method){
  ## Checks and Validations
  
  if(length(Y)!=length(X[,1])) stop("The dimension of Y and X should fit.")
  if(length(Y)<15) stop("The sample size should be at least 15 to let initialization work.")
  paralist=list(a=a, b=b, c=c, s=s)
  res=Gib_samp(listfunc=HDSSLOGITPG,Y=Y, X=X, paralist=paralist, nmc=nmc, method = "HDSSPG")
  betahat=apply(res$betamat,1 ,function(x) mean(x[x!=0]))
  lvarZ = apply(res$lvarmat,1 , mean)
  return(list(betahat=betahat, Gibbs_res=res, lvarZ= lvarZ))
}


#############################################
############
#############

#### Step function for logit Fix PG ###########

FIXPG_step<- function(X,paralist,auxlist,iterpara){
  
  n=auxlist$n;  p=auxlist$p; Yminus05 = auxlist$Yminus05; X = auxlist$X
  lambdasq = auxlist$lambdasq; const = auxlist$const; 
  beta=iterpara$beta; res=iterpara$res; w=iterpara$w;
  invlambdasq = 1/lambdasq
  
  ### UPDATING beta and lj (sequentially)
  diagXtWX = diag(t(X)%*%diag(w)%*%X)
  lvar= numeric()
  for(j in 1:p){
    C1 <- diagXtWX[j] + invlambdasq 
    res <- res +  beta[j]*X[,j] 
    C2 <- sum(res*w*X[,j])
    logcj <- const +log(sqrt(invlambdasq/C1)) + (C2^2/(2*C1)); 
    if(is.infinite(exp(logcj))){
      beta[j] = rnorm(1,C2/C1,sqrt(1/C1))
      lj  =1
    }else{
      lj = rbinom(1,1,exp(logcj)/(1+exp(logcj)) ) ## update lj
      ### update betaj
      if(lj == 0 ) beta[j] = 0
      if(lj ==1 ) beta[j] = rnorm(1,C2/C1,sqrt(1/C1))
    }
    res <- res -  beta[j]*X[,j] 
    lvar[j] <- lj
  }
  
  
  w <- BayesLogit::rpg.sp(n,1, X%*%beta)### updated pg latent variables
  res = Yminus05/w - (X%*%beta) #tranformed binary response 
  
  return(list("beta"=beta, "res"=res, "w" = w, "lvar"= lvar))
}

## list of function for method specification
LOGITPG_FIX=list(init_fun=initial_func,aux_fun=dim_func,step_fun=FIXPG_step)


## function for running ind logit pg 
FIX_LPG <- function(Y, X, a0, delta0, nmc=5000, method) {
  
  if(length(Y)!=length(X[,1])) stop("The dimension of Y and X should fit.")
  if(length(Y)<15) stop("The sample size should be at least 15 to let initialization work.")
  paralist=list(a0=a0, delta0 = delta0 )
  res=Gib_samp(listfunc=LOGITPG_FIX,Y=Y, X=X,paralist=paralist,nmc=nmc, method= "FIXPG")
  betahat=apply(res$betamat,1 ,function(x) mean(x[x!=0]))
  lvarZ = apply(res$lvarmat,1 , mean)
  return(list(betahat=betahat, Gibbs_res=res, lvarZ=lvarZ))
}


###################################################

##### Skinny Gibbs 

#################################################

HSPG_STEP<- function(X,paralist,auxlist,iterpara){
  
  Xnames = colnames(X) 
  
  n=auxlist$n;p=auxlist$p; Eminus05= auxlist$Yminus05
  c=paralist$c; s=paralist$s; rq=paralist$rq;
   
  beta=iterpara$beta;  tau1_2=iterpara$y;
  gammaZ = iterpara$lvar; w=iterpara$w; 
  
  tau0_2 = 1/n; inv_tau0_2 <- 1 / tau0_2
    
    W <- diag(w)
    Y <- Eminus05/w
    idx_active <- which(gammaZ == 1)
    idx_inactive <- which(gammaZ == 0)
    n_active <- length(idx_active)
    inv_tau1_2 <- 1 / tau1_2
    const <- (rq*sqrt(tau0_2))/((1 - rq)*sqrt(tau1_2))
    
    ##### Update beta
    beta <- numeric(p)
    if (n_active == 0) {
      beta <- rnorm(p, mean = 0, sd = 1 / sqrt(n - 1 + inv_tau0_2))
    } else {
      X_active <- as.matrix(X[, idx_active])
      XtWX <- (t(X_active)%*% W %*% X_active) + inv_tau1_2*diag(n_active)
      sigma1 <-  solve(XtWX)
      mean1 <- sigma1 %*% crossprod(X_active, W%*%Y)
      beta[idx_active] <- t(unlist(rmvn(1, mu = mean1, sigma = sigma1))) 
      beta[idx_inactive] <- rnorm(p - n_active, mean = 0, sd = 1 / sqrt(n - 1 + inv_tau0_2))
    }
    
    ##### Update gammaZ
    if (n_active > 0) {
      temp_2_12 <- Y - as.matrix(X[,idx_active]) %*% beta[idx_active]
    } else {
      temp_2_12 <- Y - numeric(n)
    }
    
    beta_active = beta[idx_active]
    
    
    for (j in 1:p) {
      X_col <- X[, j]
      temp_0 <- 0.5 * beta[j]^2 * (inv_tau0_2 - inv_tau1_2)
      temp_2_1 <- beta[j]*X[,j]%*%W
      if (j %in% idx_active){
        temp_2_11 <- (Y-as.matrix(X_active[, -which(idx_active==j)], nrow=n)%*%beta_active[-which(idx_active==j)])
        temp_2 <- temp_2_1%*% temp_2_11
      }else{
        temp_2 <- temp_2_1%*%temp_2_12
      }
      temp_3 <- 0.5*(beta[j]^2) * sum(X[,j]*diag(1-w)*X[,j])
      
      log_d_j <- log(const) + temp_0 + temp_2 + temp_3
      
      if(is.infinite(exp(log_d_j))) {
        gammaZ[j] <- 1
      } else { 
        gammaZ[j] <- rbinom(1, 1, exp(log_d_j)/(1 + exp(log_d_j)) )
      }
    }
  ##Update  tau1 square
    tau1_2 = rinvgamma(1,  shape = (sum(gammaZ)/2 +c), rate = (s + sum(gammaZ*beta*beta)/2) )
  ###Update Omega (w)
    if (n_active > 0) {
      w <- rpg.sp(n, 1,  X_active%*%beta_active)
    }
  
  return(list("beta"=beta,  "w" = w , "y"= tau1_2, "lvar" = gammaZ))
  
}


HSGPG_LST=list(init_fun=initial_func,aux_fun=dim_func,step_fun=HSPG_STEP)

HSPG_LPG <- function(Y, X, c, s,rq, nmc=5000,method = "HSGPG") {
  if(length(Y)!=length(X[,1])) stop("The dimension of Y and X should fit.")
  if(length(Y)<15) stop("The sample size should be at least 15 to let initialization work.")
  paralist=list(c=c, s=s, rq =rq)
  res=Gib_samp(listfunc=HSGPG_LST,Y=Y, X=X,paralist=paralist,nmc=nmc, method = "HSGPG")
  betahat=apply(res$betamat,1 ,function(x) mean(x[x!=0]))
  lvarZ = apply(res$lvarmat,1 , mean)
  return(list(betahat=betahat, Gibbs_res=res, lvarZ=lvarZ))
}



###################################################

##### Hierarchical  Exact Gibbs  Gibbs 

#################################################

HEGPG_STEP <- function(X,paralist,auxlist,iterpara){
    
    Xnames = colnames(X) 
    n=auxlist$n;p=auxlist$p; E_minus_0_5= auxlist$Yminus05
    c=paralist$c; s=paralist$s; rq=paralist$rq;
    
    beta=iterpara$beta;  tau1_2=iterpara$y;
    gammaZ = iterpara$lvar; w=iterpara$w; 
  
   tau0_2 = 1/n; inv_tau0_2 <- 1 / tau0_2
   eps_diag <- .Machine$double.eps*diag(p)
  
    W <- diag(w)
    Y <- E_minus_0_5 / w
    inv_tau1_2 <- 1 / tau1_2
    const <- (rq*sqrt(tau0_2))/((1 - rq)*sqrt(tau1_2))
    idx_active <- which(gammaZ == 1)
    # Update beta
    D_z <- diag((gammaZ*inv_tau1_2) + (1 - gammaZ)*inv_tau0_2)
    V <- t(X) %*% W %*% X + D_z+ eps_diag
    sigma1 <- solve(V)
    mean1 <- sigma1 %*% t(X) %*% W %*% Y
    beta <- t(unlist(rmvn(1, mu = mean1, sigma = sigma1))) 
    
    # Update gammaZ ##### update gamma
    for (j in 1:p){
      temp1 <- beta[j]^2*0.5*(inv_tau0_2-inv_tau1_2)
      log_d_j <- log(const)+temp1
      if(is.infinite(exp(log_d_j)) == T){
        gammaZ[j] <-  1
      }else{
        gammaZ[j] <- rbinom(1, 1, exp(log_d_j)/(1+exp(log_d_j)))
      }
    }
    
    ##Update  tau1 square
    tau1_2 =  rinvgamma(1, shape = (sum(gammaZ)/2 +c), rate = (s + sum(gammaZ*beta*beta)/2))
    
    ###Update Omega (w)
    eta <- X %*% beta
    w <- rpg.sp(n, 1,  X%*%beta)
    
    return(list("beta"=beta,  "w" = w , "y"= tau1_2, "lvar" = gammaZ))
    
}



HEGPG_LST=list(init_fun=initial_func,aux_fun=dim_func,step_fun=HEGPG_STEP)

HEGPG_LPG <- function(Y, X, c, s,rq, nmc=5000,method = "HEGPG") {
  if(length(Y)!=length(X[,1])) stop("The dimension of Y and X should fit.")
  if(length(Y)<15) stop("The sample size should be at least 15 to let initialization work.")
  paralist=list(c=c, s=s, rq =rq)
  res=Gib_samp(listfunc=HEGPG_LST,Y=Y, X=X,paralist=paralist,nmc=nmc, method = "HEGPG")
  betahat=apply(res$betamat,1 ,function(x) mean(x[x!=0]))
  lvarZ = apply(res$lvarmat,1 , mean)
  return(list(betahat=betahat, Gibbs_res=res, lvarZ=lvarZ))
}


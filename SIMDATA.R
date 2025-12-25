


library(mvtnorm)

corfunc = function( p,density, rho1, rho2, rho3){

  pact = (density/100)*p 
  
  covr1=(1- rho1)*diag(pact) + array(rho1,c(pact,pact)) 
  covr3=array(rho3,c(pact,p-pact))
  covr2=(1- rho2)*diag(p-pact) + array(rho2,c(p-pact,p-pact))
  covr=rbind(cbind(covr1,covr3),cbind(t(covr3),covr2))
  return(covr)
}


gen_X <- function(n, p, density, case, rho=NULL,  rho1=NULL, rho2=NULL, rho3=NULL){
  ar1func <- function(n, rho) rho^(abs(matrix(1:n-1, nrow=n, ncol=n, byrow=TRUE)-(1:n-1)))
  csfunc <- function(n, rho) matrix(data=rho, nrow=n, ncol=n)+diag(n)*(1-rho)
  if (case==1) sigm <- diag(p)
  if(case == 2 & is.null(rho)==T) stop("under case2 (ar1) correlation structure, you must specify rho")
  if (case==2) sigm <- ar1func(p, rho)
  if(case == 3 &&  (is.null(rho1)==T |  is.null(rho2)==T |  is.null(rho3)==T)  ){
    stop("under case3 (unstructured) correlation structure, you must specify rho1, rho2 and rho3")
  }
  if (case==3) sigm <- corfunc(p=p, density=density, rho1=rho1 , rho2=rho2, rho3=rho3)

  X <- rmvnorm(n=n, mean=rep(0, p), sigma=sigm)
  return(X)
}


### Generating Beta
gen_Beta <- function(p, density, set){
  set.seed(11032025)
  pact = (density/100)*p 
  beta <- rep(0, p)
  if(pact > 0){
    p0 = floor(pact/2)
    #Weak + only 
    if(set==1) beta[1:pact] <- c(runif(n=pact, min=0.5, max=1.5))
    #Strong + only 
    if (set==2) beta[1:pact] <-c(runif(n=pact, min=1.5, max=3))
    #Weak 
    if(set==3) beta[1:pact] <- c(runif(n=pact-p0, min=0.5, max=1.5), runif(n=p0, min=-1.5, max=-0.5))
    #Strong 
    if (set==4) beta[1:pact] <-c(runif(n=pact-p0, min=1.5, max=3), runif(n=p0, min=-3, max=-1.5))
  }
  return(beta)
}

simdata <- function(n,p, density, case,  set,  rho=NULL,  rho1=NULL, rho2=NULL, rho3=NULL){
 
  X0 = gen_X(n=n, p=p, density=density, case=case, rho=rho,  rho1=rho1, rho2=rho2, rho3=rho3)
  
  Bc = gen_Beta(p, density, set)
  ## Logit probabilities
  Y = plogis(X0%*%Bc) 
  ## binary response variable 
  E = rbinom(n, 1, Y) 
  Xunscaled = as.matrix(X0)
  X0= as.matrix(scale(X0)) ### scaling
  X <- as.matrix(X0)
  return(list(X=X0, E=E, B=Bc, Xunscaled = Xunscaled,n =n, p=p, density =density, 
              case = case, set = set, rho=rho,  rho1=rho1, rho2=rho2, rho3=rho3))
}





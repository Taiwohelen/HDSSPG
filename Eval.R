

library(caret)

evaluation <- function(actual, predicted, beta, beta_est, X, Y,  method = c("logit","probit")){
  
  beta_est <- ifelse(predicted == 1, beta_est, 0)
  
  true.idx <- which(actual==1)
  false.idx <- which(actual==0)
  positive.idx <- which(predicted==1)
  negative.idx <- which(predicted==0)
  
  TP <- length(intersect(true.idx, positive.idx))
  FP <- length(intersect(false.idx, positive.idx))
  FN <- length(intersect(true.idx, negative.idx))
  TN <- length(intersect(false.idx, negative.idx))
  
  
  Sensitivity <- TP/(TP+FN)
  if ((TP+FN)==0) Sensitivity <- 1
  
  Specific <- TN/(TN+FP)
  if ((TN+FP)==0) Specific <- 1
  
  MCC.denom <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  
  if(is.na(MCC.denom) == F){
    
  if (MCC.denom==0) MCC.denom <- 1
  MCC <- (TP*TN-FP*FN)/MCC.denom
  if ((TN+FP)==0) MCC <- 1
  } else{MCC <- NA}
  
  if(method == "logit") p_pred =  plogis(X%*%beta_est)
  if(method == "probit") p_pred =  pnorm(X%*%beta_est)
  y_pred = ifelse(p_pred >= 0.5, 1, 0)
  
  rMSE <- sum((beta_est-beta)^2)/sum(beta*beta)
  MSPE <- mean((p_pred-Y)^2)
  
  confusion_matrix <- confusionMatrix(as.factor(c(y_pred)), as.factor(Y))
  
  betaresults <- list("SEN"=Sensitivity, "SPE"=Specific, "MCC"=MCC,  "TP"=TP, "FP"=FP, 
                      "TN"=TN, "FN"=FN, "rMSE"=rMSE, "MSPE" = MSPE,"SEN_Y" = confusion_matrix$byClass[1], 
                      "SPE_Y" = confusion_matrix$byClass[2],  "ACC_Y" = confusion_matrix$overall[1] )
  
  
  
  return(list("Results_Beta"=betaresults, "Conftab" = confusion_matrix$table))
}
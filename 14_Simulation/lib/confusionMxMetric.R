################################################################################
# List of functions to get metrics derived from confustion matrix
# Source: https://en.wikipedia.org/wiki/Sensitivity_and_specificity
### FUNCTION ###################################################################
confMxMetric <- function(CONFMX=COMPIJMX[,c("TP","FP", "TN", "FN"),],
                         metric="TPR"){
  
  TP <- CONFMX[,"TP"]
  FP <- CONFMX[,"FP"]
  TN <- CONFMX[,"TN"]
  FN <- CONFMX[,"FN"]
  rm(CONFMX); gc()
  
  TPR=function(TP, FP, TN, FN){ TP/(TP+FN) }
  TNR=function(TP, FP, TN, FN){ TN/(TN+FP) }
  PPV=function(TP, FP, TN, FN){ TP/(TP+FP) }
  NPV=function(TP, FP, TN, FN){ TN/(TN+FN) }
  FNR=function(TP, FP, TN, FN){ FN/(FN+TP) }
  FPR=function(TP, FP, TN, FN){ FP/(FP+TN) }
  FDR=function(TP, FP, TN, FN){ FP/(FP+TP) }
  FOR=function(TP, FP, TN, FN){ FN/(FN+TN) }
  PT =function(TP, FP, TN, FN){
    tpr <- TPR(TP, FP, TN, FN) 
    tnr <- TNR(TP, FP, TN, FN)
    return( (sqrt(tpr*(-tnr+1))+tnr-1)/(tpr+tnr-1) )
  }
  TS =function(TP, FP, TN, FN){ TP/(TP+FN+FP) }
  ACC=function(TP, FP, TN, FN){
    (TP+TN)/(TP+TN+FP+FN)
  }
  BA =function(TP, FP, TN, FN){
    (TPR(TP, FP, TN, FN) + TNR(TP, FP, TN, FN))/2
  }
  F1 =function(TP, FP, TN, FN){
    (2*TP)/((2*TP)+FP+FN)
  }
  MCC=function(TP, FP, TN, FN){
    (TP*TN-FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) )
  }
  FM =function(TP, FP, TN, FN){
    sqrt( (TP/(TP+FP)) * (TP/(TP+FN)) )
  }
  BM =function(TP, FP, TN, FN){
    TPR(TP, FP, TN, FN) + TNR(TP, FP, TN, FN) - 1
  }
  MK =function(TP, FP, TN, FN){
    PPV(TP, FP, TN, FN) + NPV(TP, FP, TN, FN) - 1
  }
  
  return( eval(parse(text=paste0(metric,"(TP=TP,FP=FP,TN=TN,FN=FN)"))) )
  
}
################################################################################



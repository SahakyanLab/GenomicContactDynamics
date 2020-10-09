################################################################################
# List of functions to get metrics derived from confustion matrix
# Source: https://en.wikipedia.org/wiki/Sensitivity_and_specificity
### FUNCTION ###################################################################
confMxMetric <- function(CONFMX=COMPIJMX[,c("TP","FP", "TN", "FN", 
                                            "SP", "SN", "RP", "RN"),],
                         metric="TPR"){
  
  TP <- CONFMX[,"TP"]
  FP <- CONFMX[,"FP"]
  TN <- CONFMX[,"TN"]
  FN <- CONFMX[,"FN"]
  SP <- CONFMX[,"SP"]
  SN <- CONFMX[,"SN"]
  RP <- CONFMX[,"RP"]
  RN <- CONFMX[,"RN"]
  rm(CONFMX); gc()
  
  SPN=function(TP, FP, TN, FN, SP, SN, RP, RN){ SP/SN }
  RPN=function(TP, FP, TN, FN, SP, SN, RP, RN){ RP/RN }
  TPR=function(TP, FP, TN, FN, SP, SN, RP, RN){ TP/(TP+FN) }
  TNR=function(TP, FP, TN, FN, SP, SN, RP, RN){ TN/(TN+FP) }
  PPV=function(TP, FP, TN, FN, SP, SN, RP, RN){ TP/(TP+FP) }
  NPV=function(TP, FP, TN, FN, SP, SN, RP, RN){ TN/(TN+FN) }
  FNR=function(TP, FP, TN, FN, SP, SN, RP, RN){ FN/(FN+TP) }
  FPR=function(TP, FP, TN, FN, SP, SN, RP, RN){ FP/(FP+TN) }
  FDR=function(TP, FP, TN, FN, SP, SN, RP, RN){ FP/(FP+TP) }
  FOR=function(TP, FP, TN, FN, SP, SN, RP, RN){ FN/(FN+TN) }
  PT =function(TP, FP, TN, FN, SP, SN, RP, RN){
    tpr <- TPR(TP, FP, TN, FN, SP, SN, RP, RN) 
    tnr <- TNR(TP, FP, TN, FN, SP, SN, RP, RN)
    return( (sqrt(tpr*(-tnr+1))+tnr-1)/(tpr+tnr-1) )
  }
  TS =function(TP, FP, TN, FN, SP, SN, RP, RN){ TP/(TP+FN+FP) }
  ACC=function(TP, FP, TN, FN, SP, SN, RP, RN){
    (TP+TN)/(TP+TN+FP+FN)
  }
  BA =function(TP, FP, TN, FN, SP, SN, RP, RN){
    (TPR(TP, FP, TN, FN, SP, SN, RP, RN) + TNR(TP, FP, TN, FN, SP, SN, RP, RN))/2
  }
  F1 =function(TP, FP, TN, FN, SP, SN, RP, RN){
    (2*TP)/((2*TP)+FP+FN)
  }
  MCC=function(TP, FP, TN, FN, SP, SN, RP, RN){
    # as.numeric() to bypass integer limit and avoid integer overflow
    TP <- as.numeric(TP)
    FP <- as.numeric(FP)
    TN <- as.numeric(TN)
    FN <- as.numeric(FP)
    (TP*TN-FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) )
  }
  FM =function(TP, FP, TN, FN, SP, SN, RP, RN){
    sqrt( (TP/(TP+FP)) * (TP/(TP+FN)) )
  }
  BM =function(TP, FP, TN, FN, SP, SN, RP, RN){
    TPR(TP, FP, TN, FN, SP, SN, RP, RN) + TNR(TP, FP, TN, FN, SP, SN, RP, RN) - 1
  }
  MK =function(TP, FP, TN, FN, SP, SN, RP, RN){
    PPV(TP, FP, TN, FN, SP, SN, RP, RN) + NPV(TP, FP, TN, FN, SP, SN, RP, RN) - 1
  }
  
  return( eval(parse(text=paste0(metric,"(TP=TP,FP=FP,TN=TN,FN=FN,SP=SP,SN=SN,RP=RP,RN=RN)"))) )
  
}
################################################################################



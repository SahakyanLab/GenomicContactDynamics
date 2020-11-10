################################################################################
#                                                                              #
#        Function that returns a temperature control unit object definition    #
#        as a string to support OptimusSA and OptimusRE                        #
#                                                                              #
#******************************************************************************#
# NumofAccRatSMIdeal -- Number of times acceptance ratio is sequentially
#                       smaller than ideal.
# NumofAccRatGRIdeal -- Number of times acceptance ratio is sequentially
#                       greater than ideal.
# TSCLnum            -- Cutoff for one of the above numbers after which the
#                       adjustment step for the temperature will be twice
#                       increased or decreased.
#*******************************************************************************
################################################################################
tempControlDefinition <- function(){
  'tempControlUnit <- setRefClass("tempControlUnit",
                                  fields = list(NumofAccRatSMIdeal = "numeric",
                                                NumofAccRatGRIdeal = "numeric",
                                                t.adjstep = "numeric",
                                                T.ADJSTEP = "numeric",
                                                AccR.category = "character",
                                                new.T.INI = "numeric",
                                                instanceOFswitch = "numeric",
                                                minT = "numeric",
                                                max = "numeric",
                                                scaling = "numeric",
                                                DELTA = "numeric"),
                                  methods = list(
                                    updateTemp = function(currentT, idealRatio, currentRatio){
                                      newT <- currentT
                                      if(abs(currentRatio-idealRatio) < DELTA){
                                        return(currentT)
                                      }
                                      if(currentRatio < idealRatio){#-- CATEGORY - SMALLER
                                        NumofAccRatSMIdeal <<- NumofAccRatSMIdeal + 1
                                        NumofAccRatGRIdeal <<- 0
                                        if(NumofAccRatSMIdeal==max){
                                          t.adjstep <<- scaling*t.adjstep; NumofAccRatSMIdeal <<- 0
                                        }
                                        if(AccR.category=="GREATER"){
                                          t.adjstep <<- T.ADJSTEP; instanceOFswitch <<- instanceOFswitch + 1
                                          if(instanceOFswitch==1){ new.T.INI <<- currentT }
                                        } # resets t.adjstep, if it has just crossed the ideality line.
                                        newT <- currentT + t.adjstep
                                        AccR.category <<- "SMALLER"
                                      }
                                      if(currentRatio >= idealRatio){ #-- CATEGORY - GREATER
                                        NumofAccRatSMIdeal <<- 0
                                        NumofAccRatGRIdeal <<- NumofAccRatGRIdeal + 1
                                        if(NumofAccRatGRIdeal==max){
                                          t.adjstep <<- scaling*t.adjstep; NumofAccRatGRIdeal <<- 0
                                        }
                                        if(AccR.category=="SMALLER"){
                                          t.adjstep <<- T.ADJSTEP; instanceOFswitch <<- instanceOFswitch + 1
                                          if(instanceOFswitch==1){ new.T.INI <<- currentT }
                                        } # resets t.adjstep, if it has just crossed the ideality line.
                                        newT <- currentT - t.adjstep
                                        #-- not allowing minus temperatures (absolute scale) and 0 (DE/T=NaN)
                                        if(newT <= 0){ newT <- minT }
                                        AccR.category <<- "GREATER"
                                      }
                                      return(newT)},
                                    resetTemp = function(){
                                      instanceOFswitch <<- 0
                                      return(new.T.INI)
                                    })
  )'
}

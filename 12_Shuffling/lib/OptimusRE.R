################################################################################
#                                                                              #
#        Replica Exchange Acceptance Ratio Monte Carlo optimiser               #
#        Featuring:                                                            #
#        * an automatic adaptive pseudo-temperature bath control,              #
#        * memory friendly processing capability,                              #
#        * progress visualisation,                                             #
#                                                                              #
#*******************************************************************************
#-- NUMITER        = 1000000 # Number of model optimisation steps.
#-- EXCHANGE.FREQ  = 1000    # Frequency of exchanges in steps (NUMITER should
#                              be divisable by this number).
#*******************************************************************************
#-- ACCRATIO   = c(90, 50, 5, 1)  # Vector of Acceptance Ratios for each
#                                   replica (length of ACCRATIO must be equal to
#                                   NCPU).
#-- STATWINDOW = 70               # Number of last ongoing iterations to
#                                   calculate acceptance ratio for temperature
#                                   auto-adjustment.
#*******************************************************************************
#-- T.INI     = 0.00001      # Initial temperature (K) for Metropolis criterion.
#-- T.ADJSTEP = 0.000000005  # Temperature change step-size for temperature
#                              auto-adjustment based on the actual acceptance
#                              ratio.
#-- TSCLnum   = 2            # Cutoff for one of the NumofAccRatSMIdeal and
#                              NumofAccRatGRIdeal numbers after which the
#                              adjustment step is multiplied by T.SCALING.
#-- T.SCALING = 3            # See above.
#-- T.MIN     = 0.000000005  # Value to which the pseudo-temperature is set when
#                              the temperatures control unit makes T negative.
#-- T.DELTA   = 2            # Minimum value by which acceptance ratio in a
#                              STATWINDOW must differ from the ideal acceptance
#                              ratio for the temperature control unit to make a
#                              temperature adjustment.
#*******************************************************************************
#-- DUMP.FREQ = 10000        # The frequency (in steps) of writing the found
#                              ongoing best model.
#*******************************************************************************
#-- LIVEPLOT      = TRUE     # Plotting the optimisation process in a pdf file.
#-- LIVEPLOT.FREQ = 100000   # Frequency (in steps) of plotting the results.
#-- PDFheight     = 29       # Plot height in inches.
#-- PDFwidth      = 20       # Plot width in inches.
#*******************************************************************************
#-- NCPU    = 4              # Number of CPU cores/replicas to use (NCPU must be
#                              greater than 1 and equal to the length of
#                              ACCRATIO).
#-- LONG    = TRUE           # If TRUE, it means that a long simulation is
#                              expected to be done, hence the memory-friendly
#                              mode will be activated.
#-- SEED    = 840            # Setting the seed for the random number generator.
#-- OPTNAME = ""             # The name of the optimisation process.
#*******************************************************************************
#-- K.INITIAL                # The initial parameter configuration from which
#                              the optimisation process will begin.
#-- mDEF                     # Model function that operates on the parameters to
#                              be optimized (K) and returns an observable
#                              object O.
#-- uDEF                     # Function that evaluates the performance of a
#                              given set of parameters K.
#-- rDEF                     # Function that defines a rule by which the
#                              paramters(K) are randomly altered.
#-- DATA = NULL              # A list that holds any supplementary data that
#                              functions mDEF or uDEF need to access.
#*******************************************************************************
################################################################################
#' @import foreach
#' @import iterators
#' @import parallel
#' @importFrom doParallel registerDoParallel
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics lines par plot
#' @importFrom stats rbinom runif

OptimusRE <- function(NUMITER       = 1000000,
                      STATWINDOW    = 70,
                      T.INI         = 0.00001,
                      T.ADJSTEP     = 0.000000005,
                      TSCLnum       = 2,
                      T.SCALING     = 3,
                      T.MIN         = 0.000000005,
                      T.DELTA       = 2,
                      DUMP.FREQ     = 10000,
                      LIVEPLOT      = TRUE,
                      LIVEPLOT.FREQ = 100000,
                      PDFheight     = 29,
                      PDFwidth      = 20,
                      NCPU          = 4,
                      LONG          = TRUE,
                      SEED          = 840,
                      OPTNAME       = "",
                      EXCHANGE.FREQ = 1000,
                      ACCRATIO      = c(90, 50, 5, 1),
                      DATA          = NULL,
                      K.INITIAL     = 0,
                      rDEF,
                      mDEF,
                      uDEF
){
  ################################################################################

  #-- Check that at least 2 processors are available
  if(!(NCPU > 1))
    stop("Replica Exchange requires at least 2 processors, ideally 4 or more.")
  #-- Check that the number of processors matches the number of acceptance ratios
  if(length(ACCRATIO)!= NCPU)
    stop("The number of processors must match the length of the acceptance ratio vector")

  set.seed(SEED)
  suppressWarnings(requireNamespace("foreach"))
  suppressWarnings(requireNamespace("iterators"))
  suppressWarnings(requireNamespace("parallel"))
  suppressWarnings(requireNamespace("doParallel"))
  registerDoParallel(cores = NCPU)
  `%op%` <- `%dopar%`
  print(paste("Running ", NCPU, " replicas.", sep=""), quote=F)

  seeds <- floor(runif(n=NCPU, min=1, max=1000))

  #-- INITIALISATION # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #-- The target energies of all the NUMITER trials.
  ENERGY.TRIAL.VEC.cache      <- list()
  #-- The energy of the system, updated only if the move is accepted.
  ENERGY.DE.FACTO.cache       <- list()
  #-- Probability of acceptance for each step.
  PROB.VEC.cache              <- list()
  #-- Acceptance vector for each step: 1/0 = accept/don't accept.
  ACCEPTANCE.cache            <- list()
  #-- The real temperature for each step.
  T.DE.FACTO.cache            <- list()
  #-- The calculated de facto acceptance rate over the last STATWINDOW.
  ACC.VEC.DE.FACTO.cache      <- list()
  #-- The corresponding checkpoint steps.
  STEP4ACC.VEC.DE.FACTO.cache <- list()
  #-- The quality of the prediction.
  Q.STRG.cache                <- list()
  #-- The stored step.
  STEP.STORED.cache           <- list()
  for(i in 1:NCPU){
    ENERGY.TRIAL.VEC.cache[[i]]      <- 0
    ENERGY.DE.FACTO.cache[[i]]       <- 0
    PROB.VEC.cache[[i]]              <- 0
    ACCEPTANCE.cache[[i]]            <- 0
    T.DE.FACTO.cache[[i]]            <- 0
    ACC.VEC.DE.FACTO.cache[[i]]      <- 0
    STEP4ACC.VEC.DE.FACTO.cache[[i]] <- 0
    Q.STRG.cache[[i]]                <- 0
    STEP.STORED.cache[[i]]           <- 0
  }

  #-- Temperature control unit variables
  NumofAccRatSMIdeal.cache <- rep(0, NCPU)
  NumofAccRatGRIdeal.cache <- rep(0, NCPU)
  t.adjstep.cache          <- rep(T.ADJSTEP, NCPU)
  AccR.category.cache      <- rep("INITIAL", NCPU)
  new.T.INI.cache          <- rep(T.INI, NCPU)
  instanceOFswitch.cache   <- rep(0, NCPU)

  #-- The aimed ideal acceptance ratio.
  IDEAL.ACC.VEC <- ACCRATIO
  #-- initial big value, so that the first step will always be accepted!
  E.old.vec     <- rep(100000, NCPU)
  #-- initial T for starting further adjustments
  Tvec          <- rep(T.INI, NCPU)
  #-- initial meaningless value!
  Q.old.vec     <- rep(0, NCPU)
  #-- initial big value for the lowest energy model stored
  E.stored.vec  <- rep(100000, NCPU)
  K.stored.vec  <- list()
  for(i in 1:NCPU){K.stored.vec[[i]] <- 0}
  O.stored.vec  <- list()
  for(i in 1:NCPU){O.stored.vec[[i]] <- 0}

  Step.stored.cache   <- rep(1, NCPU)
  STEP.add.cache      <- rep(1, NCPU)
  #-- vector to store current configuration across all replicas
  configs <- list()
  for(i in 1:NCPU){configs[[i]] <- K.INITIAL}

  tempControlDefinitionAsString <- tempControlDefinition()
  r <- rDEF
  m <- mDEF
  u <- uDEF
  #^^ INITIALISATION # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  print("Replica Exchange Monte Carlo optimisation", quote=FALSE)

  for(EXCHANGE in 1:EXCHANGE.FREQ){

    #-- PARALLEL PROCESSING WRAP # # # # # # # # #
    suppressWarnings(result <- foreach(repl=1:NCPU, .inorder=FALSE, .export = ls(environment())) %op% {

      eval(parse(text = tempControlDefinitionAsString))

      #initialize a temperature control unit
      tempControl <- tempControlUnit$new(NumofAccRatSMIdeal = NumofAccRatGRIdeal.cache[repl],
                                         NumofAccRatGRIdeal = NumofAccRatSMIdeal.cache[repl],
                                         t.adjstep = t.adjstep.cache[repl],T.ADJSTEP = T.ADJSTEP,
                                         AccR.category = AccR.category.cache[repl],
                                         new.T.INI = new.T.INI.cache[repl],
                                         instanceOFswitch = instanceOFswitch.cache[repl],
                                         minT = T.MIN, max = TSCLnum, scaling = T.SCALING,
                                         DELTA = T.DELTA)
      set.seed(seeds[repl])

      #-- retrieve replica variables from caches
      ENERGY.TRIAL.VEC      <- ENERGY.TRIAL.VEC.cache[[repl]]
      ENERGY.DE.FACTO       <- ENERGY.DE.FACTO.cache[[repl]]
      PROB.VEC              <- PROB.VEC.cache[[repl]]
      ACCEPTANCE            <- ACCEPTANCE.cache[[repl]]
      T.DE.FACTO            <- T.DE.FACTO.cache[[repl]]
      ACC.VEC.DE.FACTO      <- ACC.VEC.DE.FACTO.cache[[repl]]
      STEP4ACC.VEC.DE.FACTO <- STEP4ACC.VEC.DE.FACTO.cache[[repl]]
      Q.STRG                <- Q.STRG.cache[[repl]]
      STEP.STORED           <- STEP.STORED.cache[[repl]]
      Step.stored           <- Step.stored.cache[repl]
      STEP.add              <- STEP.add.cache[repl]

      K        <- configs[[repl]]
      T        <- Tvec[repl]
      E.old    <- E.old.vec[repl]
      Q.old    <- Q.old.vec[repl]
      E.stored <- E.stored.vec[repl]
      K.stored <- K.stored.vec[[repl]]
      O.stored <- O.stored.vec[[repl]]

      #-- Execute MC iterations until next exchange occurs
      for(INT in 1:(NUMITER/EXCHANGE.FREQ)){

        #-- update the STEP parameter (for indexing)
        STEP <- EXCHANGE * (NUMITER/EXCHANGE.FREQ) + INT

        K.new <- r(K=K)
        O     <- m(K=K.new, DATA = DATA)
        snp   <- u(O=O, DATA = DATA)

        Q       <- snp$Q
        E       <- snp$E                      ; ENERGY.TRIAL.VEC[STEP.add]<- E
        DE      <- E - E.old                  ; T.DE.FACTO[STEP.add]      <- T
        #-- Metropolis criterion
        P.accept<- min(1, exp(-DE/T))         ; PROB.VEC[STEP.add]        <- P.accept
        accept  <- rbinom(1, 1, prob=P.accept); ACCEPTANCE[STEP.add]      <- accept

        ################
        if( accept==1 ){
          E.old  <- E
          K      <- K.new
          Q.old  <- Q

          #-- # # # # # # # # #
          if(E.old < E.stored){
            E.stored       <- E.old
            Step.stored    <- STEP
            K.stored       <- K
            O.stored       <- O
            DUMP.MODEL     <- NULL
            DUMP.MODEL[1]  <- "QUALITY:"
            DUMP.MODEL[2]  <- paste("E: ",round(E.old,3),sep="")
            DUMP.MODEL[3]  <- paste("Q: ",round(Q.old,3),sep="")
            DUMP.MODEL[4]  <- paste("Acceptance Ratio: ", IDEAL.ACC.VEC[repl], sep="")
            DUMP.MODEL[5]  <- "TERMS:"
            DUMP.MODEL[6]  <- paste(as.character(names(K.stored)), collapse=" ")
            DUMP.MODEL[7]  <- "COEFFICIENTS:"
            DUMP.MODEL[8]  <- paste(format(as.vector(K.stored), scientific=F, trim=T), collapse=" ")
            #        DUMP.MODEL[9]  <- "OBSERVABLES:"
            #        DUMP.MODEL[10]  <- paste(as.character(names(O.stored)), collapse=" ")
            #        DUMP.MODEL[11] <- "PREDICTIONS:"
            #        DUMP.MODEL[12] <- paste(format(as.vector(O.stored), scientific=F, trim=T), collapse=" ")
            #        DUMP.MODEL[13] <- "TARGET:"
            #        DUMP.MODEL[14] <- paste(format(as.vector(y), scientific=F, trim=T), collapse=" ")
            #        DUMP.MODEL[15] <- "################################"
          }
          #-- # # # # # # # # #

        }
        ################

        ENERGY.DE.FACTO[STEP.add] <- E.old
        Q.STRG[STEP.add]          <- Q.old
        STEP.STORED[STEP.add]     <- STEP

        #-- Adjusting the temperature every STATWINDOW STEP
        if((length(ACCEPTANCE)%%STATWINDOW) == 0){ #-- thus the number of entries is n*STATWINDOW
          #-- determine the new acceptance ratio
          accept.ratio.inlastWIN <- 100*sum(ACCEPTANCE[(length(ACCEPTANCE)-STATWINDOW+1):length(ACCEPTANCE)])/
            STATWINDOW
          ACC.VEC.DE.FACTO <- c(ACC.VEC.DE.FACTO, accept.ratio.inlastWIN)
          STEP4ACC.VEC.DE.FACTO <- c(STEP4ACC.VEC.DE.FACTO, STEP)

          #-- call the temperature control unit
          T <- tempControl$updateTemp(currentT = T, idealRatio = IDEAL.ACC.VEC[repl], currentRatio = accept.ratio.inlastWIN)
        }

        #-- # # # # # # # # # # # # # # # # # # # # # # # #

        #-- PLOTTING  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        if(LIVEPLOT==TRUE){
          if(STEP%%LIVEPLOT.FREQ == 0){

            pdf(width=PDFwidth, height=PDFheight, file=paste(OPTNAME,repl,".pdf",sep=""))
            par(mfrow=c(5,1))

            plot(y=PROB.VEC, x=STEP.STORED, ylab="Acceptance P", xlab="Step",
                 ylim=c(0,1), xlim=range(STEP.STORED), type="l", col="blue", lwd=1)
            lines(y=c(0,1), x=rep(Step.stored,2), lwd="2.5", col="green")

            plot(y=T.DE.FACTO, x=STEP.STORED, ylab="Temperature", xlab="Step",
                 xlim=range(STEP.STORED), type="l", col="orange", lwd=2)
            lines(y=range(T.DE.FACTO), x=rep(Step.stored,2), lwd="2.5", col="green")

            plot(y=rep(IDEAL.ACC.VEC[repl],length(STEP.STORED)), x=STEP.STORED,
                 ylab="Acceptance ratios (%)", xlab="Step", ylim=c(0,100),
                 xlim=range(STEP.STORED), type="l",lty="dashed", lwd=2)
            lines(y=c(0,100), x=rep(Step.stored,2), lwd="2.5", col="green")
            if(length(ACC.VEC.DE.FACTO)!=0){
              lines(y=ACC.VEC.DE.FACTO, x=STEP4ACC.VEC.DE.FACTO, ylim=c(0,100),
                    xlim=range(STEP.STORED), col="red", lwd=2)
            }

            plot(y=ENERGY.DE.FACTO, x=STEP.STORED, ylab="Energy", xlab="Step",
                 xlim=range(STEP.STORED), type="l", col="navy", lwd=2)
            lines(y=range(ENERGY.DE.FACTO), x=rep(Step.stored,2), lwd="2.5",
                  col="green")

            plot(y=Q.STRG, x=STEP.STORED, ylab="Q", xlab="Step",
                 xlim=range(STEP.STORED), type="l", col="blue", lwd=2)
            lines(y=range(Q.STRG), x=rep(Step.stored,2), lwd="2.5", col="green")

            dev.off()

          }
        }
        #^^ PLOTTING  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


        if(STEP%%DUMP.FREQ == 0 & exists("DUMP.MODEL") ){
          write(DUMP.MODEL, file=paste(OPTNAME,repl,"_model_QE.log",sep=""))
          save(K.stored,    file=paste(OPTNAME,repl,"_model_K.Rdata",sep=""))
          save(O.stored,    file=paste(OPTNAME,repl,"_model_O.Rdata",sep=""))
        }

        STEP.add <- STEP.add + 1    ########

        #-- Trimming the data holding vectors to be memory-friendly if LONG==TRUE
        if(LONG==TRUE){
          if(STEP%%10000 == 0 & STEP!=10000){
            ENERGY.DE.FACTO       <- ENERGY.DE.FACTO[(length(ENERGY.DE.FACTO)-10000+1):length(ENERGY.DE.FACTO)]
            Q.STRG                <- Q.STRG[(length(Q.STRG)-10000+1):length(Q.STRG)]
            STEP.STORED           <- STEP.STORED[(length(STEP.STORED)-10000+1):length(STEP.STORED)]
            ENERGY.TRIAL.VEC      <- ENERGY.TRIAL.VEC[(length(ENERGY.TRIAL.VEC)-10000+1):length(ENERGY.TRIAL.VEC)]
            T.DE.FACTO            <- T.DE.FACTO[(length(T.DE.FACTO)-10000+1):length(T.DE.FACTO)]
            PROB.VEC              <- PROB.VEC[(length(PROB.VEC)-10000+1):length(PROB.VEC)]
            ACCEPTANCE            <- ACCEPTANCE[(length(ACCEPTANCE)-10000+1):length(ACCEPTANCE)]
            ACC.VEC.DE.FACTO      <- ACC.VEC.DE.FACTO[(length(ACC.VEC.DE.FACTO)-(10000/STATWINDOW)+1):length(ACC.VEC.DE.FACTO)]
            STEP4ACC.VEC.DE.FACTO <- STEP4ACC.VEC.DE.FACTO[(length(STEP4ACC.VEC.DE.FACTO)-(10000/STATWINDOW)+1):length(STEP4ACC.VEC.DE.FACTO)]
            STEP.add              <- 10001
          }
        }
        #^^ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

      } ##-- for(EXCHANGE in 1:EXCHANGE.FREQ)
      #^^ Execute MC iterations until next exchange occurs  # # # # # # # # # # #

      #-- if LONG==TRUE, OUTPUT will only hold the trimmed data!
      OUTPUT                       <- NULL

      OUTPUT$E.stored              <- E.stored
      OUTPUT$E.old                 <- E.old
      OUTPUT$Q.old                 <- Q.old
      OUTPUT$K                     <- K
      OUTPUT$T                     <- T
      OUTPUT$K.stored              <- K.stored
      OUTPUT$O.stored              <- O.stored

      OUTPUT$STEP                  <- STEP
      OUTPUT$PROB.VEC              <- PROB.VEC
      OUTPUT$T.DE.FACTO            <- T.DE.FACTO
      OUTPUT$IDEAL.ACC.VEC         <- IDEAL.ACC.VEC
      OUTPUT$ACC.VEC.DE.FACTO      <- ACC.VEC.DE.FACTO
      OUTPUT$STEP4ACC.VEC.DE.FACTO <- STEP4ACC.VEC.DE.FACTO
      OUTPUT$ENERGY.DE.FACTO       <- ENERGY.DE.FACTO
      OUTPUT$Q.STRG                <- Q.STRG
      OUTPUT$ACCEPTANCE            <- ACCEPTANCE
      OUTPUT$STEP.STORED           <- STEP.STORED

      OUTPUT$Step.stored           <- Step.stored
      OUTPUT$ENERGY.TRIAL.VEC      <- ENERGY.TRIAL.VEC
      OUTPUT$STEP.add              <- STEP.add

      OUTPUT$NumofAccRatSMIdeal    <- tempControl$NumofAccRatSMIdeal
      OUTPUT$NumofAccRatGRIdeal    <- tempControl$NumofAccRatGRIdeal
      OUTPUT$t.adjstep             <- tempControl$t.adjstep
      OUTPUT$AccR.category         <- tempControl$AccR.category
      OUTPUT$new.T.INI             <- tempControl$new.T.INI
      OUTPUT$instanceOFswitch      <- tempControl$instanceOFswitch

      save(OUTPUT, file=paste(OPTNAME,repl,"_model_ALL.Rdata", sep=""))

      OUTPUT
    })
    #^^ PARALLEL PROCESSING WRAP # # # # # # # # #

    #-- Store the output from parallel workers in the respective caches
    for(i in 1:NCPU){
      data <- result[[i]]

      E.old.vec[i] <- data$E.old
      Q.old.vec[i] <- data$Q.old
      configs[[i]] <- data$K
      Tvec[i]      <- data$T

      E.stored.vec[i]   <- data$E.stored
      K.stored.vec[[i]] <- data$K.stored
      O.stored.vec[[i]] <- data$O.stored

      ENERGY.TRIAL.VEC.cache[[i]]      <- data$ENERGY.TRIAL.VEC
      ENERGY.DE.FACTO.cache[[i]]       <- data$ENERGY.DE.FACTO
      PROB.VEC.cache[[i]]              <- data$PROB.VEC
      ACCEPTANCE.cache[[i]]            <- data$ACCEPTANCE
      T.DE.FACTO.cache[[i]]            <- data$T.DE.FACTO
      ACC.VEC.DE.FACTO.cache[[i]]      <- data$ACC.VEC.DE.FACTO
      STEP4ACC.VEC.DE.FACTO.cache[[i]] <- data$STEP4ACC.VEC.DE.FACTO
      Q.STRG.cache[[i]]                <- data$Q.STRG
      STEP.STORED.cache[[i]]           <- data$STEP.STORED
      STEP.add.cache[i]                <- data$STEP.add
      Step.stored.cache[i]             <- data$Step.stored

      NumofAccRatSMIdeal.cache[i] <- data$NumofAccRatSMIdeal
      NumofAccRatGRIdeal.cache[i] <- data$NumofAccRatGRIdeal
      t.adjstep.cache[i]          <- data$t.adjstep
      AccR.category.cache[i]      <- data$AccR.category
      new.T.INI.cache[i]          <- data$new.T.INI
      instanceOFswitch.cache[i]   <- data$instanceOFswitch
    }

    ######-- Choose two adjacent replicas at random and swap them
    index1 <- sample(1:NCPU,1); index2 <- NULL
    if(index1 == 1)
      index2 <- 2
    else{
      if(index1 == NCPU)
        index2 <- NCPU - 1
      else
      {
        x <- sample(c(-1, 1),1)
        index2 <- index1 + x
      }
    }

    temp <- NULL

    temp$E.old    <- E.old.vec[index1]
    temp$Q.old    <- Q.old.vec[index1]
    temp$K        <- configs[[index1]]
    temp$E.stored <- E.stored.vec[index1]
    temp$K.stored <- K.stored.vec[[index1]]

    E.old.vec[index1]      <- E.old.vec[index2]
    Q.old.vec[index1]      <- Q.old.vec[index2]
    configs[[index1]]      <- configs[[index2]]
    E.stored.vec[index1]   <- E.stored.vec[index2]
    K.stored.vec[[index1]] <- K.stored.vec[[index2]]

    E.old.vec[index2]      <- temp$E.old
    Q.old.vec[index2]      <- temp$Q.old
    configs[[index2]]      <- temp$K
    E.stored.vec[index2]   <- temp$E.stored
    K.stored.vec[[index2]] <- temp$K.stored
    #####

    #-- Reset the corresponding temp control unit values
    NumofAccRatSMIdeal.cache[index1] <- 0
    NumofAccRatGRIdeal.cache[index1] <- 0
    t.adjstep.cache[index1]          <- T.ADJSTEP
    AccR.category.cache[index1]      <- "INITIAL"
    new.T.INI.cache[index1]          <- T.INI
    instanceOFswitch.cache[index1]   <- 0

    NumofAccRatSMIdeal.cache[index2] <- 0
    NumofAccRatGRIdeal.cache[index2] <- 0
    t.adjstep.cache[index2]          <- T.ADJSTEP
    AccR.category.cache[index2]      <- "INITIAL"
    new.T.INI.cache[index2]          <- T.INI
    instanceOFswitch.cache[index2]   <- 0
    #####
  }
  print("An optimal model is obtained and the data are saved !!!", quote=FALSE)
  if(NCPU==1){ return(OUTPUT) } else { return(0) }

}

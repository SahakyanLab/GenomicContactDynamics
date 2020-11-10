################################################################################
#                                                                              #
#        Interface supporting Monte Carlo Optimisation with Simulated          #
#        Annealing and Replica Exchange Acceptance Ratio Monte Carlo           #
#        Optimisation                                                          #
#                                                                              #
################################################################################

#' Interface for Acceptance Ratio Annealing Monte Carlo Optimisation
#' and Replica Exchange Acceptance Ratio Monte Carlo Optimisation.
#'
#' @param NUMITER       Number of model optimisation steps.
#' @param STATWINDOW    Number of last ongoing iterations to calculate
#'                      acceptance ratio for temperature auto-adjustment.
#' @param T.INI         Initial temperature (K) for Metropolis criterion.
#' @param T.ADJSTEP     Temperature change step-size for temperature
#'                      auto-adjustment based on the actual acceptance ratio.
#' @param TSCLnum       Cutoff for one of the NumofAccRatSMIdeal and
#'                      NumofAccRatGRIdeal numbers after which the adjustment
#'                      step is multiplied by T.SCALING.
#' @param T.SCALING     See above.
#' @param T.MIN         Value to which the pseudo-temperature is set when the
#                       temperatures control unit makes T negative.
#' @param T.DELTA       Minimum value by which acceptance ratio in a STATWINDOW
#'                      must differ from the ideal acceptance ratio for the
#'                      temperature control unit to make a temperature
#'                      adjustment.
#' @param DUMP.FREQ     The frequency (in steps) of writing the found ongoing
#'                      best model.
#' @param LIVEPLOT      Plotting the optimisation process in a pdf file.
#' @param LIVEPLOT.FREQ Frequency (in steps) of plotting the results.
#' @param PDFheight     Plot height in inches.
#' @param PDFwidth      Plot width in inches.
#' @param NCPU          Number of CPU cores to use, by running more replicas
#'                      of the optimisation. The usage of more than 1 cores
#'                      will attempt to load the foreach and doParallel
#'                      libraries in the case of SA Optimus
#' @param LONG          If TRUE, it means that a long simulation is expected to
#'                      be done, hence the memory-friendly mode will be
#'                      activated.
#' @param SEED          Setting the seed for the random number generator.
#' @param OPTNAME       The name of the optimisation process.
#' @param DATA          A list that holds any supplementary data that
#'                      functions mDEF or uDEF need to access.
#' @param K.INITIAL     The initial parameter configuration from which the
#'                      optimisation process will begin.
#' @param rDEF          Function that defines a rule by which the
#'                      paramters(K) are randomly altered.
#' @param mDEF          Model function that operates on the parameters to be
#'                      optimized (K) and returns an observable object O.
#' @param uDEF          Function that evaluates the performance of a given set
#'                      of parameters K.
#' @param EXCHANGE.FREQ Frequency of exchanges in steps (NUMITER should be
#'                      divisable by this number, for RE Optimus).
#' @param ACCRATIO      Vector of Acceptance Ratios for each replica (length of
#'                      ACCRATIO must be equal to NCPU, for RE Optimus).
#' @param CYCLES        Number of annealing cycles (NUMITER should be divisable
#'                      by this number, for SA Optimus).
#' @param ACCRATIO.IN   Initial acceptance ratio (%) at the beginning of each
#'                      constituent annealing cycles (for SA Optimus).
#' @param ACCRATIO.FIN  Final acceptance ratio (%) at the end of each
#'                      constituent annealing cycles (for SA Optimus).
#' @param OPT.TYPE      String specifying which optimisation protocol to use.
#'                      Enter "SA" for Simulated Annealing or "RE" for Replica
#'                      Exchange (default value is "SA")
#' @return A probabilistic optimal parameter configuration K.
#'
#' @export

Optimus <- function(NUMITER       = 1000000,
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
                    DATA          = NULL,
                    K.INITIAL     = 0,
                    rDEF,
                    mDEF,
                    uDEF,
                    EXCHANGE.FREQ = 1000,
                    ACCRATIO      = c(90, 50, 5, 1),
                    CYCLES        = 10,
                    ACCRATIO.IN   = 90,
                    ACCRATIO.FIN  = 0.5,
                    OPT.TYPE = "SA"
                    ){
  #-- Call Simulated Annealing Optimus
  if(OPT.TYPE == "SA"){
    OptimusSA(NUMITER       = NUMITER,
              STATWINDOW    = STATWINDOW,
              T.INI         = T.INI,
              T.ADJSTEP     = T.ADJSTEP,
              TSCLnum       = TSCLnum,
              T.SCALING     = T.SCALING,
              T.MIN         = T.MIN,
              T.DELTA       = T.DELTA,
              DUMP.FREQ     = DUMP.FREQ,
              LIVEPLOT      = LIVEPLOT,
              LIVEPLOT.FREQ = LIVEPLOT.FREQ,
              PDFheight     = PDFheight,
              PDFwidth      = PDFwidth,
              NCPU          = NCPU,
              LONG          = LONG,
              SEED          = SEED,
              OPTNAME       = OPTNAME,
              CYCLES        = CYCLES,
              ACCRATIO.IN   = ACCRATIO.IN,
              ACCRATIO.FIN  = ACCRATIO.FIN,
              DATA          = DATA,
              K.INITIAL     = K.INITIAL,
              rDEF          = rDEF,
              mDEF          = mDEF,
              uDEF          = uDEF
              )
  } else {
    #-- Call Replica Exchange Optimus
    if(OPT.TYPE == "RE")
      OptimusRE(NUMITER       = NUMITER,
                STATWINDOW    = STATWINDOW,
                T.INI         = T.INI,
                T.ADJSTEP     = T.ADJSTEP,
                TSCLnum       = TSCLnum,
                T.SCALING     = T.SCALING,
                T.MIN         = T.MIN,
                T.DELTA       = T.DELTA,
                DUMP.FREQ     = DUMP.FREQ,
                LIVEPLOT      = LIVEPLOT,
                LIVEPLOT.FREQ = LIVEPLOT.FREQ,
                PDFheight     = PDFheight,
                PDFwidth      = PDFwidth,
                NCPU          = NCPU,
                LONG          = LONG,
                SEED          = SEED,
                OPTNAME       = OPTNAME,
                EXCHANGE.FREQ = EXCHANGE.FREQ,
                ACCRATIO      = ACCRATIO,
                DATA          = DATA,
                K.INITIAL     = K.INITIAL,
                rDEF          = rDEF,
                mDEF          = mDEF,
                uDEF          = uDEF
                )
    else
      #-- If OPT.TYPE is neither "SA" nor "RE"
      print(paste(OPT.TYPE," is not a recognized optimisation protocol.", sep = ""))
  }
}

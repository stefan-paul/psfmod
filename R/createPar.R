##' Function that mixes parameter with defaults
##'
##' @param init named list with initial values; alternatively all initial values can be passed in a vector (in the right order!).
##' @param pars named list with parameter values; alternatively all parameters can be passed in a vector (in the right order!).
##' @return Complete set of parameters that can be used in \code{\link{psfmod}}
##' @details The init list can have the following elements: \cr
##'  \itemize{
##'  \item{"BA"}{ Biomass species A}
##'  \item{"BB"}{ Biomass species A}
##'  \item{"S"}{ Nutrient availability in absence of plants}
##'  \item{"DA"}{ Litter biomass species A}
##'  \item{"DB"}{ Litter biomass species B}
##'  \item{"L0"}{ Light supply}
##'  \item{"Sa"}{ Proportional abundance of soil community attributed to species A}
##'  \item{"Dca"}{ Proportional effciency of decomposer community compared to monoculture species A}
##'  \item{"Dcb"}{ Proportional effciency of decomposer community compared to monoculture species A}
##'   }
##'  The pars list contains the model parameters. Parameters marked with AB are a vector with two entries,
##'  whereas the first one corresponds to species A and the second to species B.
##'  The list can have the following elements: \cr
##'    \itemize{
##'  \item{"gLC"}{ Maximum growth rate AB under light limitation}
##'  \item{"kLC"}{ Light availability at which half the maximal growth rates of AB are reached}
##'  \item{"gNC"}{ Maximum growth rate AB under nutrient limitation}
##'  \item{"kNC"}{ Nitrogen availability at which half the maximal growth rates of AB are reached}
##'  \item{"mC"}{ Mortality rate AB}
##'  \item{"a"}{ Turnover rate of nutrient supply}
##'  \item{"qNC"}{ Nitrogen content of tissue AB}
##'  \item{"rho"}{ Soil bulk density}
##'  \item{"l_Root"}{ Rooting depth}
##'  \item{"QNC"}{ Nitrogen content of AB litter at which it decomposes at rate dC}
##'  \item{"alphaNC"}{ Light-litter feedback coefficient AB}
##'  \item{"dC"}{ Litter decomposition rate AB}
##'  \item{"cRate"}{ Parameter determining the characteristic timescale of light dynamics}
##'  \item{"gamma_LC"}{ Light interception coefficient AB}
##'  \item{"alpha_LC"}{ Light-litter feedback coefficient AB}
##'  \item{"alpha_a"}{ Influence of Sa on A}
##'  \item{"alpha_b"}{ Influence of Sa on B}
##'  \item{"beta_a"}{ Influence of Sb on A}
##'  \item{"beta_b"}{ Influence of Sb on B}
##'  \item{"sccr"}{Vector with two vaues indicating the speed with which the soil cummunities react to changes
##'  in the conditions. The first value represents the decomposers, the second the soil competitors. A value of 1 means
##'  immediate reaction, whereas values < 1 indicate delayed reaction.}
##'  \item{"hfa"}{Determines whether the efiiciency of the soil composers is based on litter quality (value 0) or
##'  "home field advantage". All values != 0 determine the strength of the home field advantage.}
##'   }
##'
##' @export

createPar <- function(init = NULL, pars = NULL){

  yini <-  list("BA" = 1, "BB" = 1, "S" = 60, "DA" = 0, "DB" = 0, "L0" = 55,
                "Sa" = 0.5, "Dca" = 1, "Dcb" = 1)

  if(is.null(init)) init <- yini
  if(is.numeric(init)){
    if(length(init) != length(yini)) stop("Length of 'init' does not match number of inital values")
    for(i in 1:length(init)) yini[i] <- init[i]
    init <- yini
  }

  y <- modifyList(yini, init)

  S <- y$S
  L0 <- y$L0

  defaults <- list("gLC" = c(0.25, 0.25),
            "kLC" = c(50, 21),
            "gNC" = c(0.25, 0.25),
            "kNC" = c(30, 35),
            "mC" = c(0.005, 0.01),
            "a" = 0.005,
            "S" = S,
            "qNC" = c(15, 15),
            "rho" = 530,
            "l_Root" = 1,
            "QNC" = c(15, 15),
            "alphaNC" = c(0.7, 0.7),
            "dC" = c(0.003, 0.003),
            "cRate" = 10,
            "L0" = L0,
            "gamma_LC" = c(0.03, 0.04),
            "alpha_LC" = c(0.01, 0.013),
            "alpha_a" = 0,
            "alpha_b" = 0,
            "beta_a" = 0,
            "beta_b" = 0,
            "sccr" = c(1,1),
            "hfa" = 0
  )


   if(is.null(pars)) pars <- defaults
   if(is.numeric(pars)){

     if(length(pars) != length(unlist(defaults))) stop("Length of 'pars' does not match
                                               number of parameters")
     for(i in 1:length(defaults)){
       for(k in 1:length(defaults[[i]])) defaults[[i]][k] <- pars[i]

     }
       pars <- defaults
   }

    parNew <- modifyList(defaults, pars)

    # check for consitency

    pars <- list(yini = unlist(y), par = unlist(parNew))
    # gNinit
    if((pars$par[30] + pars$par[32]) < (- pars$par[1])) warning("Feedback parameter could lead to negative growth rates, gN. Will be corrected to 0.")
    else if(pars$par[30] < (- pars$par[1]) | pars$par[32] < (- pars$par[1])  ) warning("Feedback parameter could lead to negative growth rates, gN. Will be corrected to 0.")
    else if((pars$par[31] + pars$par[33]) < (-pars$par[2])) warning("Feedback parameter could lead to negative growth rates, gN. Will be corrected to 0.")
    else if(pars$par[31] < -(pars$par[2]) | pars$par[33] < -(pars$par[2])  ) warning("Feedback parameter could lead to negative growth rates, gN. Will be corrected to 0.")

    # gLinit
    if((pars$par[30] + pars$par[32]) < (-pars$par[5])) warning("Feedback parameter could lead to negative growth rates, gL. Will be corrected to 0.")
    else if(pars$par[30] < (-pars$par[5]) | pars$par[32] < (-pars$par[5])  ) warning("Feedback parameter could lead to negative growth rates, gL. Will be corrected to 0.")
    else if((pars$par[31] + pars$par[33]) < (-pars$par[2])) warning("Feedback parameter could lead to negative growth rates, gL. Will be corrected to 0.")
    else if(pars$par[31] < (-pars$par[6]) | pars$par[33] < (-pars$par[6] ) ) warning("Feedback parameter could lead to negative growth rates, gL. Will be corrected to 0.")

   class(pars) <- "psfpars"

  return(pars)


}

##' Function to update par
##'
##' @param out output of simulation
##' @param par parameter object as created with \code{\link{createPar}}
##' @param newVal optional named list of new parameter values. This can e.g. be used to simulatie mowing by removing biomass.
##' @export

updatePar <- function(out, par, newVal =  NULL){

  # Get size of out
  iter <- nrow(out)

  # Get current values
  yini <- list("BA" = out[iter,2], "BB" = out[iter,3],"S" =  out[iter,4],"DA" =  out[iter,5],
               "DB" = out[iter,6],"L0"= out[iter,7] ,
               "Sa"= out[iter,8],"Dca"= out[iter,9],"Dcb" = out[iter,10])

  if(is.null(newVal)) newVal <- list()
  yini <- modifyList(yini, newVal)

  # Update parameters
  par$yini <- unlist(yini)
  return(par)
}

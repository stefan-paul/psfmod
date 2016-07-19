##' Convenient wrapper function to call the model
##'
##' @param days number of days to run the modle
##' @param pars list with parameter values as created by \code{\link{createPar}}
##' @param out output of previous simulation. If provided function will restart the simulation
##' @param newVal named list with new parameter values. This can e.g. be used to simulatie mowing by removing biomass.
##' @return matrix with the results. Can be diplayed using \code{\link{psfplot}}
##' @export
##' @references Eppinga, Maarten B., et al. "Litter feedbacks, evolutionary change and exotic plant invasion." Journal of Ecology 99.2 (2011): 503-514. \cr\cr Bever, James D. "Soil community feedback and the coexistence of competitors: conceptual frameworks and empirical tests." New Phytologist 157.3 (2003): 465-473.
##' @examples
##' # First definition of a paramater set. I choose to only modify som
##' # of the starting values
##' pars <- createPar(init = list("S" = 20, "L0" = 60))
##'
##' # Now the simulation can be run
##' out <- psfmod(pars = pars, days = 5000)
##'
##' # And the results plotted
##' psfplot(out)

psfmod <- function(pars, days, newVal = NULL, out = NULL){
  times <- 1:days

  if(!is.null(out)) restart <- TRUE
  else restart <- FALSE

  if(restart) pars <- updatePar(out = out, par = pars, newVal = newVal)

  y <- pars[[1]]
  pars <- pars[[2]]
  res <- deSolve::ode(y = y, func = feedbackModel,
           times = times, parms = pars)

  if(restart) {
    res <- rbind(out, res)
    res[,1] <- 1:nrow(res)
  }
  class(res) <- "psfmod"
  return(res)
}


##' Main feedback model
##'
##' @param y initial values (see Details)
##' @param t t
##' @param mu parameters
##' @details For further explanation of the parameter or the model see \code{\link{psfmod}} or the package vignette.
##' @return matrix with results.
##'
feedbackModel <- function (t, y, mu) {

  # Extract parameter objects from parameter vector
  gLinit <- mu[1:2]
  kL <- mu[3:4]
  gNinit <- mu[5:6]
  kN <- mu[7:8]
  m <- mu[9:10]
  a <- mu[11]
  S <- mu[12]
  qN <- mu[13:14]
  rho <- mu[15]
  l_Root <- mu[16]
  QN <- mu[17:18]
  alphaN <- mu[19:20]
  d <- mu[21:22]
  cRate <- mu[23]
  L0 <- mu[24]
  gamma_L <- mu[25:26]
  alpha_L <- mu[27:28]
  alpha <- mu[29:30] # see Bever 2003
  beta <- mu[31:32] # see Bever 2003
  sccr <- mu[33:34] # soil community change rate for decomposers and competitors
  qmax <- QN # plant tissue N content which has maximum decomposition rate
  hfa <- mu[35]

  gL <- gN <- numeric(2)



  # Extract state variables from state variable vector
  B <- y[1:2] # Biomass: Carex and Phalaris
  N <- y[3] # Nutrients
  D <- y[4:5] # Litter mass:  Carex and Phalaris
  L <- y[6] # Light
  S_a <- y[7] #S_a in Bever 2003
  dc <- y[8:9]

  gL[1] <- max(0,(gLinit[1] + alpha[1]*S_a + beta[1]*(1-S_a)))
  gL[2] <- max(0,(gLinit[2] + alpha[2]*S_a + beta[2]*(1-S_a)))

  gN[1] <- max(0,(gNinit[1] + alpha[1]*S_a + beta[1]*(1-S_a)))
  gN[2] <- max(0,(gNinit[2] + alpha[2]*S_a + beta[2]*(1-S_a)))


  # Check for Hme Field Advantage
  if(hfa == 0) dc_dot <- - (dc- ((qN[1]/qmax[1])*(D[1]/(sum(D))) + (qN[2]/qmax[2])*(D[2]/(sum(D)))))*sccr[1]
  else dc_dot <- c(- (dc[1]- (B[1]/sum(B)))*sccr[3], - (dc[2]- (B[2]/sum(B)))*sccr[1])

  Growth <- pmin(gL * L /(kL + L), gN * N / (kN + N)) * B
  B_dot <-  Growth - m * B

  N_dot <- a * (S - N) - sum(qN / (rho * l_Root) * Growth) +
    sum(qN[1]^2 / (QN[1] * rho * l_Root) * alphaN[1] * d[1] * D[1] * dc[1], qN[2]^2 / (QN[2] * rho * l_Root) * alphaN[2] * d[2] * D[2] * dc[2])

  D_dot <- c(m[1] * B[1] - qN[1]/QN[1] * d[1] * D[1] * dc[1],
             m[2] * B[2] - qN[2]/QN[2] * d[2] * D[2] * dc[2])


  L_dot <- cRate * (L0 - sum(gamma_L * B * L) - sum(alpha_L * D * L) - L)

  # S_dot <- S_a*(1-S_a)*(B[1]/(B[1]+B[2]) - v*(B[2]/(B[1]+B[2])))*sccr  # From Bever
  # S_dot <- S_a*(1-S_a)* ((B[1] - v*B[2])/sum(B)  )   # From Revilla


  S_dot <- - (S_a- (B[1]/sum(B)))*sccr[2] # simple version by me



  list(c(B_dot, N_dot, D_dot, L_dot, S_dot, dc_dot))
}


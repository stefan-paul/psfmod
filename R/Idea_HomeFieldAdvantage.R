
feedbackModel2 <- function (t, y, mu) {

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
  sccr <- mu[33:35] # soil community change rate for decomposers and competitors
  qmax <- QN # plant tissue N content which has maximum decomposition rate
  gL <- gN <- numeric(2)

  hfa <- mu[36]

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


S = L0 = 10

# Parameters
refpars <- list("gLC" = c(0.030, 0.028),  # Maximum growth rate under light limitation #DONE
                "kLC" = c(20, 20),  #      # Light availability under which 1/2 gLC is reached
                "gNC" = c(0.028, 0.030), # Maximum growth rate under nutrient limitation #DONE
                "kNC" = c(20, 20), #   # N availability under which 1/2 gNC is reached
                "mC" = c(0.001, 0.001), #   # per capita mortality rate
                "a" = 0.05, #   # turnover rate of nutrient supply. Can later be determined
                # because it determiens the biomass in equilibrium (see Eppinga 2011 Appendix)
                "S" = S, #DONE
                "qNC" = c(18, 18), #  DONE # N content of tissue
                "rho" = 1200, # DONE # soil bulk density
                "l_Root" = 1, #DONE # rooting depyth
                "QNC" = c(18, 18), #  DONE
                "alphaNC" = c(0.7, 0.7), #DONE # Nutrient litter feedback
                "dC" = c(0.015, 0.015), #DONE # Litter decomposition rate
                "cRate" = 10, #DONE
                "L0" = L0, #
                "gamma_LC" = c(0.001, 0.001), #  # light interception coefficient
                "alpha_LC" = c(0.001, 0.001), #  # light litter feedback coefficient
                "alpha_a" = 0, # ?
                "alpha_b" = 0, # ?
                "beta_a" = 0, # partially DONE
                "beta_b" = -0.02, #DONE
                "sccr" = c(1,1,1), #DONE
                "hfa" = 1
)



yini <- list("BA" = 100, "BB" = 10,"S" =  S,"DA" =  10, "DB" = 100,"L0"=  L0,
             "Sa"= 1,"Dca"= 1,"Dcb" = 1)



out <- deSolve::ode(y = unlist(yini), times = 0:5000, parms = unlist(refpars), func = feedbackModel2)


#' Visualization of temporal dynamics
#' ================================
mycol <- rainbow(6)
matplot(out[, 1], out[,c(-1,-8,-9,-10)], type = "l", col = mycol, lty = 1,
        xlab = "Time (years)", ylab = "Values")


out[1:10,]
plot(out[,9], type = "l")  # TODO home field advantage is not working yet....


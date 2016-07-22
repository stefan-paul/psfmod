S <- 60
L0 <- 55
yini <-  c(1, 1, S, 0, 0, L0,0.5,1,1)


pars <- c("gLC" = 0.25, 0.25,
          "kLC" = 50, 21,
          "gNC" = 0.25, 0.25,
          "kNC" = 30, 35,
          "mC" = 0.005, 0.01,
          "a" = 0.005,
          "S" = S,
          "qNC" = 15, 15,
          "rho" = 530,
          "l_Root" = 1,
          "QNC" = 15, 15,
          "alphaNC" = 0.7, 0.7,
          "dC" = 0.003, 0.003,
          "cRate" = 10,
          "L0" = L0,
          "gamma_LC" = 0.03, 0.04,
          "alpha_LC" = 0.01, 0.013,
          "v" = 1,
          "alpha_a" = 0,
          "alpha_b" = -0.2,
          "beta_a" = 0,
          "beta_b" = 0,
          "sccr" = 1, 1,
          "hfa" = 0
)

indices <- list("gLC" = 1:2, "kLC" = 3:4, "gNC" = 5:6, "kNC" = 7:8,
               "mC" = 9:10, "a" = 11, "S" = 12, "qNC" = 13:14, "rho" = 15,
               "l_Root" = 16, "QNC" = 17:18, "alphaNC" = 19:20, "dC" = 21:22,
               "cRate" = 23, "L0" = 24, "gamma_LC" = 25:26, "alpha_LC" = 27:28,
               "alpha_a" = 29, "alpha_b" = 30, "beta_a" = 31,
               "beta_b" = 32, "sccr" = 33:34, "hfa" = 35)

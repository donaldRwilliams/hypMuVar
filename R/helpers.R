# helpers
hypSD_helper <- function(a0_spike, a1_spike,
                         a0_slab, a1_slab, pi){

  # holder for random effects SD
  a0 <- a1 <- matrix(0, 1, 2)

  # holder for rho sigma
  rho_sd <- matrix(0, 1, 3)

  # holder for truncation
  T1 <- T2 <- matrix(0, 1, 3)

  # covariance matrix
  Sigma <- matrix(0, 2, 2)

  # spike and slab components (RE sd)
  ## shape
  a0[1] <- a0_spike
  a0[2] <- a0_slab
  ## rate
  a1[1] <- a1_spike
  a1[2] <- a1_slab

  pi <- rbinom(n = 1, size = 1, prob = pi)
  piI <- pi + 1

  re_sd <- rgamma(1, shape =  a0[piI], rate = a1[piI])

  data.frame(Indicator = piI, RE_sd = re_sd)


}

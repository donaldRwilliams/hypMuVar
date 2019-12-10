re_int_both <- "model{
  for(i in 1:N){
y[i] ~ dnorm(yhat[i], 1/exp(shat[i])^2)
yhat[i] <- inprod(X_loc[i,], beta) + u[ID[i],1]
shat[i] <- inprod(X_scl[i,], eta)  + u[ID[i], 2]
}

for(j in 1:J){
u[j, 1:2] ~ dmnorm.vcov(bhat[j,], Sigma[1:2,1:2])
bhat[j,1] <- 0
bhat[j,2] <- 0
}

Sigma[1:2,1:2] <- Tau %*% Omega %*% Tau

Omega[1,1] <- 1
Omega[2,2] <- 1
Omega[1,2] <- rho[1]
Omega[2,1] <- Omega[1,2]
"

re_slp_loc <- "
for(j in 1:J){
u[j, 1:3] ~ dmnorm.vcov(bhat[j,], Sigma[1:3,1:3])
bhat[j,1] <- 0
bhat[j,2] <- 0
bhat[j,3] <- 0
}

Sigma[1:3,1:3] <- Tau %*% Omega %*% Tau

Omega[1,1] <- 1
Omega[2,2] <- 1
Omega[3,3] <- 1

# location int and slp (rho_01)
Omega[1,2] <- rho[1]
Omega[2,1] <- Omega[1,2]

# location int and scale int (rho_02)
Omega[1,3] <- rho[2]
Omega[3,1] <- Omega[1,3]

# location slp and scale int (rho_12)
Omega[2,3] <- rho[3]
Omega[3,2] <- Omega[2,3]
"



re_slp_scl <- "
for(j in 1:J){
u[j, 1:3] ~ dmnorm.vcov(bhat[j,], Sigma[1:3,1:3])
bhat[j,1] <- 0
bhat[j,2] <- 0
bhat[j,3] <- 0
}

Sigma[1:3,1:3] <- Tau %*% Omega %*% Tau

Omega[1,1] <- 1
Omega[2,2] <- 1
Omega[3,3] <- 1

# location int and scale int (rho_01)
Omega[1,2] <- rho[1]
Omega[2,1] <- Omega[1,2]

# location int and scale slope (rho_02)
Omega[1,3] <- rho[2]
Omega[3,1] <- Omega[1,3]

# scale slp and scale int (rho_12)
Omega[2,3] <- rho[3]
Omega[3,2] <- Omega[2,3]
"




Omega_re_slp_both <- "
for(j in 1:J){
u[j, 1:4] ~ dmnorm.vcov(re_mu[], Sigma[1:4,1:4])
# bhat[j,1] <- 0
# bhat[j,2] <- 0
# bhat[j,3] <- 0
}
re_mu[1] <- 0
re_mu[2] <- 0
re_mu[3] <- 0
re_mu[4] <- 0

Sigma[1:4,1:4] <- Tau %*% Omega %*% Tau

Omega[1,1] <- 1
Omega[2,2] <- 1
Omega[3,3] <- 1
Omega[4,4] <- 1

# location int and location slp (rho_01)
Omega[1,2] <- rho[1]
Omega[2,1] <- Omega[1,2]

# location int and scale int (rho_02)
Omega[1,3] <- rho[2]
Omega[3,1] <- Omega[1,3]

# location int and scale slope (rho_03)
Omega[1,4] <- rho[3]
Omega[4,1] <- Omega[1,4]

# location slope and scale int (rho_12)
Omega[2,3] <- rho[4]
Omega[3,2] <- Omega[2,3]

# location slope and scale slope (rho_13)
Omega[2,4] <- rho[5]
Omega[4,2] <- Omega[2,4]

# scale int and scale slope (rho_23)
Omega[3,4] <- rho[6]
Omega[4,3] <- Omega[3,4]
"

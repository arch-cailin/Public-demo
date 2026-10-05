s# ==================================
# Load libraries
# ==================================
library(nimble)
library(nimbleMacros)
library(parallel)

#=================================
# Load the data
# ================================
load("df_model.Rdata")
source("parallelMCMC.r") # wrapper function written by the Nimble dev team


# ====================================
# Nimble code
# ====================================
oce_m3 <- nimbleCode({
  
  # ==============================
  # Priors
  # ==============================
  
  # --- Mu and Sigma for the MVN y (layer 1) --- 
  
  # LKJ prior for the cholesky factor (correlation structure)
  # (as an upper triangular matrix)
  Ustar[1:p,1:p] ~ dlkj_corr_cholesky(1.3, p)
  
  # Priors for the variances
  for(j in 1:p){
    sds[j] ~dunif(0,10) # non-informative
  }
  
  # --- Beta 0 ---
  beta0[K] <- 0 # for softmax identifiability (setting K as reference class)
  # --- Beta1 (autumn) ---
  beta_aut[K] <- 0
  # --- Beta2 (winter) ---
  beta_win[K] <- 0
  # --- Beta3 (spring) ---
  beta_spr[K] <- 0
  
  for(k in 1:(K-1)){
    beta0[k] ~ dnorm(0, sd = 5) 
    beta_aut[k]  ~ dnorm(0, sd = 5) 
    beta_win[k]  ~ dnorm(0, sd = 5) 
    beta_spr[k]  ~ dnorm(0, sd = 5) 
  }
  
  # --- Spatial effects --- 
  # For softmax identifiability:
  u[1:L, K] <- 0  
  v[1:L, K] <- 0
  
  # Priors for the structured and unstructured spatial variance:
  sigma.u ~ dunif(0,5)
  tau.u <- 1/(sigma.u^2)
  
  sigma.v ~ dunif(0, 5)
  tau.v <- 1 / (sigma.v^2)
  
  
  for(k in 1:(K-1)){
    u[1:L, k] ~ dcar_normal(adj[1:Nadj],
                            weights[1:Nadj],
                            num[1:L],
                            tau.u,
                            zero_mean = 1) # ensure identifiability between
    # beta0 and the u's
    
    for(l in 1:L){
      v[l, k] ~ dnorm(0, tau.v)
    }
  }
  
  
  # ====================================
  # Deterministic quantities
  # ====================================
  
  # --- Scaled cholesky, construction of Sigma--- 
  
  # Scales the cholesky factor (Ustar) by sds (Ustar %*% diag(sds)) 
  # uppertri_mult_diag is a nimble function supplied by the nimbleMacros library
  # for this purpose
  U[1:p,1:p] <- uppertri_mult_diag(Ustar[1:p, 1:p], sds[1:p])
  
  # --- Proportions (softmax) and likelihood ---
  for(i in 1:N){
    for(k in 1:K){
      eta[i, k] <- beta0[k] +
        beta_aut[k]  * aut[i] +
        beta_win[k]  * win[i] +
        beta_spr[k]  * spr[i] +
        u[location[i], k] + 
        v[location[i], k] 
    }
    maxeta[i] <- max(eta[i, 1:K])
    denom[i]  <- sum(exp(eta[i, 1:K] - maxeta[i]))
    for(k in 1:K){
      x[i, k] <- exp(eta[i, k] - maxeta[i]) / denom[i]
    }
    
    # --- Observed means (sums of proportions * source water type means) ---
    
    for(j in 1:p){
      mu.obs[i, j] <- inprod(x[i, 1:K], mu[1:K, j])
    }
    
    y[i, 1:p] ~ dmnorm(mean = mu.obs[i, 1:p],
                       cholesky = U[1:p, 1:p], prec_param = 0)
  }
})
# =====================================
# Transform the data
# =======================================

# Data in matrix form (observations) 
y=as.matrix(df_model[,c("p.temp","salinity")])


# ===================================
# Source water types (means)
# ===================================

# Source water types (means):

Gtable <- read.table(text = "
id type            temp   salinity
1  ICW_27.0        9.00   34.69
2  siAAIW_27.25    6.80   34.57
3  siAAIW_27.40    4.50   34.36
4  siAAIW_27.55    3.51   34.40
5  dAAIW_27.25     4.40   34.17
6  dAAIW_27.40     3.46   34.20
7  dAAIW_27.55     2.78   34.30
8  RSIW_27.25     12.01   35.75
9  RSIW_27.55      9.09   35.48
10 IIW_27.25       6.89   34.61
11 IIW_27.55       4.65   34.60
12 UCDW_Atl_27.80  2.53   34.55
13 UCDW_Ind_27.80  2.67   34.58
14 NIDW_27.80      5.10   35.03", header = TRUE)


# Function to create 3 water types from this table:

avg <- function(rows, name) data.frame(
  water_type = name,
  averaged   = paste(Gtable$type[rows], collapse = " + "),
  temp       = mean(Gtable$temp[rows]),
  salinity   = mean(Gtable$salinity[rows]))

# From Claude: (trying to find the best polygon that will 
# make sense from an oceanography perspective and still 
# contain many of the points)

# 10 of 14 rows, 98.20% of points inside
G <- rbind(
  avg(c(8, 9),          "RSIW"),
  avg(c(12, 13),        "UCDW"),   
  avg(c(2, 3, 4, 5, 6,7), "AAIW"))  

mu=as.matrix(G[,-(1:2)])


# ===========================================
# Dimensions 
# ===========================================

K=3 # number of water classes
N=nrow(y) 
nlocations <- N
p=2 # number of covariates/variables

# =================================================
# 2.5 degree grid, rook adjacency on occupied cells
# =================================================

# Region boundaries
minlat <- -40; maxlat <- -15
minlong <- 25; maxlong <- 50
res <- 2.5                             # grid resolution (degrees)

nx <- (maxlong - minlong) / res        # 10 columns (longitude)
ny <- (maxlat  - minlat)  / res        # 10 rows (latitude)

# Assign each profile to a grid column/row 
lon_bin <- floor((df_model$longitude - minlong) / res) + 1
lat_bin <- floor((df_model$latitude  - minlat)  / res) + 1

# avoid points exactly on the upper boundary
lon_bin <- pmin(lon_bin, nx)
lat_bin <- pmin(lat_bin, ny)
stopifnot(all(lon_bin >= 1), all(lat_bin >= 1))

# Full-grid cell id (row by row from the south-west corner)
cell <- (lat_bin - 1) * nx + lon_bin

# Keep only cells that contain data, renumbered 1..L
occ <- sort(unique(cell))
L   <- length(occ)
df_model$location <- match(cell, occ)

# Column/row of each occupied cell
occ_lon <- (occ - 1) %% nx + 1
occ_lat <- (occ - 1) %/% nx + 1

# Rook neighbours: occupied cells sharing an edge
nb <- lapply(seq_len(L), function(l) {
  d <- abs(occ_lon - occ_lon[l]) + abs(occ_lat - occ_lat[l])
  as.integer(which(d == 1))
})

# Checks: one connected component
stopifnot(all(lengths(nb) > 0))

# Vectors for dcar_normal
num     <- lengths(nb)
adj     <- unlist(nb)
weights <- rep(1, length(adj))
Nadj    <- length(adj)


# ============================================
# Lists of constants, data, initial values for nimble:
# ============================================

constants_m3<-list(N = N,
                   K = K,
                   p = p,
                   L = L,                  
                   adj = adj,              
                   num = num,
                   weights = weights,
                   Nadj = Nadj,
                   dist_main = as.numeric(scale(df_model$main_dist)),
                   aut = df_model$Autumn,
                   win = df_model$Winter,
                   spr = df_model$Spring,
                   location = df_model$location
)


data_list_m3 <- list(y = y, mu = mu)


inits_m3 <- lapply(1:4, function(i) {
  list(
    Ustar = diag(p),
    sds = runif(p, 0.1, 5),
    sigma.u = runif(1, 0.1, 5),
    beta0 = c(rnorm(K-1, 0, 0.5), 0),
    beta_main = c(rnorm(K-1, 0, 0.5), 0),
    beta_aut = c(rnorm(K-1, 0, 0.5), 0),
    beta_win = c(rnorm(K-1, 0, 0.5), 0),
    beta_spr = c(rnorm(K-1, 0, 0.5), 0),
    u = matrix(0, nrow = L, ncol = K),
    v = matrix(0, nrow = L, ncol = K),
    sigma.v = runif(1, 0.1, 5)
  )
})


# ============================================
# Run NIMBLE MCMC
# ============================================

# --- for debugging: do this before running the rest: --- 
m_m3 <- nimbleModel(
  oce_m3,
  constants=constants_m3,
  data=data_list_m3,
  inits=inits_m3[[1]]
)
m_m3$calculate()


mcmc.output_m3 <- nimbleParallelMCMC(
  code = oce_m3, 
  constants = constants_m3,
  data = data_list_m3, 
  inits = inits_m3,
  nchains = 4, 
  niter = 50000, 
  nburnin = 25000,
  summary = TRUE, 
  WAIC = TRUE,
  ncores = 4,
  monitors = c(
    'beta0', 
    'beta_aut', 'beta_win', 'beta_spr',
    'sigma.u', 'sds', 'Ustar', 'u', 'v', 'sigma.v')
)


# ============================================
# Extract Results
# ============================================

# Save complete MCMC output
saveRDS(
  mcmc.output_m3,
  "m3_3waters_grid_unstructured.rds"
)
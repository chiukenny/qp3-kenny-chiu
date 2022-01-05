# Simulation parameters
interf = c(4, 6, 8) # {4,6,8} mean
#interf = 0.5 # sum

sims = 500


# Initialization
# ---------------------------------------------------------

# Load packages
library(rjson)
library(tidyverse)
library(Matrix)

# Load functions from separate file
source("functions.R")

# Set seed for reproducibility
set.seed(1)



# Data loading and cleaning
# ---------------------------------------------------------

# Load nodes and features as data frame
dat_json = fromJSON(file="Data/musae_ENGB_features.json")
dat_mat = matrix(0, 7126, 3)
for (user_str in names(dat_json))
{
  # Adjust IDs 1-based rather than 0-based
  id = as.numeric(user_str) + 1
  dat_mat[id,1] = id
  
  # Extract features 224 and 569
  for (g in unlist(dat_json[user_str]))
  {
    if (g == 224) {dat_mat[id,2] = 1}
    else if (g == 569) {dat_mat[id,3] = 1}
  }
}
dat = as.data.frame(dat_mat)
colnames(dat) = c("id", "game1", "game2")

# Load edge data as sparse matrix
edge_csv = read.csv("Data/musae_ENGB_edges.csv", header=T) + 1
dat_edMat = sparseMatrix(i=edge_csv$from, j=edge_csv$to, symmetric=T)

# TODO
# Add neighbourhood features
# Mean
dat$Ngame1 = as.vector(dat_edMat %*% dat$game1) / rowSums(dat_edMat)
dat$Ngame2 = as.vector(dat_edMat %*% dat$game2) / rowSums(dat_edMat)
# Sum
#dat$Ngame1 = as.vector(dat_edMat %*% dat$game1)
#dat$Ngame2 = as.vector(dat_edMat %*% dat$game2)
dat$N = rowSums(dat_edMat)

# Clean up large files
remove(dat_json)
remove(dat_mat)
remove(edge_csv)



# Simulation: Table 1 - Scenario 1
# ---------------------------------------------------------

res_tab1 = data.frame(beta=integer(),
                      bias_none=double(),
                      bias_ind=double(),
                      bias_all=double())

sim = dat
# TODO
# Simulate individual and neighbourhood treatments Z and G
sim$Z = rbernoulli(nrow(sim), ind_prop(sim$game1, sim$game2))*1
sim$G = as.vector(dat_edMat %*% sim$Z) / rowSums(dat_edMat) # Mean
#sim$G = as.vector(dat_edMat %*% sim$Z) # Sum

for (b in interf)
{
  # Bias when not adjusting for covariates
  bias_t2b = bias_T2B(sim, b)
  
  # Bias when adjusting for individual covariates
  bias_c2i = bias_C2I(sim, b)
  
  # Bias when adjusting for individual and neighbourhood covariates
  bias_c2n = bias_C2N(sim, b)
  
  res_tab1 = rbind(res_tab1, c(b,bias_t2b,bias_c2i,bias_c2n))
}
colnames(res_tab1) = c("beta","bias_none","bias_ind","bias_all")



# Simulation: Table 2 - Scenario 1
# ---------------------------------------------------------

naive_err = matrix(0, sims, length(interf))
reg_err = matrix(0, sims, length(interf))
for (n in 1:sims)
{
  # Use same network across methods in each simulation
  sim = dat
  # TODO
  # Simulate individual and neighbourhood treatments Z and G
  sim$Z = rbernoulli(nrow(sim), ind_prop(sim$game1, sim$game2))*1
  sim$G = as.vector(dat_edMat %*% sim$Z) / rowSums(dat_edMat) # Mean
  #sim$G = as.vector(dat_edMat %*% sim$Z) # Sum
  
  for (j in 1:length(interf))
  {
    b = interf[j]
    
    # Simulate outcomes
    sim$Y = rnorm(nrow(sim), outcome_mean(sim$Z,sim$G,sim$game1,sim$game2,b))
    
    # Compute expected treatment effect for each node
    taus = outcome_mean(1,sim$G,sim$game1,sim$game2,b) -
      outcome_mean(0,sim$G,sim$game1,sim$game2,b)
    
    # Compute error of naive estimator
    naive_est = mean(sim$Y[which(sim$Z==1)]) - mean(sim$Y[which(sim$Z==0)])
    naive_err[n,j] = naive_est - mean(taus)
    
    # Compute error of regression estimator
    reg_err[n,j] = lm(Y~Z+game1+game2, data=sim)$coefficients["Z"] - mean(taus)
  }
}

# Compute bias and RMSE
bias_naive = colMeans(naive_err)
rmse_naive = sqrt(colMeans(naive_err^2))
bias_reg = colMeans(reg_err)
rmse_reg = sqrt(colMeans(reg_err^2))

res_tab2 = data.frame(beta=interf,
                      bias_naive=bias_naive,
                      rmse_naive=rmse_naive,
                      bias_reg=bias_reg,
                      rmse_reg=rmse_reg)
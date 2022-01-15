# Simulation parameters
sum_exposure = T

if (sum_exposure)
{
  # Sum
  spillover = 0.8
} else {
  # Mean
  spillover = 8
}

sims = 25
subclasses = 2:10


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

# Add neighbourhood features
dat$Ngame1 = as.vector(dat_edMat %*% dat$game1) / rowSums(dat_edMat)
dat$Ngame2 = as.vector(dat_edMat %*% dat$game2) / rowSums(dat_edMat)
dat$N = rowSums(dat_edMat)

# Clean up large files
remove(dat_json)
remove(dat_mat)
remove(edge_csv)


# Simulation
# ---------------------------------------------------------

sub_ind_err = matrix(0, sims, length(subclasses))
sub_all_err = matrix(0, sims, length(subclasses))
sub_gps_err = matrix(0, sims, length(subclasses))

system.time(
for (n in 1:sims)
{
  # Use same network across methods in each simulation
  sim = dat
  
  # Simulate individual and neighbourhood treatments Z and G
  sim$Z = rbernoulli(nrow(sim), ind_prop(sim$game1, sim$game2))*1
  if (sum_exposure)
  {
    # Sum
    sim$G = as.vector(dat_edMat %*% sim$Z)
  } else {
    # Mean
    sim$G = as.vector(dat_edMat %*% sim$Z) / rowSums(dat_edMat)
  }
  
  # Simulate outcomes
  sim$Y = rnorm(nrow(sim), outcome_mean(sim$Z,sim$G,sim$game1,sim$game2,b))
  
  # Compute expected treatment effect for each node
  tau = mean(outcome_mean(1,sim$G,sim$game1,sim$game2,spillover) -
               outcome_mean(0,sim$G,sim$game1,sim$game2,spillover))
  
  for (i in 1:length(subclasses))
  {
    j = subclasses[i]
    
    # Compute error of subclass estimators
    sub_ind_err[n,i] = subcl_est(sim, subcl_ind(sim,j)) - tau
    sub_all_err[n,i] = subcl_est(sim, subcl_all(sim,j)) - tau
    
    # Compute error of Forastiere (2021) proposed estimator
    sub_gps_err[n,i] = subcl_gps(sim,j,sum_exposure=sum_exposure) - tau
  }
}
)

# Compute bias and RMSE
res_tab = data.frame(subclass=subclasses,
                     estimator=rep(c("Ind","All","GPS"),each=length(subclasses)),
                     bias=c(colMeans(sub_ind_err),colMeans(sub_all_err),colMeans(sub_gps_err)),
                     rmse=c(sqrt(colMeans(sub_ind_err^2)),sqrt(colMeans(sub_all_err^2)),sqrt(colMeans(sub_gps_err^2))))

res_tab %>%
  ggplot(aes()) +
  geom_line(aes(x=subclass, y=rmse, color=estimator), size=1) +
  labs(x="Subclasses", y="RMSE", color="Estimator") +
  theme_bw()
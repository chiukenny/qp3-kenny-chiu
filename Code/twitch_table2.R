sims = 500

#interf = 4 # {4,6,8} mean
interf = 0.5 # sum


library(rjson)
library(tidyverse)
library(Matrix)


set.seed(1)


dat_json = fromJSON(file = "musae_ENGB_features.json")
dat_mat = matrix(0, 7126, 3)

for (user_str in names(dat_json))
{
  id = as.numeric(user_str) + 1
  dat_mat[id,1] = id
  #feat_i = unname(unlist(dat_json[user]))
  for (g in unlist(dat_json[user_str]))
  {
    if (g == 224) {dat_mat[id,2] = 1}
    else if (g == 569) {dat_mat[id,3] = 1}
  }
  #dat = bind_rows(dat, c(id=id, game1=game1, game2=game2))
  #feats = c(feats, feat_i)
  #users = c(users, rep(user_i, length(feat_i)))
}
dat = as.data.frame(dat_mat)
colnames(dat) = c("id", "game1", "game2")

remove(dat_json)
remove(dat_mat)

edge_csv = read.csv("musae_ENGB_edges.csv", header=T) + 1
dat_edMat = sparseMatrix(i=edge_csv$from, j=edge_csv$to, symmetric=T)
remove(edge_csv)
#degrees = as.data.frame(table(c(dat_csv$from,dat_csv$to))) %>%
#  arrange(desc(Freq))
#summary(degrees$Freq)



# Functions
expit = function(x) {1 / (1+exp(-x))}
ind_prop = function(x1, x2)
{
  return(expit(-3 + 3*x1 + 5*x2))
}
outcome_mean = function(z,g,x1,x2,beta)
{
  prop = ind_prop(x1,x2)
  return(5 + 6*(prop>=0.7) + 10*z - 3*z*(prop>=0.7) + beta*g)
}


# Simulation

naive_err = rep(0, sims)
reg_err = rep(0, sims)

for (n in 1:sims)
{
  sim = dat
  sim$Ngame1 = as.vector(dat_edMat %*% sim$game1)
  sim$Ngame2 = as.vector(dat_edMat %*% sim$game2)
  sim$N = rowSums(dat_edMat)
  
  # Scenario 1
  prop = ind_prop(sim$game1, sim$game2)
  sim$Z = rbernoulli(nrow(sim), prop)*1
  
  #sim$G = as.vector(dat_edMat %*% sim$Z) / rowSums(dat_edMat) # Mean
  sim$G = as.vector(dat_edMat %*% sim$Z) # Sum
  
  mu = outcome_mean(sim$Z,sim$G,sim$game1,sim$game2,interf)
  sim$Y = rnorm(nrow(sim), mu)
  
  taus = outcome_mean(1,sim$G,sim$game1,sim$game2,interf)-outcome_mean(0,sim$G,sim$game1,sim$game2,interf)
  # for (x1 in 0:1)
  # {
  #   for (x2 in 0:1)
  #   {
  #     i_xx = which(sim$game1==x1 & sim$game2==x2)
  #     tau = tau + 
  #   }
  # }
  
  # Naive estimator
  naive_est = mean(sim$Y[which(sim$Z==1)]) - mean(sim$Y[which(sim$Z==0)])
  naive_err[n] = naive_est - mean(taus)
  
  # Regression estimator
  lin_fit = lm(Y~Z+game1+game2, data=sim)
  z1 = data.frame(Z=1, game1=sim$game1, game2=sim$game2)
  z0 = data.frame(Z=0, game1=sim$game1, game2=sim$game2)
  reg_err[n] = mean(predict(lin_fit,z1) - predict(lin_fit,z0) - taus)
}

bias_naive = mean(naive_err)
rmse_naive = sqrt(mean(naive_err^2))

bias_reg = mean(reg_err)
rmse_reg = sqrt(mean(reg_err^2))
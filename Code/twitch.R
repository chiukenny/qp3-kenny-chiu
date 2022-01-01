# Games -> Self-promotion -> Views?

# https://snap.stanford.edu/data/twitch_gamers.html
# https://arxiv.org/pdf/1909.13021.pdf
# Nodes: 7126
# Features: 3169

# Top features:
#   feats Freq
#1   920 6742
#2   224 3816
#3   569 3710
#4   861 3242
#5  3152 3180
#6  2645 2461

library(rjson)
library(tidyverse)
library(Matrix)
#library(MASS)

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


# Simulation

sim = dat
sim$Ngame1 = 0
sim$Ngame2 = 0
sim$N = 0

for (i in 1:nrow(sim))
{
  id = sim$id[i]
  ii = which(dat_edMat[id,])
  sim$Ngame1[i] = mean(sim$game1[ii])
  sim$Ngame2[i] = mean(sim$game2[ii])
  sim$N[i] = length(ii)
}

delta = 2

# Scenario 1
expit = function(x) {1 / (1+exp(-x))}
prop1 = expit(-3 + 3*sim$game1 + 5*sim$game2)
sim$Z = rbernoulli(nrow(sim), prop1)*1

# Function to compute neighbourhood covariates
nbr_trt = function(x)
{
  #return(mean(sim$Z[which(dat_edMat[x,])]))
  return(sum(sim$Z[which(dat_edMat[x,])]))
}
sim$G = sapply(sim$id, nbr_trt)

mu = 5 + 6*(prop1>=0.7) + 10*sim$Z - 3*sim$Z*(prop1>=0.7) + delta*sim$G
sim$Y = rnorm(nrow(sim), mu)

#i00 = which(sim$game1==0 & sim$game2==0)
#i01 = which(sim$game1==0 & sim$game2==1)
#i10 = which(sim$game1==1 & sim$game2==0)
#i11 = which(sim$game1==1 & sim$game2==1)

# Table 1

trt_est = function(sim)
{
  obs_g = sort(unique(sim$G))
  
  sim_1 = sim[which(sim$Z==1),]
  sim_0 = sim[which(sim$Z==0),]
  sim_1g = sim_1
  sim_0g = sim_0
  
  tau = 0
  for (g in obs_g)
  {
    P_1g = length(which(sim_1$G==g)) / nrow(sim_1)
    P_0g = length(which(sim_0$G==g)) / nrow(sim_0)
    sim_1g = sim_1g[which(sim_1g$N >= g),]
    sim_0g = sim_0g[which(sim_0g$N >= g),]
    
    Y1 = sim_1g$Y[which(sim_1g$G==g)]
    Y0 = sim_0g$Y[which(sim_0g$G==g)]
    if (length(Y1)>0) {tau = tau + mean(Y1)*P_1g}
    if (length(Y0)>0) {tau = tau - mean(Y0)*P_0g}
  }
  return(tau)
}

trt_est_Xind = function(sim)
{
  obs_g = sort(unique(sim$G))
  
  tau = 0
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      i = which(sim$game1==x1 & sim$game2==x2)
      sim_x = sim[i,]
      sim_1x = sim_x[which(sim_x$Z==1),]
      sim_0x = sim_x[which(sim_x$Z==0),]
      sim_1xg = sim_1x
      sim_0xg = sim_0x
      
      trt_x = 0
      for (g in obs_g)
      {
        P_1xg = length(which(sim_1x$G==g)) / nrow(sim_1x)
        P_0xg = length(which(sim_0x$G==g)) / nrow(sim_0x)
        sim_1xg = sim_1xg[which(sim_1xg$N >= g),]
        sim_0xg = sim_0xg[which(sim_0xg$N >= g),]
        
        Y1 = sim_1xg$Y[which(sim_1xg$G==g)]
        Y0 = sim_0xg$Y[which(sim_0xg$G==g)]
        if (length(Y1)>0) {trt_x = trt_x + mean(Y1)*P_1xg}
        if (length(Y0)>0) {trt_x = trt_x - mean(Y0)*P_0xg}
      }
      tau = tau + trt_x*length(i)/nrow(sim)
    }
  }
  return(tau)
}

tau1 = 0
for (x1 in 0:1)
{
  for (x2 in 0:1)
  {
    tau1 = tau1 + (10-3*expit(-3+3*x1+5*x2)) * length(which(sim$game1==x1 & sim$game2==x2))/nrow(sim)
  }
}

# Bias of estimator that does not adjust for any covariates
tau_nocov = trt_est(sim)
bias_nocov = tau_nocov - tau1

# Bias of estimator that adjusts for X_ind
tau_Xind = trt_est_Xind(sim)
bias_Xind = tau_Xind - tau1




# Proposed method (not required for tasks?)
# 
# ind_prop_fit = glm(Z~game1*game2, data=sim, family=binomial)
# #summary(ind_prop_fit)
# sim$ind_prop = ind_prop_fit$fitted.values
# # Subclasses are just each level of (game1,game2)
# 
# nbr_prop00_fit = glm(G~Z+Ngame1*Ngame2*N, data=sim[i00,], family=quasibinomial, weights=N)
# #summary(nbr_prop00_fit)
# nbr_prop01_fit = glm(G~Z+Ngame1*Ngame2*N, data=sim[i01,], family=quasibinomial, weights=N)
# #summary(nbr_prop01_fit)
# nbr_prop10_fit = glm(G~Z+Ngame1*Ngame2*N, data=sim[i10,], family=quasibinomial, weights=N)
# #summary(nbr_prop10_fit)
# nbr_prop11_fit = glm(G~Z+Ngame1*Ngame2*N, data=sim[i11,], family=quasibinomial, weights=N)
# #summary(nbr_prop11_fit)
# sim$nbr_prop = 0
# sim$nbr_prop[i00] = nbr_prop00_fit$fitted.values
# sim$nbr_prop[i01] = nbr_prop01_fit$fitted.values
# sim$nbr_prop[i10] = nbr_prop10_fit$fitted.values
# sim$nbr_prop[i11] = nbr_prop11_fit$fitted.values
# 
# # TODO: fix/adjust for v_G
# outcome00_fit = lm(Y~Z*G*nbr_prop, data=sim[i00,])
# #summary(outcome00_fit)
# mu00 = mean(outcome00_fit$fitted.values)
# outcome01_fit = lm(Y~Z*G*nbr_prop, data=sim[i01,])
# #summary(outcome01_fit)
# mu01 = mean(outcome01_fit$fitted.values)
# outcome10_fit = lm(Y~Z*G*nbr_prop, data=sim[i10,])
# #summary(outcome10_fit)
# mu10 = mean(outcome10_fit$fitted.values)
# outcome11_fit = lm(Y~Z*G*nbr_prop, data=sim[i11,])
# #summary(outcome11_fit)
# mu11 = mean(outcome11_fit$fitted.values)
# 
# mu_hat = mu00*length(i00)

# Extension: G as sum rather than mean?
#nbr_prop00_fit = glm.nb(G~Z+Ngame1*Ngame2*N, data=sim[i00,])
#summary(nbr_prop00_fit)
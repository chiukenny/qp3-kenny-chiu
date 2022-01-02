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

delta = 2 #0.5

# Scenario 1
expit = function(x) {1 / (1+exp(-x))}
prop1 = expit(-3 + 3*sim$game1 + 5*sim$game2)
sim$Z = rbernoulli(nrow(sim), prop1)*1

# Function to compute neighbourhood covariates
nbr_trt = function(x)
{
  return(mean(sim$Z[which(dat_edMat[x,])]))
  #return(sum(sim$Z[which(dat_edMat[x,])]))
}
sim$G = sapply(sim$id, nbr_trt)

mu = 5 + 6*(prop1>=0.7) + 10*sim$Z - 3*sim$Z*(prop1>=0.7) + delta*sim$G
sim$Y = rnorm(nrow(sim), mu)

#i00 = which(sim$game1==0 & sim$game2==0)
#i01 = which(sim$game1==0 & sim$game2==1)
#i10 = which(sim$game1==1 & sim$game2==0)
#i11 = which(sim$game1==1 & sim$game2==1)


# Table 1

# trt_est = function(sim)
# {
#   obs_g = sort(unique(sim$G))
#   
#   sim_1 = sim[which(sim$Z==1),]
#   sim_0 = sim[which(sim$Z==0),]
#   n_Z1 = nrow(sim_1)
#   n_Z0 = nrow(sim_0)
#   
#   tau = 0
#   sim_V1 = sim_1
#   sim_V0 = sim_0
#   for (g in obs_g)
#   {
#     sim_V1 = sim_V1[which(sim_V1$N >= g),]
#     sim_V0 = sim_V0[which(sim_V0$N >= g),]
#     Y_1g = sim_V1$Y[which(sim_V1$G==g)]
#     Y_0g = sim_V0$Y[which(sim_V0$G==g)]
#     P_1g = length(Y_1g) / n_Z1
#     P_0g = length(Y_0g) / n_Z0
#     
#     if (length(Y_1g)>0) {tau = tau + mean(Y_1g)*P_1g}
#     if (length(Y_0g)>0) {tau = tau - mean(Y_0g)*P_0g}
#   }
#   return(tau)
# }

# Theorem 2B (eqn. 12)
bias_T2B = function(sim)
{
  # Take g' = 0, u' = (0,0)
  obs_g = sort(unique(sim$G))
  
  Y_1000 = sim$Y[which(sim$Z==1 & sim$G==0 & sim$game1==0 & sim$game2==0)]
  Y_0000 = sim$Y[which(sim$Z==0 & sim$G==0 & sim$game1==0 & sim$game2==0)]
  if (length(Y_1000)>0) {Y_1000 = mean(Y_1000)} else {Y_1000 = 0}
  if (length(Y_0000)>0) {Y_0000 = mean(Y_0000)} else {Y_0000 = 0}
  
  n = nrow(sim)
  n_Z1 = length(which(sim$Z==1))
  n_Z0 = length(which(sim$Z==0))
  
  bias = 0
  sim_V = sim
  for (g in obs_g)
  {
    sim_V = sim_V[which(sim_V$N>=g),]
    sim_Vg = sim_V[which(sim_V$G==g),]
    sim_V1g = sim_Vg[which(sim_Vg$Z==1),]
    sim_V0g = sim_Vg[which(sim_Vg$Z==0),]
    
    Pg = nrow(sim_Vg) / n
    Pg_1 = nrow(sim_V1g) / n_Z1
    Pg_0 = nrow(sim_V0g) / n_Z0
    
    for (u1 in 0:1)
    {
      for (u2 in 0:1)
      {
        if (g==0 & u1==0 & u2==0) {next}
        Y_V1gu = sim_V1g$Y[which(sim_V1g$game1==u1 & sim_V1g$game2==u2)]
        Y_V0gu = sim_V0g$Y[which(sim_V0g$game1==u1 & sim_V0g$game2==u2)]
        Pu_1g = length(Y_V1gu) / nrow(sim_V1g)
        Pu_0g = length(Y_V0gu) / nrow(sim_V0g)
        Pu_g = length(which(sim_Vg$game1==u1 & sim_Vg$game2==u2)) / nrow(sim_Vg)
        if (length(Y_V1gu)>0)
        {
          bias = bias + (mean(Y_V1gu)-Y_1000)*(Pu_1g*Pg_1-Pu_g*Pg)
        }
        if (length(Y_V0gu)>0)
        {
          bias = bias - (mean(Y_V0gu)-Y_0000)*(Pu_0g*Pg_0-Pu_g*Pg)
        }
      }
    }
  }
  return(bias)
}

# trt_est_Xind = function(sim)
# {
#   obs_g = sort(unique(sim$G))
#   
#   tau = 0
#   for (x1 in 0:1)
#   {
#     for (x2 in 0:1)
#     {
#       sim_x = sim[which(sim$game1==x1 & sim$game2==x2),]
#       sim_1x = sim_x[which(sim_x$Z==1),]
#       sim_0x = sim_x[which(sim_x$Z==0),]
#       sim_V1x = sim_1x
#       sim_V0x = sim_0x
#       
#       tau_x = 0
#       for (g in obs_g)
#       {
#         sim_V1x = sim_V1x[which(sim_V1x$N >= g),]
#         sim_V0x = sim_V0x[which(sim_V0x$N >= g),]
#         Y_V1gx = sim_V1x$Y[which(sim_V1x$G==g)]
#         Y_V0gx = sim_V0x$Y[which(sim_V0x$G==g)]
#         Pg_1x = length(Y_V1gx) / nrow(sim_1x)
#         Pg_0x = length(Y_V0gx) / nrow(sim_0x)
#         
#         if (length(Y_V1gx)>0) {tau_x = tau_x + mean(Y_V1gx)*Pg_1x}
#         if (length(Y_V0gx)>0) {tau_x = tau_x - mean(Y_V0gx)*Pg_0x}
#       }
#       tau = tau + tau_x*nrow(sim_x)/nrow(sim)
#     }
#   }
#   return(tau)
# }

bias_C2 = function(sim)
{
  # Take g' = 0
  obs_g = sort(unique(sim$G))
  
  bias = 0
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      sim_x = sim[which(sim$game1==x1 & sim$game2==x2),]
      
      Y_0x = sim_x$Y[which(sim_x$G==0)]
      if (length(Y_0x)>0) {Y_0x = mean(Y_0x)} else {Y_0x = 0}
      
      n_x = nrow(sim_x)
      n_1 = length(which(sim_x$Z==1))
      n_0 = n_x - n_1
      
      bias_x = 0
      sim_Vx = sim_x
      for (g in obs_g[-1])
      {
        sim_Vx = sim_Vx[which(sim_Vx$N>=g),]
        sim_Vgx = sim_Vx[which(sim_Vx$G==g),]
        Y_Vgx = sim_Vgx$Y
        
        if (length(Y_Vgx)>0)
        {
          n_1g = length(which(sim_Vgx$Z==1))
          Pg_1x = n_1g / n_1
          Pg_0x = (nrow(sim_Vgx)-n_1g) / n_0
          bias_x = bias_x + (mean(Y_Vgx)-Y_0x)*(Pg_1x-Pg_0x)
        }
      }
      bias = bias + bias_x*n_x/nrow(sim)
    }
  }
  return(bias)
}


# tau1 = 0
# for (x1 in 0:1)
# {
#   for (x2 in 0:1)
#   {
#     #tau1 = tau1 + (10-3*(expit(-3+3*x1+5*x2)>=0.7))*length(which(sim$game1==x1 & sim$game2==x2))
#     for (g in unique(sim$G))
#     {
#       Y1 = sim$Y[which(sim$game1==x1 & sim$game2==x2 & sim$G==g & sim$Z==1)]
#       Y0 = sim$Y[which(sim$game1==x1 & sim$game2==x2 & sim$G==g & sim$Z==0)]
#       if (length(Y1)>0)
#       {
#         tau1 = tau1 + mean(Y1)*length(which(sim$G==g))*length(which(sim$game1==x1 & sim$game2==x2)) / nrow(sim)^2
#       }
#       if (length(Y0)>0)
#       {
#         tau1 = tau1 - mean(Y0)*length(which(sim$G==g))*length(which(sim$game1==x1 & sim$game2==x2)) / nrow(sim)^2
#       }
#     }
#   }
# }
#tau1 = tau1 / nrow(sim)^2

# Bias of estimator that does not adjust for any covariates
#tau_nocov = trt_est(sim)
#bias_nocov = tau_nocov - tau1
bias_t2b = bias_T2B(sim)

# Bias of estimator that adjusts for X_ind
#tau_Xind = trt_est_Xind(sim)
#bias_Xind = tau_Xind - tau1
bias_c2 = bias_C2(sim)



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
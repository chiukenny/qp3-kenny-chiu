# Games -> Self-promotion -> Views?

# https://snap.stanford.edu/data/twitch-social-networks.html
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


# Parameters
interf = 4 # {4,6,8} mean
#interf = 0.5 # sum


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

sim = dat
# Mean
#sim$Ngame1 = as.vector(dat_edMat %*% sim$game1) / rowSums(dat_edMat)
#sim$Ngame2 = as.vector(dat_edMat %*% sim$game2) / rowSums(dat_edMat)
# Sum
sim$Ngame1 = as.vector(dat_edMat %*% sim$game1)
sim$Ngame2 = as.vector(dat_edMat %*% sim$game2)
sim$N = rowSums(dat_edMat)

# Scenario 1
prop = ind_prop(sim$game1, sim$game2)
sim$Z = rbernoulli(nrow(sim), prop)*1

sim$G = as.vector(dat_edMat %*% sim$Z) / rowSums(dat_edMat) # Mean
#sim$G = as.vector(dat_edMat %*% sim$Z) # Sum

#mu = 5 + 6*(prop1>=0.7) + 10*sim$Z - 3*sim$Z*(prop1>=0.7) + delta*sim$G
mu = outcome_mean(sim$Z,sim$G,sim$game1,sim$game2,interf)
sim$Y = rnorm(nrow(sim), mu)


# Table 1

# Theorem 2B (eqn. 12)

# Computed using generating model
# bias_T2B = function(sim, beta)
# {
#   # Take g' = 0, u' = (0,0)
#   n = nrow(sim)
#   n_Z1 = length(which(sim$Z==1))
#   n_Z0 = length(which(sim$Z==0))
#   
#   bias = 0
#   sim_V = sim
#   obs_g = sort(unique(sim$G))
#   for (g in obs_g)
#   {
#     sim_V = sim_V[which(sim_V$N>=g),]
#     sim_Vg = sim_V[which(sim_V$G==g),]
#     sim_V1g = sim_Vg[which(sim_Vg$Z==1),]
#     sim_V0g = sim_Vg[which(sim_Vg$Z==0),]
#     
#     Pg = nrow(sim_Vg) / n
#     Pg_1 = nrow(sim_V1g) / n_Z1
#     Pg_0 = nrow(sim_V0g) / n_Z0
#     
#     for (u1 in 0:1)
#     {
#       for (u2 in 0:1)
#       {
#         Pu_g = length(which(sim_Vg$game1==u1 & sim_Vg$game2==u2)) / nrow(sim_Vg)
#         if (Pu_g == 0) {next}
#         Pu_1g = length(which(sim_V1g$game1==u1 & sim_V1g$game2==u2)) / nrow(sim_V1g)
#         Pu_0g = length(which(sim_V0g$game1==u1 & sim_V0g$game2==u2)) / nrow(sim_V0g)
#         if (is.nan(Pu_1g)) {Pu_1g = 0}
#         if (is.nan(Pu_0g)) {Pu_0g = 0}
#         bias = bias + (outcome_mean(1,g,u1,u2,beta)-outcome_mean(1,0,0,0,beta))*(Pu_1g*Pg_1-Pu_g*Pg)
#         bias = bias - (outcome_mean(0,g,u1,u2,beta)-outcome_mean(0,0,0,0,beta))*(Pu_0g*Pg_0-Pu_g*Pg)
#       }
#     }
#   }
#   return(bias)
# }

bias_T2B = function(sim, beta)
{
  # Take g' = 0, u' = (0,0)
  n = nrow(sim)
  i_1 = which(sim$Z==1)
  n_1 = length(i_1)
  i_0 = sim$id[-i_1]
  n_0 = n - n_1
  bias = 0
  for (g in unique(sim$G))
  {
    i_g = which(sim$G==g)
    n_g = length(i_g)
    
    i_1g = intersect(i_1, i_g)
    n_1g = length(i_1g)
    i_0g = intersect(i_0, i_g)
    n_0g = n_g - n_1g
    for (u1 in 0:1)
    {
      for (u2 in 0:1)
      {
        i_uu = which(sim$game1==u1 & sim$game2==u2)
        
        Pug_1 = length(intersect(i_1g, i_uu)) / n_1 #/ n_1g
        #Pg_1 = n_1g / n_1
        Pug_0 = length(intersect(i_0g, i_uu)) / n_0
        Pug = length(intersect(i_g, i_uu)) / n
        bias = bias + (outcome_mean(1,g,u1,u2,beta)-outcome_mean(1,0,0,0,beta))*(Pug_1-Pug) - (outcome_mean(0,g,u1,u2,beta)-outcome_mean(0,0,0,0,beta))*(Pug_0-Pug)
      }
    }
  }
  return(bias)
}


# Corollary 2
# Computed using generating model
# bias_C2 = function(sim, beta)
# {
#   # Take g' = 0
#   bias = 0
#   #counter = 0
#   total_p = 0
#   for (x1 in 0:1)
#   {
#     for (x2 in 0:1)
#     {
#       sim_x = sim[which(sim$game1==x1 & sim$game2==x2),]
#       
#       n_x = nrow(sim_x)
#       n_1x = length(which(sim_x$Z==1))
#       n_0x = n_x - n_1x
#       
#       bias_x = 0
#       obs_g = sort(unique(sim_x$G))
#       sim_Vx = sim_x
#       for (g in obs_g)
#       {
#         #spl = split(sim_Vx, sim_Vx$G==g)
#         #sim_Vx = spl["FALSE"]
#         #sim_Vgx = spl["TRUE"]
#         sim_Vx = sim_Vx[which(sim_Vx$N>=g),]
#         Z_gx = sim_Vx$Z[which(sim_Vx$G==g)]
#         
#         n_1gx = length(which(Z_gx==1))
#         Pg_1x = n_1gx / n_1x
#         Pg_0x = (length(Z_gx)-n_1gx) / n_0x
#         #bias_x = bias_x + beta*g*(Pg_1x-Pg_0x)
#         # Z irrelevant in bias computation
#         bias_x = bias_x + (outcome_mean(0,g,x1,x2,beta)-outcome_mean(0,0,x1,x2,beta))*(Pg_1x-Pg_0x)
#         #counter = counter + nrow(sim_Vgx)
#       }
#       total_p = total_p + n_x/nrow(sim)
#       bias = bias + bias_x*n_x#/nrow(sim)
#     }
#   }
#   #message(counter)
#   message(total_p)
#   return(bias/nrow(sim))
# }
bias_C2 = function(sim, beta)
{
  # Take g' = 0
  i_1 = which(sim$Z==1)
  i_0 = sim$id[-i_1]
  bias = 0
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      i_xx = which(sim$game1==x1 & sim$game2==x2)
      n_xx = length(i_xx)
      
      i_1xx = intersect(i_1, i_xx)
      n_1xx = length(i_1xx)
      i_0xx = intersect(i_0, i_xx)
      n_0xx = n_xx - n_1xx
      
      bias_xx = 0
      for (g in unique(sim$G))
      {
        Pg_1xx = length(which(sim$G[i_1xx]==g)) / n_1xx
        Pg_0xx = length(which(sim$G[i_0xx]==g)) / n_0xx
        bias_xx = bias_xx + (outcome_mean(0,g,x1,x2,beta)-outcome_mean(0,0,x1,x2,beta))*(Pg_1xx-Pg_0xx)
      }
      bias = bias + bias_xx*n_xx
    }
  }
  return(bias/nrow(sim))
}


# Include Ngame1, Ngame2, and N as covariates
# bias_C2_Xn = function(sim, beta)
# {
#   # Take g' = 0
#   bias = 0
#   for (x1 in 0:1)
#   {
#     for (x2 in 0:1)
#     {
#       sim_x = sim[which(sim$game1==x1 & sim$game2==x2),]
# 
#       obs_N = sort(unique(sim_x$N))
#       for (N in obs_N)
#       {
#         sim_Nx = sim_x[which(sim_x$N==N),]
# 
#         obs_n1 = sort(unique(sim_Nx$Ngame1))
#         for (n1 in obs_n1)
#         {
#           #sim_xN = sim_xN[which(sim_xN$N>=n1),]
#           sim_xn = sim_Nx[which(sim_Nx$Ngame1==n1),]
# 
#           obs_n2 = sort(unique(sim_xn$Ngame2))
#           for (n2 in obs_n2)
#           {
#             #sim_xn = sim_xn[which(sim_xn$N>=n2),]
#             sim_xnn = sim_xn[which(sim_xn$Ngame2==n2),]
#             #if (nrow(sim_xnn)==0) {next}
# 
#             n_xnn = nrow(sim_xnn)
#             n_1xnn = length(which(sim_xnn$Z==1))
#             n_0xnn = n_xnn - n_1xnn
# 
#             bias_xnn = 0
#             obs_g = sort(unique(sim_xnn$G))
#             #sim_Vxnn = sim_xnn
#             for (g in obs_g)
#             {
#               #sim_Vxnn = sim_Vxnn[which(sim_Vxnn$N>=g),]
#               #sim_Vgxnn = sim_Vxnn[which(sim_Vxnn$G==g),]
#               sim_Vgxnn = sim_xnn[which(sim_xnn$G==g),]
# 
#               n_1gxnn = length(which(sim_Vgxnn$Z==1))
#               p = 0
#               if (n_1xnn>0) {p = p + n_1gxnn/n_1xnn}
#               if (n_0xnn>0) {p = p - (nrow(sim_Vgxnn)-n_1gxnn)/n_0xnn}
#               bias_xnn = bias_xnn + (outcome_mean(0,g,x1,x2,beta)-outcome_mean(0,0,x1,x2,beta))*p
#             }
#             bias = bias + bias_xnn*n_xnn#/nrow(sim)
#           }
#         }
#       }
#     }
#   }
#   return(bias/nrow(sim))
# }

bias_C2_Xn = function(sim, beta)
{
  # Take g' = 0
  bias = 0
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      i_xx = which(sim$game1==x1 & sim$game2==x2)

      for (N in unique(sim$N[i_xx]))
      {
        i_xxN = which(sim$game1==x1 & sim$game2==x2 & sim$N==N)
        sim_xxN = sim[i_xxN,]
        for (n1 in unique(sim_xxN$Ngame1))
        {
          for (n2 in unique(sim_xxN$Ngame2))
          {
            i_xxnnN = which(sim_xxN$Ngame1==n1 & sim_xxN$Ngame2==n2)
            n_xxnnN = length(i_xxnnN)
            if (n_xxnnN==0) {next}
            i_1xxnnN = intersect(i_xxnnN,which(sim_xxN$Z==1))
            n_1xxnnN = length(i_1xxnnN)
            i_0xxnnN = intersect(i_xxnnN,which(sim_xxN$Z==0))
            n_0xxnnN = length(i_0xxnnN)

            bias_xxnnN = 0
            for (g in unique(sim_xxN$G[i_xxnnN]))
            {
              i_gxxnnN = intersect(i_xxnnN,which(sim_xxN$G==g))
              p = 0
              if (n_1xxnnN>0) {p = p+length(intersect(i_1xxnnN,i_gxxnnN))/n_1xxnnN}
              if (n_0xxnnN>0) {p = p-length(intersect(i_0xxnnN,i_gxxnnN))/n_0xxnnN}
              bias_xxnnN = bias_xxnnN + (outcome_mean(0,g,x1,x2,beta)-outcome_mean(0,0,x1,x2,beta))*p
            }
            bias = bias + bias_xxnnN*n_xxnnN
          }
        }
      }
    }
  }
  return(bias/nrow(sim))
}

bias_C2N = function(dat, beta)
{
  bias = 0
  # Subset by (game1, game2, N, Ngame1, Ngame2)
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      for (N in unique(dat$N[dat$game1==x1 & dat$game2==x2]))
      {
        i_xxN = dat$game1==x1 & dat$game2==x2 & dat$N==N
        dat_xxN = dat[i_xxN,] # Save subsetted data for performance
        i_1xxN = dat_xxN$Z==1
        
        for (n1 in unique(dat_xxN$Ngame1))
        {
          for (n2 in unique(dat_xxN$Ngame2))
          {
            i_xxnnN = dat_xxN$Ngame1==n1 & dat_xxN$Ngame2==n2
            n_xxnnN = sum(i_xxnnN)
            if (n_xxnnN==0) {next}
            i_1xxnnN = i_xxnnN & i_1xxN
            n_1xxnnN = sum(i_1xxnnN)
            i_0xxnnN = i_xxnnN & !i_1xxN
            n_0xxnnN = n_xxnnN - n_1xxnnN
            
            bias_xxnnN = 0
            for (g in unique(dat_xxN$G[i_xxnnN]))
            {
              i_gxxN = dat_xxN$G==g
              
              # Compute empirical probabilities
              pg_1xxnnN = 0
              pg_0xxnnN = 0
              pg_xxnnN = sum(i_xxnnN & i_gxxN) / n_xxnnN
              if (n_1xxnnN>0) {pg_1xxnnN = sum(i_1xxnnN & i_gxxN)/n_1xxnnN}
              if (n_0xxnnN>0) {pg_0xxnnN = sum(i_0xxnnN & i_gxxN)/n_0xxnnN}
              
              # Take g'=0
              bias_xxnnN = bias_xxnnN +
                (outcome_mean(1,g,x1,x2,beta)-outcome_mean(1,0,x1,x2,beta))*(pg_1xxnnN-pg_xxnnN) - (outcome_mean(0,g,x1,x2,beta)-outcome_mean(0,0,x1,x2,beta))*(pg_0xxnnN-pg_xxnnN)
            }
            bias = bias + bias_xxnnN*n_xxnnN
          }
        }
      }
    }
  }
  return(bias/nrow(dat))
}


# n_1xxnnN = 0
# n_0xxnnN = 0
# for (g in unique(sim$G[i_xx]))
# {
#   i_gxx = which(sim$game1==x1 & sim$game2==x2 & sim$G==g)
#   
#   i_1gxx = intersect(i_1, i_gxx)
#   i_0gxx = intersect(i_0, i_gxx)
#   
#   bias_g = 0
#   sim_1gxx = sim[i_1gxx,]
#   for (N in unique(sim$N[i_1gxx]))
#   {
#     for (n1 in unique(sim$Ngame1[i_1gxx]))
#     {
#       for (n2 in unique(sim$Ngame2[i_1gxx]))
#       {
#         i_1gxxnnN = which(sim_1gxx$N==N & sim_1gxx$Ngame==n1 & sim_1gxx$Ngame==n2)
#         if (length(i_1gxxnnN)==0) {next}
#         
#       }
#     }
#   }
# }


# Bias of estimator that does not adjust for any covariates
#tau_nocov = trt_est(sim)
#bias_nocov = tau_nocov - tau1
bias_t2b = bias_T2B(sim, interf)

# Bias of estimator that adjusts for X_ind
#tau_Xind = trt_est_Xind(sim)
#bias_Xind = tau_Xind - tau1
bias_c2 = bias_C2(sim, interf)

# Bias of estimator that adjusts for X_z
bias_c2_xn = bias_C2_Xn(sim, interf)
# TODO: should equal bias_c2?




# TODO: DELETE

# g vs g'
# bias_T2B = function(sim)
# {
#   # Take g' = 0, u' = (0,0)
#   obs_g = sort(unique(sim$G))
# 
#   Y_1000 = sim$Y[which(sim$Z==1 & sim$G==0 & sim$game1==0 & sim$game2==0)]
#   Y_0000 = sim$Y[which(sim$Z==0 & sim$G==0 & sim$game1==0 & sim$game2==0)]
#   if (length(Y_1000)>0) {Y_1000 = mean(Y_1000)} else {Y_1000 = 0}
#   if (length(Y_0000)>0) {Y_0000 = mean(Y_0000)} else {Y_0000 = 0}
# 
#   n = nrow(sim)
#   n_Z1 = length(which(sim$Z==1))
#   n_Z0 = length(which(sim$Z==0))
# 
#   bias = 0
#   sim_V = sim
#   for (g in obs_g)
#   {
#     sim_V = sim_V[which(sim_V$N>=g),]
#     sim_Vg = sim_V[which(sim_V$G==g),]
#     sim_V1g = sim_Vg[which(sim_Vg$Z==1),]
#     sim_V0g = sim_Vg[which(sim_Vg$Z==0),]
# 
#     Pg = nrow(sim_Vg) / n
#     Pg_1 = nrow(sim_V1g) / n_Z1
#     Pg_0 = nrow(sim_V0g) / n_Z0
# 
#     for (u1 in 0:1)
#     {
#       for (u2 in 0:1)
#       {
#         if (g==0 & u1==0 & u2==0) {next}
#         Y_V1gu = sim_V1g$Y[which(sim_V1g$game1==u1 & sim_V1g$game2==u2)]
#         Y_V0gu = sim_V0g$Y[which(sim_V0g$game1==u1 & sim_V0g$game2==u2)]
#         Pu_1g = length(Y_V1gu) / nrow(sim_V1g)
#         Pu_0g = length(Y_V0gu) / nrow(sim_V0g)
#         Pu_g = length(which(sim_Vg$game1==u1 & sim_Vg$game2==u2)) / nrow(sim_Vg)
#         if (length(Y_V1gu)>0)
#         {
#           bias = bias + (mean(Y_V1gu)-Y_1000)*(Pu_1g*Pg_1-Pu_g*Pg)
#         }
#         if (length(Y_V0gu)>0)
#         {
#           bias = bias - (mean(Y_V0gu)-Y_0000)*(Pu_0g*Pg_0-Pu_g*Pg)
#         }
#       }
#     }
#   }
#   return(bias)
# }

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

# Keeping dependence on g
# bias_C2 = function(sim)
# {
#   bias = 0
#   for (x1 in 0:1)
#   {
#     for (x2 in 0:1)
#     {
#       sim_x = sim[which(sim$game1==x1 & sim$game2==x2),]
#       
#       obs_g = sort(unique(sim_x$G))
#       
#       n_x = nrow(sim_x)
#       n_1x = length(which(sim_x$Z==1))
#       n_0x = n_x - n_1x
# 
#       bias_x = 0
#       sim_Vx = sim_x
#       for (g in obs_g)
#       {
#         sim_Vx = sim_Vx[which(sim_Vx$N>=g),]
#         sim_Vgx = sim_Vx[which(sim_Vx$G==g),]
#         Y_V1gx = sim_Vgx$Y[which(sim_Vgx$Z==1)]
#         Y_V0gx = sim_Vgx$Y[which(sim_Vgx$Z==0)]
#         
#         bias_g = 0
#         if (length(Y_V1gx)>0) {bias_g = bias_g + mean(Y_V1gx)}
#         if (length(Y_V0gx)>0) {bias_g = bias_g - mean(Y_V0gx)}
#         Pg_1x = length(Y_V1gx) / n_1x
#         Pg_0x = length(Y_V0gx) / n_0x
#         #message(paste(bias_g,Pg_1x-Pg_0x,bias_g*(Pg_1x-Pg_0x)))
#         bias_x = bias_x + bias_g*(Pg_1x-Pg_0x)
#       }
#       bias = bias + bias_x*n_x/nrow(sim)
#     }
#   }
#   return(bias)
# }

# Ignores dependence on G
# bias_C2 = function(sim)
# {
#   obs_g = sort(unique(sim$G))
#   
#   bias = 0
#   for (x1 in 0:1)
#   {
#     for (x2 in 0:1)
#     {
#       sim_1x = sim[which(sim$game1==x1 & sim$game2==x2 & sim$Z==1),]
#       sim_0x = sim[which(sim$game1==x1 & sim$game2==x2 & sim$Z==0),]
#       
#       n_1 = nrow(sim_1x)
#       n_0 = nrow(sim_0x)
#       
#       bias_x = 0
#       sim_V1x = sim_1x
#       sim_V0x = sim_0x
#       for (g in obs_g)
#       {
#         sim_V1x = sim_V1x[which(sim_V1x$N>=g),]
#         sim_V0x = sim_V0x[which(sim_V0x$N>=g),]
#         
#         Pg_1x = length(which(sim_V1x$G==g)) / n_1
#         Pg_0x = length(which(sim_V0x$G==g)) / n_0
#         bias_x = bias_x + Pg_1x - Pg_0x
#       }
#       bias = bias + bias_x*(mean(sim_1x$Y)-mean(sim_0x$Y))*(n_1+n_0)/nrow(sim)
#     }
#   }
#   return(bias)
# }

# g vs g'
# This computation seems most likely following Theorem 2B
# Computation of expectation ignores value of Z
# bias_C2 = function(sim)
# {
#   # Take g' = 0
#   obs_g = sort(unique(sim$G))
# 
#   bias = 0
#   for (x1 in 0:1)
#   {
#     for (x2 in 0:1)
#     {
#       sim_x = sim[which(sim$game1==x1 & sim$game2==x2),]
# 
#       Y_0x = sim_x$Y[which(sim_x$G==0)]
#       if (length(Y_0x)>0) {Y_0x = mean(Y_0x)} else {Y_0x = 0}
# 
#       n_x = nrow(sim_x)
#       n_1 = length(which(sim_x$Z==1))
#       n_0 = n_x - n_1
# 
#       bias_x = 0
#       sim_Vx = sim_x
#       for (g in obs_g[-1])
#       {
#         sim_Vx = sim_Vx[which(sim_Vx$N>=g),]
#         sim_Vgx = sim_Vx[which(sim_Vx$G==g),]
#         Y_Vgx = sim_Vgx$Y
# 
#         if (length(Y_Vgx)>0)
#         {
#           n_1g = length(which(sim_Vgx$Z==1))
#           Pg_1x = n_1g / n_1
#           Pg_0x = (nrow(sim_Vgx)-n_1g) / n_0
#           bias_x = bias_x + (mean(Y_Vgx)-Y_0x)*(Pg_1x-Pg_0x)
#         }
#       }
#       bias = bias + bias_x*n_x/nrow(sim)
#     }
#   }
#   return(bias)
# }

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





bias_C2N_noN = function(dat, beta)
{
  bias = 0
  # Subset by (game1, game2, N, Ngame1, Ngame2)
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      i_xx = dat$game1==x1 & dat$game2==x2
      dat_xx = dat[i_xx,] # Save subsetted data for performance
      i_1xx = dat_xx$Z==1
      i_0xx = dat_xx$Z==0
      
      for (n1 in unique(dat_xx$Ngame1))
      {
        for (n2 in unique(dat_xx$Ngame2))
        {
          i_xxnn = dat_xx$Ngame1==n1 & dat_xx$Ngame2==n2
          n_xxnn = sum(i_xxnn)
          if (n_xxnn==0) {next}
          i_1xxnn = i_xxnn & i_1xx
          n_1xxnn = sum(i_1xxnn)
          i_0xxnn = i_xxnn & i_0xx
          n_0xxnn = n_xxnn - n_1xxnn
          
          bias_xxnn = 0
          for (g in unique(dat_xx$G[i_xxnn]))
          {
            if (g==0) {next}
            i_gxx = dat_xx$G==g
            
            # Compute empirical probabilities
            p = 0
            if (n_1xxnn>0) {p = p+sum(i_1xxnn & i_gxx)/n_1xxnn}
            if (n_0xxnn>0) {p = p-sum(i_0xxnn & i_gxx)/n_0xxnn}
            
            # Take g'=0
            bias_xxnn = bias_xxnn +
              (outcome_mean(0,g,x1,x2,beta)-outcome_mean(0,0,x1,x2,beta))*p
          }
          bias = bias + bias_xxnn*n_xxnn
        }
      }
    }
  }
  return(bias/nrow(dat))
}

bias_C2N_onlyN = function(dat, beta)
{
  bias = 0
  # Subset by (game1, game2, N, Ngame1, Ngame2)
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      for (N in unique(dat$N[dat$game1==x1 & dat$game2==x2]))
      {
        i_xxN = dat$game1==x1 & dat$game2==x2 & dat$N==N
        dat_xxN = dat[i_xxN,] # Save subsetted data for performance
        i_1xxN = dat_xxN$Z==1
        i_0xxN = dat_xxN$Z==0
        n_1xxN = sum(i_1xxN)
        n_0xxN = sum(i_0xxN)
        
        bias_xxN = 0
        for (g in unique(dat_xxN$G))
        {
          if (g==0) {next}
          i_gxxN = dat_xxN$G==g
          
          # Compute empirical probabilities
          p = 0
          if (n_1xxN>0) {p = p+sum(i_1xxN & i_gxxN)/n_1xxN}
          if (n_0xxN>0) {p = p-sum(i_0xxN & i_gxxN)/n_0xxN}
          
          # Take g'=0
          bias_xxN = bias_xxN +
            (outcome_mean(0,g,x1,x2,beta)-outcome_mean(0,0,x1,x2,beta))*p
        }
        bias = bias + bias_xxN*nrow(dat_xxN)
      }
    }
  }
  return(bias/nrow(dat))
}
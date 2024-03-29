# This file contains self-defined functions that are used
# in the simulation study.
# ---------------------------------------------------------

library(MASS)



# Helper functions
# ---------------------------------------------------------

# Compute expit
expit = function(x) {1 / (1+exp(-x))}

# Compute individual propensity score
ind_prop = function(x1, x2) {expit(-3 + 3*x1 + 4*x2)}

# Compute expected outcome based on model
outcome_mean = function(z, g, x1, x2, gamma)
{
  prop = ind_prop(x1, x2)
  return(5 + 6*(prop>=0.7) + 10*z - 3*z*(prop>=0.7) + gamma*g)
}

# Predicts negative binomial responses given fitted model and new data
predict_nb = function(fit, dat)
{
  if (is.null(fit)) {return(rep(0,nrow(dat)))}
  return(dnbinom(dat$N,
                 mu=predict(fit,newdata=dat,type="response"),
                 size=dat$G))
}

# Compute means and standardized difference of two groups
std_diff = function(v1, v0, binary=T)
{
  mean1 = mean(v1)
  mean0 = mean(v0)
  if (binary)
  {
    sdiff = (mean1-mean0) / sqrt((mean1*(1-mean1)+mean0*(1-mean0))/2)
  } else {
    sdiff = (mean1-mean0) / sqrt((var(v1)+var(v0))/2)
  }
  return(c(mean1, mean0, sdiff))
}



# Bias functions
# ---------------------------------------------------------

# Variable naming:
#   - i_X: indices of nodes that satisfy X
#   - n_X: number of nodes that satisfy X
#   - PY_X: (proportional) probability of Y given X

# Compute bias due to unmeasured confounders (Theorem 2B, eqn.12)
bias_T2B = function(dat, gamma)
{
  n = nrow(dat)
  i_1 = dat$Z==1
  n_1 = sum(i_1)
  i_0 = !i_1
  n_0 = n - n_1
  
  bias = 0
  for (g in unique(dat$G))
  {
    i_g = dat$G==g
    n_g = sum(i_g)
    i_1g = i_1 & i_g
    n_1g = sum(i_1g)
    i_0g = i_0 & i_g
    n_0g = n_g - n_1g
    
    # Subset by unmeasured confounders (game1, game2)
    for (u1 in 0:1)
    {
      for (u2 in 0:1)
      {
        i_uu = dat$game1==u1 & dat$game2==u2
        
        # Compute empirical probabilities
        Pug_1 = sum(i_1g & i_uu) / n_1
        Pug_0 = sum(i_0g & i_uu) / n_0
        Pug = sum(i_g & i_uu) / n
        
        # Take g'=0, u'=(0,0)
        bias = bias +
          (outcome_mean(1,g,u1,u2,gamma)-outcome_mean(1,0,0,0,gamma))*(Pug_1-Pug) -
          (outcome_mean(0,g,u1,u2,gamma)-outcome_mean(0,0,0,0,gamma))*(Pug_0-Pug)
      }
    }
  }
  return(bias)
}

# Compute bias given individual covariates (Corollary 2, eqn.11)
bias_C2I = function(dat, gamma)
{
  i_1 = dat$Z==1
  i_0 = !i_1
  
  bias = 0
  # Subset by individual covariates (game1, game2)
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      i_xx = dat$game1==x1 & dat$game2==x2
      n_xx = sum(i_xx)
      i_1xx = i_1 & i_xx
      n_1xx = sum(i_1xx)
      i_0xx = i_0 & i_xx
      n_0xx = n_xx - n_1xx
      
      bias_xx = 0
      for (g in unique(dat$G))
      {
        if (g==0) {next}
        
        # Compute empirical probabilities
        Pg_1xx = sum(dat$G[i_1xx]==g) / n_1xx
        Pg_0xx = sum(dat$G[i_0xx]==g) / n_0xx
        
        # Take g'=0
        bias_xx = bias_xx +
          (outcome_mean(0,g,x1,x2,gamma)-outcome_mean(0,0,x1,x2,gamma))*(Pg_1xx-Pg_0xx)
      }
      bias = bias + bias_xx*n_xx
    }
  }
  return(bias/nrow(dat))
}

# Compute bias given individual and neighbourhood covariates (Corollary 2, eqn.11)
bias_C2N = function(dat, gamma)
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
        i_0xxN = !i_1xxN
        
        for (n1 in unique(dat_xxN$Ngame1))
        {
          for (n2 in unique(dat_xxN$Ngame2))
          {
            i_xxnnN = dat_xxN$Ngame1==n1 & dat_xxN$Ngame2==n2
            n_xxnnN = sum(i_xxnnN)
            if (n_xxnnN==0) {next}
            i_1xxnnN = i_xxnnN & i_1xxN
            n_1xxnnN = sum(i_1xxnnN)
            i_0xxnnN = i_xxnnN & i_0xxN
            n_0xxnnN = n_xxnnN - n_1xxnnN
            
            bias_xxnnN = 0
            for (g in unique(dat_xxN$G[i_xxnnN]))
            {
              if (g==0) {next}
              i_gxxN = dat_xxN$G==g
              
              # Compute empirical probabilities
              p = 0
              if (n_1xxnnN>0) {p = p+sum(i_1xxnnN & i_gxxN)/n_1xxnnN}
              if (n_0xxnnN>0) {p = p-sum(i_0xxnnN & i_gxxN)/n_0xxnnN}
              
              # Take g'=0
              bias_xxnnN = bias_xxnnN +
                (outcome_mean(0,g,x1,x2,gamma)-outcome_mean(0,0,x1,x2,gamma))*p
            }
            bias = bias + bias_xxnnN*n_xxnnN
          }
        }
      }
    }
  }
  return(bias/nrow(dat))
}

# Compute bias contribution for each level of covariates (Corollary 2, eqn.11)
bias_C2N_contr = function(dat, gamma)
{
  classes = data.frame(game1=integer(), game2=integer(),
                       Ngame1=double(), Ngame2=double(),
                       N=integer(),     G=double(),
                       nZ1=integer(),   nZ0=integer(),
                       bias=double())
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
        i_0xxN = !i_1xxN
        
        for (n1 in unique(dat_xxN$Ngame1))
        {
          for (n2 in unique(dat_xxN$Ngame2))
          {
            i_xxnnN = dat_xxN$Ngame1==n1 & dat_xxN$Ngame2==n2
            n_xxnnN = sum(i_xxnnN)
            if (n_xxnnN==0) {next}
            i_1xxnnN = i_xxnnN & i_1xxN
            n_1xxnnN = sum(i_1xxnnN)
            i_0xxnnN = i_xxnnN & i_0xxN
            n_0xxnnN = n_xxnnN - n_1xxnnN
            
            for (g in unique(dat_xxN$G[i_xxnnN]))
            {
              i_gxxN = dat_xxN$G==g
              i_1gxxnnN = i_1xxnnN & i_gxxN
              i_0gxxnnN = i_0xxnnN & i_gxxN
              
              # Compute empirical probabilities
              p = 0
              if (n_1xxnnN>0) {p = p+sum(i_1gxxnnN)/n_1xxnnN}
              if (n_0xxnnN>0) {p = p-sum(i_0gxxnnN)/n_0xxnnN}
              
              # Take g'=0
              classes = rbind(classes,
                              c(x1,x2,n1,n2,N,g,sum(i_1gxxnnN),sum(i_0gxxnnN),
                                (outcome_mean(0,g,x1,x2,gamma)-outcome_mean(0,0,x1,x2,gamma))*p*n_xxnnN/nrow(dat)))
            }
          }
        }
      }
    }
  }
  colnames(classes) = c("game1","game2","Ngame1","Ngame2","N","G","nZ1","nZ0","bias")
  return(classes)
}



# Subclassification functions
# ---------------------------------------------------------

# Classifies units into j roughly equal sized subclasses
# based on individual covariates
subcl_ind = function(dat, j)
{
  if (j > 1)
  {
    # Estimate propensity score via logistic regression
    # and subclassify according to propensity score
    classes = cut(glm(Z~game1+game2,data=dat,family=binomial)$fitted.values,
                  j, labels=F)
  } else {
    classes = rep(1, nrow(dat))
  }
  return(classes)
}

# Classifies units into j roughly equal sized subclasses
# based on all covariates
subcl_all = function(dat, j)
{
  if (j > 1)
  {
    # Estimate propensity score via logistic regression
    # and subclassify according to propensity score
    classes = cut(glm(Z~game1+game2+Ngame1+Ngame2+N,
                      data=dat,family=binomial)$fitted.values,
                  j, labels=F)
  } else {
    classes = rep(1, nrow(dat))
  }
  return(classes)
}

# Estimates treatment effect based on given subclasses
subcl_est = function(dat, classes)
{
  i_z = dat$Z==1
  tau_est = 0
  for (c in unique(classes))
  {
    # Compute contrast within each subclass
    i_c = classes==c
    tau_c = 0
    Y_c1 = mean(dat$Y[i_c & i_z])
    Y_c0 = mean(dat$Y[i_c & !i_z])
    if (!is.na(Y_c1)) {tau_c = tau_c + Y_c1}
    if (!is.na(Y_c0)) {tau_c = tau_c - Y_c0}
    tau_est = tau_est + tau_c*sum(i_c)
  }
  return(tau_est/nrow(dat))
}

# Estimates treatment effect based on individual and neighbourhood
# propensity scores (Forastiere, 2021)
subcl_gps = function(dat, j, sum_exposure=F)
{
  # Subclassify units based on individual propensity
  ind_classes = subcl_all(dat, j)
  
  n = nrow(dat)
  obs_g = sort(unique(dat$G))
  n_g = length(obs_g)
  tau_g = rep(0, n_g)
  
  # If sum neighbourhood treatments, compute size of V_g's
  if (sum_exposure)
  {
    v_g = rep(0, n_g)
    for (i in 1:n_g) {v_g[i] = sum(dat$N>=obs_g[i])}
  }
  
  for (c in unique(ind_classes))
  {
    dat_c = dat[ind_classes==c,]
    n_c = nrow(dat_c)
    
    if (sum_exposure)
    {
      # Fit a negative binomial model on neighbourhood treatment
      # Note: may throw error if there are too few units in subclass
      #       These subclasses typically correspond to units with
      #       large N. May be reasonable to assume any estimated
      #       propensity score will be very small, so just set
      #       the scores to 0.
      prop_fit = tryCatch(glm.nb(G~Z+Ngame1+Ngame2+N,data=dat_c),
                          error=function(e) {return(NULL)})
      dat_c$nbr_prop = predict_nb(prop_fit, dat_c)
      lin_fit = lm(Y~Z+G+nbr_prop, data=dat_c)
      
      for (i in 1:n_g)
      {
        g = obs_g[i]
        i_g = dat_c$N >= g
        n_cg = sum(i_g)
        if (n_cg==0) {break} # Remaining G values impossible for units in subclass
        dat_c$G = g
        
        # Estimate neighbourhood propensity scores and outcomes
        dat_c$Z = 1
        dat_c$nbr_prop[i_g] = predict_nb(prop_fit,dat_c[i_g,])
        Y_c1g = mean(predict(lin_fit,newdata=dat_c[i_g,]))
        dat_c$Z = 0
        dat_c$nbr_prop[i_g] = predict_nb(prop_fit,dat_c[i_g,])
        Y_c0g = mean(predict(lin_fit,newdata=dat_c[i_g,]))
        
        tau_g[i] = tau_g[i] + (Y_c1g-Y_c0g)*n_cg/v_g[i]
      }
    } else {
      # Fit a logistic model on neighbourhood treatment
      # Note: assumes V_g = all units in subclass
      prop_fit = glm(G~Z+Ngame1+Ngame2+N, data=dat_c, family=binomial, weights=N)
      dat_c$nbr_prop = prop_fit$fitted.values
      lin_fit = lm(Y~Z+G+nbr_prop, data=dat_c)
      
      for (i in 1:n_g)
      {
        g = obs_g[i]
        dat_c$G = g
        
        # Estimate neighbourhood propensity scores and outcomes
        dat_c$Z = 1
        dat_c$nbr_prop = predict(prop_fit, newdata=dat_c)
        Y_c1g = mean(predict(lin_fit,newdata=dat_c))
        dat_c$Z = 0
        dat_c$nbr_prop = predict(prop_fit, newdata=dat_c)
        Y_c0g = mean(predict(lin_fit,newdata=dat_c))
        
        tau_g[i] = tau_g[i] + (Y_c1g-Y_c0g)*n_c/n
      }
    }
  }
  return(sum(tau_g * table(dat$G)/n))
}
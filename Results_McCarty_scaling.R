################################################################################
################################################################################
##
## Purpose and functions: Implements the various degree ratio adjustments and
## scaling procedures on the McCarty dataset.
##
## Input Files: McCarty_LS.RData
##
## Figures made: (1) Figure 3a of manuscript: McCarty_meal_dist_class_subset.pdf
##               (2) Figure 4 of manuscript: McCarty_meal_results.pdf
##
## Files created: None
##
################################################################################
################################################################################

library(ggplot2)
library(reshape2)
library(Metrics)


load("./Results/McCarty_LS.RData")

# Preparation -------------------------------------------------------------

latent.1 = res.post$est.latent.pos.1
latent.2 = res.post$est.latent.pos.2
latent.3 = res.post$est.latent.pos.3

latent.group.1 = res.post$est.group.pos.1
latent.group.2 = res.post$est.group.pos.2
latent.group.3 = res.post$est.group.pos.3


## Perform thinning and burn-in
## Sampler also thinned by 2 by default
latent.1 = latent.1[-c(1:1001), ]
latent.2 = latent.2[-c(1:1001), ]
latent.3 = latent.3[-c(1:1001), ]
latent.group.1 = latent.group.1[-c(1:1001), ]
latent.group.2 = latent.group.2[-c(1:1001), ]
latent.group.3 = latent.group.3[-c(1:1001), ]

latent.1 = latent.1[seq(1, nrow(latent.1), by = 5), ]
latent.2 = latent.2[seq(1, nrow(latent.2), by = 5), ]
latent.3 = latent.3[seq(1, nrow(latent.3), by = 5), ]
latent.group.1 = latent.group.1[seq(1, nrow(latent.group.1), by = 5), ]
latent.group.2 = latent.group.2[seq(1, nrow(latent.group.2), by = 5), ]
latent.group.3 = latent.group.3[seq(1, nrow(latent.group.3), by = 5), ]


## Calculate distance to groups
dist.ls.mat = matrix(NA, nrow = N.i, ncol = N.k)
for (k in 1:N.k) {
  dist.to.htr = matrix(NA, nrow = nrow(latent.group.1) - 1, ncol = N.i)
  for (i in 2:nrow(latent.group.1)) {
    ## Find center
    dist.to.htr[i - 1,] = acos(
      latent.1[i,] * latent.group.1[i, k] +
        latent.2[i,] * latent.group.2[i, k] +
        latent.3[i,] * latent.group.3[i, k]
    )
  }
  dist.ls.mat[, k] = colMeans(dist.to.htr)
}

colnames(dist.ls.mat) = names(known)
dist.ls.df = data.frame(degree = degrees.known, dist.ls.mat)
names(dist.ls.df) = c("degree", names(known))
dist.ls.melt = reshape2::melt(dist.ls.df, id.vars = c("degree"))
dist.ls.melt$group = paste0(dist.ls.melt$variable)

y.known = y[known.ind]
known.known = known[known.ind]

killworth.est = killworth.error = rep(NA, n.known)
for (k in 1:n.known) {
  degrees.est = N * rowSums(y.known[, -k]) / sum(known.known[-k])
  killworth.est[k] = N * sum(y.known[, k]) / sum(degrees.est)
  killworth.error[k] = as.numeric(killworth.est[k] / known.known[k])
}

## Calculate killlworth relative error
killworth.rel.error = as.numeric((known.known - killworth.est) / known.known)
names(killworth.error) = names(known)[known.ind]

min.distances = apply(dist.ls.mat, 2, min)[known.ind]
max.distances = apply(dist.ls.mat, 2, max)[known.ind]
ranges.distances = max.distances - min.distances

known.known = known[known.ind]
dist.ls.mat.known = dist.ls.mat[, known.ind]

# Just Degree Ratio -------------------------------------------------------
## Perform only the degree ratio adjustment, but not the scaling

## For each subpopulations:
x.fact.vec = rep(NA, n.known)
adj.factor = killworth.dr = rep(NA, n.known)
resp.use = list()
for (outer.k in 1:n.known) {
  known.ind.sub = (1:n.known)[-outer.k]
  ## Step 1) Estimate the unknown subpopulation and degrees
  ## Estimate and stored in killworth.est
  degrees.known = N * rowSums(y.known[-outer.k]) / sum(known.known[-outer.k])
  
  ## Step 2) Get LOO estimates for others using degrees.known
  loo.sub = N * colSums(y.known[, known.ind.sub]) / sum(degrees.known)
  
  ## Step 3) Estimate degree ratio
  opt.func = function(x.fact) {
    resp.use = list()
    adj.factor = killworth.est.adj = rep(NA, length(known.ind.sub))
    for (k in 1:length(known.ind.sub)) {
      resp.use[[k]] = which(dist.ls.mat.known[, known.ind.sub[k]] <
                              (min.distances[known.ind.sub[k]] +
                                 x.fact * ranges.distances[known.ind.sub[k]]))
      
      adj.factor[k] = mean(degrees.known[resp.use[[k]]]) / mean(degrees.known)
      ## Adjust current k
      killworth.est.adj[k] = loo.sub[k] / adj.factor[k]
    }
    
    adj.mape = Metrics::mape(known.known[known.ind.sub], killworth.est.adj)
    return(adj.mape)
  }
  
  opt.func(0.2)
  
  optim.fit = optimize(opt.func,
                       interval = c(0, 1),
                       maximum = FALSE)
  x.fact.vec[outer.k] = optim.fit$minimum
  
  resp.use[[outer.k]] = which(dist.ls.mat.known[, outer.k] <
                                (min.distances[outer.k] +
                                   x.fact.vec[outer.k] * ranges.distances[outer.k]))
  
  ## Adjust estimate
  adj.factor[outer.k] = mean(degrees.known[resp.use[[outer.k]]]) / mean(degrees.known)
  killworth.dr[outer.k] = killworth.est[outer.k] / adj.factor[outer.k]
}






# Just Scaling ------------------------------------------------------------
## Perform just the scaling, but not the degree ratio adjustment

## For the hidden population, use all other subpopulations to estimate C

## Step 0) Estimate the unknown subpopulation
## Done and stored in killworth.est

killworth.scale = rep(NA, n.known)
for (outer.k in 1:n.known) {
  ## Step 1) Perform LOO for all other subpopulations using only known subpopulation
  known.ind.sub = (1:n.known)[-outer.k]
  degrees.known = N * rowSums(y.known[-outer.k]) / sum(known.known[-outer.k])
  
  ## Step 2) Get LOO estimates for others using degrees.known
  loo.sub = N * colSums(y.known[, known.ind.sub]) / sum(degrees.known)
  
  ## Step 3) Calculate C
  C = mean(loo.sub / (known.known[-outer.k]))
  
  ## Step 4) Scale the unknown size estimate
  killworth.scale[outer.k] = killworth.est[outer.k] / C
}


# Degree Ratio then Scaling -----------------------------------------------
## Perform the degree ratio adjustment first, and then scale

x.fact.vec = rep(NA, n.known)
dr.scale.factor = killworth.dr.scale = rep(NA, n.known)
for (outer.k in 1:n.known) {
  all.tmp = rep(NA, n.known)
  ## Step 1) Perform LOO for all other subpopulations using only known subpopulation
  known.ind.sub = (1:n.known)[-outer.k]
  degrees.known = N * rowSums(y.known[-outer.k]) / sum(known.known[-outer.k])
  
  ## Step 2) Get LOO estimates for others using degrees.known
  loo.sub = N * colSums(y.known[, known.ind.sub]) / sum(degrees.known)
  
  ## Step 3) Estimate degree ratio
  ## Adjust unknown estimate
  opt.func = function(x.fact, known.ind.sub.in) {
    resp.use = list()
    adj.factor.tmp = killworth.est.adj = rep(NA, length(known.ind.sub.in))
    for (k in 1:length(known.ind.sub.in)) {
      resp.use[[k]] = which(dist.ls.mat.known[, known.ind.sub.in[k]] <
                              (min.distances[known.ind.sub.in[k]] +
                                 x.fact * ranges.distances[known.ind.sub.in[k]]))
      
      adj.factor.tmp[k] = mean(degrees.known[resp.use[[k]]]) / mean(degrees.known)
      ## Adjust current k
      killworth.est.adj[k] = loo.sub[k] / adj.factor.tmp[k]
    }
    
    adj.mape = Metrics::mape(known.known[known.ind.sub.in], killworth.est.adj)
    return(adj.mape)
  }
  
  opt.func(0.2, known.ind.sub)
  
  optim.fit = optimize(
    opt.func,
    interval = c(0, 1),
    known.ind.sub.in = known.ind.sub,
    maximum = FALSE
  )
  x.fact.vec[outer.k] = optim.fit$minimum
  resp.use = which(dist.ls.mat.known[, outer.k] <
                     (min.distances[outer.k] +
                        x.fact.vec[outer.k] * ranges.distances[outer.k]))
  adj.factor[outer.k] = mean(degrees.known[resp.use]) / mean(degrees.known)
  all.tmp[outer.k] = killworth.est[outer.k] / adj.factor[outer.k]
  
  ## Adjust known estimates using remaining subpopulations
  for (k in known.ind.sub) {
    k.ind = which(known.ind.sub == k)
    known.ind.sub.inner = known.ind.sub[-k.ind]
    
    opt.func = function(x.fact, known.ind.sub.in) {
      resp.use = list()
      adj.factor.tmp = killworth.est.adj = rep(NA, length(known.ind.sub.in))
      for (k in 1:length(known.ind.sub.in)) {
        resp.use[[k]] = which(dist.ls.mat.known[, known.ind.sub.in[k]] <
                                (min.distances[known.ind.sub.in[k]] +
                                   x.fact * ranges.distances[known.ind.sub.in[k]]))
        
        adj.factor.tmp[k] = mean(degrees.known[resp.use[[k]]]) / mean(degrees.known)
        ## Adjust current k
        killworth.est.adj[k] = loo.sub[k] / adj.factor.tmp[k]
      }
      
      adj.mape = Metrics::mape(known.known[known.ind.sub.in], killworth.est.adj)
      return(adj.mape)
    }
    
    opt.func(0.2, known.ind.sub.inner)
    
    optim.fit = optimize(
      opt.func,
      interval = c(0, 1),
      known.ind.sub.in = known.ind.sub.inner,
      maximum = FALSE
    )
    x.fact.vec.tmp = optim.fit$minimum
    resp.use = which(dist.ls.mat.known[, k] <
                       (min.distances[k] + x.fact.vec.tmp * ranges.distances[k]))
    adj.factor.tmp = mean(degrees.known[resp.use]) / mean(degrees.known)
    all.tmp[k] = loo.sub[k.ind] / adj.factor.tmp
  }
  
  
  
  ## Step 3) Calculate C
  C = mean(all.tmp[-outer.k] / (known.known[-outer.k]))
  
  ## Step 4) Scale all size estimates
  all.tmp = all.tmp / C
  
  killworth.dr.scale[outer.k] = all.tmp[outer.k]
  
  
  
  print(outer.k)
}


# Scaling then Degree Ratio -----------------------------------------------
## Perform the scaling first, and then adjust for the degree ratio


x.fact.vec = rep(NA, n.known)
scale.dr.factor = killworth.scale.dr = rep(NA, n.known)
for (outer.k in 1:n.known) {
  ## Step 1) Perform LOO for all other subpopulations using only known subpopulation
  known.ind.sub = (1:n.known)[-outer.k]
  degrees.known = N * rowSums(y.known[-outer.k]) / sum(known.known[-outer.k])
  
  ## Step 2) Get LOO estimates for others using degrees.known
  loo.sub = N * colSums(y.known[, known.ind.sub]) / sum(degrees.known)
  
  ## Step 3) Calculate C
  C = mean(loo.sub / (known.known[-outer.k]))
  
  ## Step 4) Scale all size estimates
  all.tmp = rep(NA, n.known)
  all.tmp[outer.k] = killworth.est[outer.k]
  all.tmp[known.ind.sub] = loo.sub
  killworth.scale.tmp = all.tmp / C
  degrees.known = degrees.known * C
  
  loo.sub.scale = killworth.scale.tmp[-outer.k]
  
  ## Step 3) Estimate degree ratio
  opt.func = function(x.fact) {
    resp.use = list()
    adj.factor.tmp = killworth.est.adj = rep(NA, length(known.ind.sub))
    for (k in 1:length(known.ind.sub)) {
      resp.use[[k]] = which(dist.ls.mat.known[, known.ind.sub[k]] <
                              (min.distances[known.ind.sub[k]] +
                                 x.fact * ranges.distances[known.ind.sub[k]]))
      
      adj.factor.tmp[k] = mean(degrees.known[resp.use[[k]]]) / mean(degrees.known)
      ## Adjust current k
      killworth.est.adj[k] = loo.sub.scale[k] / adj.factor.tmp[k]
    }
    
    adj.mape = Metrics::mape(known.known[known.ind.sub], killworth.est.adj)
    return(adj.mape)
  }
  
  opt.func(0.2)
  
  optim.fit = optimize(opt.func,
                       interval = c(0, 1),
                       maximum = FALSE)
  x.fact.vec[outer.k] = optim.fit$minimum
  
  resp.use = list()
  resp.use[[outer.k]] = which(dist.ls.mat.known[, outer.k] <
                                (min.distances[outer.k] +
                                   x.fact.vec[outer.k] * ranges.distances[outer.k]))
  
  ## Adjust estimate
  adj.factor[outer.k] = mean(degrees.known[resp.use[[outer.k]]]) / mean(degrees.known)
  killworth.scale.dr[outer.k] = killworth.scale.tmp[outer.k] / adj.factor[outer.k]
  
  print(outer.k)
}


# Numerical Results -------------------------------------------------------


## DR to original
(
  Metrics::mape(known.known, killworth.est) -
    Metrics::mape(known.known, killworth.dr)
) /
  Metrics::mape(known.known, killworth.est) * 100

## Scaled to original
(
  Metrics::mape(known.known, killworth.est) -
    Metrics::mape(known.known, killworth.scale)
) /
  Metrics::mape(known.known, killworth.est) * 100

## DR.scaled
(
  Metrics::mape(known.known, killworth.est) -
    Metrics::mape(known.known, killworth.dr.scale)
) /
  Metrics::mape(known.known, killworth.est) * 100

## Scaled.DR
(
  Metrics::mape(known.known, killworth.est) -
    Metrics::mape(known.known, killworth.scale.dr)
) /
  Metrics::mape(known.known, killworth.est) * 100





















# Repeat without names ----------------------------------------------------
## Repeat the entire analysis, but only adjusting and evaluating the non-name
## subpopulations



y.orig = y.known
known.orig = known.known

## Trim values

non.names.ind = c(13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 29, 30, 31, 32)
non.names.ind.sub = 13:29

known.known = known[non.names.ind]
known.ind = non.names.ind
dist.ls.mat.known = dist.ls.mat[, non.names.ind]
n.known = length(non.names.ind)
killworth.est = killworth.est[non.names.ind.sub]
y = y.known = y[, non.names.ind]


min.distances = apply(dist.ls.mat, 2, min)[known.ind]
max.distances = apply(dist.ls.mat, 2, max)[known.ind]
ranges.distances = max.distances - min.distances

known.known = known[known.ind]
dist.ls.mat.known = dist.ls.mat[, known.ind]

# Just Degree Ratio -------------------------------------------------------


## For each subpopulations:

x.fact.vec = rep(NA, n.known)
adj.factor = killworth.dr = rep(NA, n.known)
resp.use = list()
for (outer.k in 1:n.known) {
  known.ind.sub = (1:n.known)[-outer.k]
  ## Step 1) Estimate the unknown subpopulation and degrees
  ## Estimate and stored in killworth.est
  degrees.known = N * rowSums(y.orig[-non.names.ind[outer.k]]) /
    sum(known.orig[-non.names.ind[outer.k]])
  
  ## Step 2) Get LOO estimates for others using degrees.known
  loo.sub = N * colSums(y.known[, known.ind.sub]) / sum(degrees.known)
  
  ## Step 3) Estimate degree ratio
  opt.func = function(x.fact) {
    resp.use = list()
    adj.factor = killworth.est.adj = rep(NA, length(known.ind.sub))
    for (k in 1:length(known.ind.sub)) {
      resp.use[[k]] = which(dist.ls.mat.known[, known.ind.sub[k]] <
                              (min.distances[known.ind.sub[k]] +
                                 x.fact * ranges.distances[known.ind.sub[k]]))
      
      adj.factor[k] = mean(degrees.known[resp.use[[k]]]) / mean(degrees.known)
      ## Adjust current k
      killworth.est.adj[k] = loo.sub[k] / adj.factor[k]
    }
    
    adj.mape = Metrics::mape(known.known[known.ind.sub], killworth.est.adj)
    return(adj.mape)
  }
  
  opt.func(0.2)
  
  optim.fit = optimize(opt.func,
                       interval = c(0, 1),
                       maximum = FALSE)
  x.fact.vec[outer.k] = optim.fit$minimum
  
  resp.use[[outer.k]] = which(dist.ls.mat.known[, outer.k] <
                                (min.distances[outer.k] +
                                   x.fact.vec[outer.k] * ranges.distances[outer.k]))
  
  ## Adjust estimate
  adj.factor[outer.k] = mean(degrees.known[resp.use[[outer.k]]]) / mean(degrees.known)
  killworth.dr[outer.k] = killworth.est[outer.k] / adj.factor[outer.k]
}






# Just Scaling ------------------------------------------------------------
## For the hidden population, use all other subpopulations to estimate C

## Step 0) Estimate the unknown subpopulation
## Done and stored in killworth.est

killworth.scale = rep(NA, n.known)
for (outer.k in 1:n.known) {
  ## Step 1) Perform LOO for all other subpopulations using only known subpopulation
  known.ind.sub = (1:n.known)[-outer.k]
  degrees.known = N * rowSums(y.orig[-non.names.ind[outer.k]]) /
    sum(known.orig[-non.names.ind[outer.k]])
  
  ## Step 2) Get LOO estimates for others using degrees.known
  loo.sub = N * colSums(y.known[, known.ind.sub]) / sum(degrees.known)
  
  ## Step 3) Calculate C
  C = mean(loo.sub / (known.known[-outer.k]))
  
  ## Step 4) Scale the unknown size estimate
  killworth.scale[outer.k] = killworth.est[outer.k] / C
}


# Degree Ratio then Scaling -----------------------------------------------

x.fact.vec = rep(NA, n.known)
dr.scale.factor = killworth.dr.scale = rep(NA, n.known)
for (outer.k in 1:n.known) {
  all.tmp = rep(NA, n.known)
  ## Step 1) Perform LOO for all other subpopulations using only known subpopulation
  known.ind.sub = (1:n.known)[-outer.k]
  degrees.known = N * rowSums(y.orig[-non.names.ind[outer.k]]) /
    sum(known.orig[-non.names.ind[outer.k]])
  
  ## Step 2) Get LOO estimates for others using degrees.known
  loo.sub = N * colSums(y.known[, known.ind.sub]) / sum(degrees.known)
  
  ## Step 3) Estimate degree ratio
  ## Adjust unknown estimate
  opt.func = function(x.fact, known.ind.sub.in) {
    resp.use = list()
    adj.factor.tmp = killworth.est.adj = rep(NA, length(known.ind.sub.in))
    for (k in 1:length(known.ind.sub.in)) {
      resp.use[[k]] = which(dist.ls.mat.known[, known.ind.sub.in[k]] <
                              (min.distances[known.ind.sub.in[k]] +
                              x.fact * ranges.distances[known.ind.sub.in[k]]))
      
      adj.factor.tmp[k] = mean(degrees.known[resp.use[[k]]]) / mean(degrees.known)
      ## Adjust current k
      killworth.est.adj[k] = loo.sub[k] / adj.factor.tmp[k]
    }
    
    adj.mape = Metrics::mape(known.known[known.ind.sub.in], killworth.est.adj)
    return(adj.mape)
  }
  
  opt.func(0.2, known.ind.sub)
  
  optim.fit = optimize(
    opt.func,
    interval = c(0, 1),
    known.ind.sub.in = known.ind.sub,
    maximum = FALSE
  )
  x.fact.vec[outer.k] = optim.fit$minimum
  resp.use = which(dist.ls.mat.known[, outer.k] <
                     (min.distances[outer.k] +
                        x.fact.vec[outer.k] * ranges.distances[outer.k]))
  adj.factor[outer.k] = mean(degrees.known[resp.use]) / mean(degrees.known)
  all.tmp[outer.k] = killworth.est[outer.k] / adj.factor[outer.k]
  
  ## Adjust known estimates using remaining subpopulations
  for (k in known.ind.sub) {
    k.ind = which(known.ind.sub == k)
    known.ind.sub.inner = known.ind.sub[-k.ind]
    
    opt.func = function(x.fact, known.ind.sub.in) {
      resp.use = list()
      adj.factor.tmp = killworth.est.adj = rep(NA, length(known.ind.sub.in))
      for (k in 1:length(known.ind.sub.in)) {
        resp.use[[k]] = which(dist.ls.mat.known[, known.ind.sub.in[k]] <
                                (min.distances[known.ind.sub.in[k]] +
                                x.fact * ranges.distances[known.ind.sub.in[k]]))
        
        adj.factor.tmp[k] = mean(degrees.known[resp.use[[k]]]) / mean(degrees.known)
        ## Adjust current k
        killworth.est.adj[k] = loo.sub[k] / adj.factor.tmp[k]
      }
      
      adj.mape = Metrics::mape(known.known[known.ind.sub.in], killworth.est.adj)
      return(adj.mape)
    }
    
    opt.func(0.2, known.ind.sub.inner)
    
    optim.fit = optimize(
      opt.func,
      interval = c(0, 1),
      known.ind.sub.in = known.ind.sub.inner,
      maximum = FALSE
    )
    x.fact.vec.tmp = optim.fit$minimum
    resp.use = which(dist.ls.mat.known[, k] <
                       (min.distances[k] + x.fact.vec.tmp * ranges.distances[k]))
    adj.factor.tmp = mean(degrees.known[resp.use]) / mean(degrees.known)
    all.tmp[k] = loo.sub[k.ind] / adj.factor.tmp
  }
  
  
  
  ## Step 3) Calculate C
  C = mean(all.tmp[-outer.k] / (known.known[-outer.k]))
  
  ## Step 4) Scale all size estimates
  all.tmp = all.tmp / C
  
  killworth.dr.scale[outer.k] = all.tmp[outer.k]
  
  
  
  print(outer.k)
}


# Scaling then Degree Ratio -----------------------------------------------

x.fact.vec = rep(NA, n.known)
scale.dr.factor = killworth.scale.dr = rep(NA, n.known)
resp.use = list()
for (outer.k in 1:n.known) {
  ## Step 1) Perform LOO for all other subpopulations using only known subpopulation
  known.ind.sub = (1:n.known)[-outer.k]
  degrees.known = N * rowSums(y.orig[-non.names.ind[outer.k]]) /
    sum(known.orig[-non.names.ind[outer.k]])
  
  ## Step 2) Get LOO estimates for others using degrees.known
  loo.sub = N * colSums(y.known[, known.ind.sub]) / sum(degrees.known)
  
  ## Step 3) Calculate C
  C = mean(loo.sub / (known.known[-outer.k]))
  
  ## Step 4) Scale all size estimates
  all.tmp = rep(NA, n.known)
  all.tmp[outer.k] = killworth.est[outer.k]
  all.tmp[known.ind.sub] = loo.sub
  killworth.scale.tmp = all.tmp / C
  degrees.known = degrees.known * C
  
  loo.sub.scale = killworth.scale.tmp[-outer.k]
  
  ## Step 3) Estimate degree ratio
  opt.func = function(x.fact) {
    resp.use = list()
    adj.factor.tmp = killworth.est.adj = rep(NA, length(known.ind.sub))
    for (k in 1:length(known.ind.sub)) {
      resp.use[[k]] = which(dist.ls.mat.known[, known.ind.sub[k]] <
                              (min.distances[known.ind.sub[k]] +
                                 x.fact * ranges.distances[known.ind.sub[k]]))
      
      adj.factor.tmp[k] = mean(degrees.known[resp.use[[k]]]) / mean(degrees.known)
      ## Adjust current k
      killworth.est.adj[k] = loo.sub.scale[k] / adj.factor.tmp[k]
    }
    
    adj.mape = Metrics::mape(known.known[known.ind.sub], killworth.est.adj)
    return(adj.mape)
  }
  
  opt.func(0.2)
  
  optim.fit = optimize(opt.func,
                       interval = c(0, 1),
                       maximum = FALSE)
  x.fact.vec[outer.k] = optim.fit$minimum
  
  
  resp.use[[outer.k]] = which(dist.ls.mat.known[, outer.k] <
                                (min.distances[outer.k] +
                                   x.fact.vec[outer.k] * ranges.distances[outer.k]))
  
  ## Adjust estimate
  scale.dr.factor[outer.k] = mean(degrees.known[resp.use[[outer.k]]]) /
    mean(degrees.known)
  killworth.scale.dr[outer.k] = killworth.scale.tmp[outer.k] /
    scale.dr.factor[outer.k]
  
  print(outer.k)
}





# Numerical Results -------------------------------------------------------


## DR to original
(
  Metrics::mape(known.known, killworth.est) -
    Metrics::mape(known.known, killworth.dr)
) /
  Metrics::mape(known.known, killworth.est) * 100

## Scaled to original
(
  Metrics::mape(known.known, killworth.est) -
    Metrics::mape(known.known, killworth.scale)
) /
  Metrics::mape(known.known, killworth.est) * 100

## DR.scaled
(
  Metrics::mape(known.known, killworth.est) -
    Metrics::mape(known.known, killworth.dr.scale)
) /
  Metrics::mape(known.known, killworth.est) * 100

## Scaled.DR
(
  Metrics::mape(known.known, killworth.est) -
    Metrics::mape(known.known, killworth.scale.dr)
) /
  Metrics::mape(known.known, killworth.est) * 100


# Plotting ----------------------------------------------------------------



degrees.known = N * rowSums(y.orig) / sum(known.orig)

class.mat = matrix("No", nrow = N.i, ncol = n.known)
for (k in 1:n.known) {
  class.mat[resp.use[[k]], k] = "Yes"
}


class.vec = c(class.mat)
dist.ls.df = data.frame(degree = degrees.known, dist.ls.mat.known)
dist.ls.melt.class = reshape2::melt(dist.ls.df, id.vars = c("degree"))
dist.ls.melt.class$class.member = class.vec
dist.ls.melt.class$class.member = factor(dist.ls.melt.class$class.member)
dist.ls.melt.class$class.member = factor(dist.ls.melt.class$class.member,
                                         levels = c("Yes", "No"))


df.class = subset(
  dist.ls.melt.class,
  variable %in%
    c("CAR.ACCIDENT", "DIABETES", "NATIVE.AMERICAN", "AIDS")
)
df.class$variable = as.character(df.class$variable)
df.class$variable[df.class$variable == "NATIVE.AMERICAN"] = "NATIVE AMERICAN"
df.class$variable[df.class$variable == "CAR.ACCIDENT"] = "CAR ACCIDENT"
df.class$variable = factor(df.class$variable)
gg.class = ggplot(df.class) +
  geom_point(aes(x = degree, y = value, col = class.member), size = 2) +
  facet_wrap( ~ variable, scales = "free") +
  xlab("Estimated Degree") + ylab("Estimated Distance") +
  theme_grey(base_size = 22) +
  labs(color = "Close") +
  theme(
    strip.text.x = element_text(size = 28, face = "bold"),
    axis.text.y = element_text(size = 28, face = "bold"),
    axis.title.x = element_text(size = 28, face = "bold"),
    axis.title.y = element_text(size = 28, face = "bold"),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 28,
      face = "bold"
    ),
    legend.text = element_text(size = 28, face = "bold"),
    legend.title = element_text(size = 28, face = "bold")
  )
gg.class

ggsave(
  "./Figures/McCarty_dist_class_subset.pdf",
  gg.class,
  width = 20,
  height = 10,
  units = "in"
)






rel.error.killworth = (known.known - killworth.est) / known.known
rel.error.threshold = (known.known - killworth.scale.dr) / known.known

df.rel.error = data.frame(
  name = names(known)[non.names.ind],
  true = known.known,
  killworth = rel.error.killworth,
  threshold = rel.error.threshold
)
df.rel.error$name = toupper(df.rel.error$name)
df.rel.error$better = NA
df.rel.error$better[which(abs(df.rel.error$threshold) <
                            abs(df.rel.error$killworth))] = "Yes"

df.melt.rel.error = reshape2::melt(df.rel.error,
                                   id.vars = c("name", "true", "better"))
df.melt.rel.error = subset(df.melt.rel.error,!is.na(true))
df.melt.rel.error$name = factor(df.melt.rel.error$name,
                                levels = names(known.known)[order(known.known)])


gg.res = ggplot(df.melt.rel.error,
                aes(
                  x = name,
                  y = value,
                  col = variable,
                  group = variable,
                  label = name
                )) +
  geom_path(
    aes(
      x = name,
      y = value,
      group = name,
      col = better
    ),
    arrow = arrow(length = unit(0.55, "cm")),
    lwd = 2,
    show.legend = FALSE
  ) +
  geom_point(size = 6) +
  scale_color_discrete(
    na.value = "transparent",
    breaks = c("killworth", "threshold"),
    labels = c("\nMLE\n", "\nScaled and \nAdjusted MLE\n")
  ) +
  guides(color = guide_legend(override.aes = list(shape = 19))) +
  geom_abline(intercept = 0,
              slope = 0,
              linewidth = 1) +
  xlab("Subpopulation") + ylab("Relative Error") +
  theme_grey(base_size = 28) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.y = element_text(size = 28, face = "bold"),
    axis.title.x = element_text(size = 28, face = "bold"),
    axis.title.y = element_text(size = 28, face = "bold"),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 28,
      face = "bold"
    ),
    legend.text = element_text(size = 28, face = "bold"),
    legend.title = element_text(size = 28, face = "bold")
  ) +
  labs(color = "Model")
gg.res

ggsave(
  "./Figures/McCarty_results.pdf",
  gg.res,
  width = 20,
  height = 10,
  units = "in"
)

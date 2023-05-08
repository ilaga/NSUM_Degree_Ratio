################################################################################
################################################################################
##
## Purpose and functions: Simulated data from a stochastic block model and
## fits the Latent Surface Model
##
## Input Files: Tyler_latent_surface_code_MH.R, main.R, 
##
## Figures made: None
##
## Files created: (1) SBM_LS.RData: Contains posterior samples
## from model
##
################################################################################
################################################################################

library(igraph)
library(movMF)

## Source latent surface code
source('./Latent_surface_code/Tyler_latent_surface_code_MH.R')
source('./Latent_surface_code/main.R')


set.seed(52)

## Initialize Parameters
N.i = 800
N.k = 8
means = diag(seq(0, 0.35, length.out = N.k))
means = means + 0.05

## Simulate SBM data
sbm.sim = sample_sbm(N.i, means, rep(100, N.k), directed = F, loops = F)
node.groups = rep(1:N.k, each = 100)
node.groups = model.matrix(~factor(node.groups) - 1)
colnames(node.groups) = 1:N.k
adj.mat = as.matrix(get.adjacency(sbm.sim))

## Calculate subpopulation sizes
known = known.sizes = colSums(node.groups)
known / N.i

which.group = apply(node.groups, 1, function(x){which(x == 1)})

known.degree = rowSums(adj.mat)

## Get ARD
ard = adj.mat %*% node.groups

N.k = ncol(ard)
N.i = nrow(ard)
N = nrow(ard)

killworth.degree.est = N * rowSums(ard) / sum(known)



# Fit latent surface ------------------------------------------------------

m.fix = 3
p.dim = 3
N = N.i
n.group = N.k
muk.fix.ind = c(2, 5, 7)
muk.fix = rbind(c(0, -1, 0),
                c(sqrt(2/3), sqrt(1/3), 0),
                c(-sqrt(2/3), sqrt(1/3), 0))


known.ind = 1:N.k
known.sizes = known


n.group.known = 1:N.k
res.post = main.posterior(
  y = ard,
  total.prop = sum(known[n.group.known] / N),
  muk.fix = muk.fix,
  n.iter = 30000,
  m.iter = 1,
  n.thin = 1,
  ls.dim = p.dim,
  is.sample = F
)

## Save Results
save.image("./Results/SBM_LS.RData")
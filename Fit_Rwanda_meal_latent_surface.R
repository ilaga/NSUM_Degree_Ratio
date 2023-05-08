################################################################################
################################################################################
##
## Purpose and functions: Fits the Latent Surface Model to the Rwanda Meal Data
##
## Input Files: (1) Latent_surface_code/Tyler_latent_surface_code_MH.R
##              (2) Latent_surface_code/main.R
##              (3) Data/rwanda-varnames.csv
##              (4) Data/qvsq_meal.RDS
##
## Figures made: None
##
## Files created: (1) Rwanda_meal_latent_fit.RData: Contains posterior samples
## from model
##
################################################################################
################################################################################

library(movMF)

## Source latent surface code
source('./Latent_surface_code/Tyler_latent_surface_code_MH.R')
source('./Latent_surface_code/main.R')



## Read in data
knownpop = read.csv("./Quality_vs_quantity_code/data/rwanda-varnames.csv")
nsum = readRDS("./Data/qvsq_meal.RDS")

subpop.names = knownpop$varname


y = nsum
names(y) = subpop.names


subpop.names = subpop.names[!is.na(knownpop$known.size)]
y = y[,which(!is.na(knownpop$known.size))]

y = y[complete.cases(y),]
N.i = nrow(y)
N.k = ncol(y)

N = 5e6
known = knownpop$known.size
known = known[!is.na(known)]
n.known = ncol(y)
known.ind = 1:n.known
known.sizes = known

degrees.known = N * rowSums(y) / sum(known)


# Fit latent surface ------------------------------------------------------

p.dim = 3
m.fix = 6
muk.fix.ind = sample(1:N.k, size = m.fix, replace = F)
muk.fix = matrix(rnorm(m.fix * p.dim), nrow = m.fix, ncol = p.dim)
muk.fix = sweep(muk.fix, MARGIN = 1, 1 / sqrt(rowSums(muk.fix ^ 2)), `*`)

n.group.known = 1:N.k
n.group = N.k
res.post = main.posterior(
  y = y,
  total.prop = sum(known / N, na.rm = T),
  muk.fix = muk.fix,
  n.iter = 30000,
  m.iter = 1,
  n.thin = 1,
  ls.dim = p.dim,
  is.sample = F
)

## Save Results
save.image("./Results/Rwanda_meal_latent_fit.RData")
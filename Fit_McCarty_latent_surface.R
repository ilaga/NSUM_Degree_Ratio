################################################################################
################################################################################
##
## Purpose and functions: Fits the Latent Surface Model to the McCarty Data
##
## Input Files: (1) Latent_surface_code/Tyler_latent_surface_code_MH.R
##              (2) Latent_surface_code/main.R
##              (3) twomwth2.csv
##
## Figures made: None
##
## Files created: (1) McCarty_LS.RData: Contains posterior samples from model
##
################################################################################
################################################################################


library(movMF)

## Source latent surface code
source('./Latent_surface_code/Tyler_latent_surface_code_MH.R')
source('./Latent_surface_code/main.R')

## Read in data
dat = read.csv("./Data/twomwth2.csv", header = T)

y = dat[,c(5:37)]
y = y[,-2]

known = data.frame(subn1a     = 3187000,
                   subn2a     = 351000,
                   subn3a     = 1255000,
                   subn4a     = 291000,
                   subn5a     = 4023000,
                   subn6a     = 1118000,
                   subn7a     = 874000,
                   subn8a     = 642000,
                   subn9a     = 3811000,
                   subn10a    = 510000,
                   subn11a    = 2865000,
                   subn12a    = 358000,
                   subtr1     = 2000000,
                   subtr2     = 4000000,
                   subtr3     = 118000,
                   subtr4     = 3335000,
                   subtr5     = 170000,
                   subtr6     = 810000,
                   subtr7     = 534000,
                   subtr8     = 190000,
                   subtr9     = 6500000,
                   subtr10    = 630000,
                   subtr11    = 5300000,
                   subtr12    = 240000,
                   subtr13 = NA,
                   subtr14    = 207000,
                   subtr15 = NA,
                   subtr16 = NA,
                   subtr17    = 892000,
                   subtr18    = 27400,
                   subtr19    = 30200,
                   subtr20    = 45200)
known = as.numeric(known)

names(known) = c("MICHAEL",
                 "CHRISTINA",
                 "CHRISTOPHER",
                 "JAQUELINE",
                 "JAMES",
                 "JENNIFER",
                 "ANTHONY",
                 "KIMBERLY",
                 "ROBERT",
                 "STEPHANIE",
                 "DAVID",
                 "NICOLE",
                 "NATIVE AMERICAN",
                 "GAVE BIRTH",
                 "ADOPTED",
                 "OVER 65",
                 "DIALYSIS",
                 "POSTAL WORKER",
                 "PILOT",
                 "JAYCEE",
                 "DIABETES",
                 "OPEN BUSINESS",
                 "TWIN",
                 "GUN DEALER",
                 "HIV",
                 "AIDS",
                 "HOMELESS",
                 "RAPED",
                 "JAILED",
                 "MURDERED",
                 "SUICIDE",
                 "CAR ACCIDENT")


y = y[complete.cases(y),]
N.i = nrow(y)
N.k = ncol(y)
known.ind = which(!is.na(known))
unknown.ind = which(is.na(known))
n.known = length(known.ind)
n.unknown = N.k - n.known
N = 250000000


degrees.known = N * rowSums(y[,known.ind]) / sum(known[known.ind])
log.degrees.known = log(degrees.known)

# Fit latent surface ------------------------------------------------------

p.dim = 3
m.fix = 6
muk.fix.ind = sample(1:N.k, size = m.fix, replace = F)
names(known)[muk.fix.ind]
muk.fix = matrix(rnorm(m.fix * p.dim), nrow = m.fix, ncol = p.dim)
muk.fix = sweep(muk.fix, MARGIN = 1, 1 / sqrt(rowSums(muk.fix ^ 2)), `*`)

known.sizes = known

n.group.known = known.ind
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
save.image("./Results/McCarty_LS.RData")
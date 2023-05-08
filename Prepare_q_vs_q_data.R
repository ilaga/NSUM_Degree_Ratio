################################################################################
################################################################################
##
## Purpose and functions: Read and process raw Rwanda data into the prepared
##
## Input Files: (1) Latent_surface_code/Tyler_latent_surface_code_MH.R
##              (2) Latent_surface_code/main.R
##              (3) Quality_vs_quantity_code/data/rwanda-varnames.csv
##              (4) Quality_vs_quantity_data/RWIQ6ASD/RWIQ6AFL.SAS7BDAT
##
## Rwanda Meal ARD dataset for analysis
##
## Figures made: None
##
## Files created: (1) qvsq_meal_cov.RDS: Contains Rwanda meal ARD dataset
##                (2) qvsq_standard_cov.RDS: Contains Rwanda standard ARD
##                    dataset (not used in this)
##
################################################################################
################################################################################

library(haven)

## Read in data
known.mat = read.csv("./Data/rwanda-varnames.csv")
dat = read_sas("./Data/RWIQ6AFL.SAS7BDAT")

keep.ind = which(
  names(dat) %in% paste0(
    "Q",
    known.mat$qnum
  )
)

dat.ard = dat[,keep.ind]

dat.meal = dat.ard[which(dat$QHQTYPE == 2),]
dat.standard = dat.ard[which(dat$QHQTYPE == 1),]

saveRDS(dat.meal, file = "./Data/qvsq_meal.RDS")
saveRDS(dat.meal, file = "./Data/qvsq_standard.RDS")

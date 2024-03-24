library(ggplot2)
library(reshape2)
library(ggpmisc)

library(dplyr)






# McCarty -----------------------------------------------------------------

## Read in data
dat = read.csv("./Data/twomwth2.csv", header = T)

y = dat[,c(5:37)]
y = y[,-2]

sizes = data.frame(subn1a     = 3187000,
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
sizes = as.numeric(sizes)

names(sizes) = c("MICHAEL",
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
known.ind = which(!is.na(sizes))
unknown.ind = which(is.na(sizes))
n.known = length(known.ind)
n.unknown = N.k - n.known
N = 250000000




## Subset to known
y = y[,known.ind]
known = sizes[known.ind]

y = y[-which(rowSums(y) == 0),]

N.i = nrow(y)
N.k = ncol(y)


degrees.est.all = N * rowSums(y) / sum(known)


## Define weights
degrees.loo.est = matrix(NA, nrow = N.i, ncol = n.known)
ratio.mat = ratio.mat.2 = matrix(NA, nrow = N.i, ncol = n.known)
for(k in 1:n.known){
  degrees.loo.est[,k] = N * rowSums(y[,-k]) / sum(known[-k])
  ratio.mat[,k] = y[,k]
}



## Scale
# ratio.mat = scale(ratio.mat, center = FALSE)
ratio.mat = scale(ratio.mat, center = FALSE, scale = apply(ratio.mat, 2, mean, na.rm = TRUE))




## Plot ratio against estimated loo degree
df.deg = data.frame(degrees.loo.est)
names(df.deg) = names(known)
df.ratio = data.frame(ratio.mat)
names(df.ratio) = names(known)

degrees.melt = melt(df.deg)
logratio.melt = melt(df.ratio)
logratio.melt$degrees = degrees.melt$value
names(logratio.melt)[2] = "LogRatio"

## Check slopes
slopes = rep(NA, N.k)
for(k in 1:N.k){
  lm.fit.2 = lm(df.ratio[,k] ~ df.deg[,k])
  slopes[k] = coef(lm.fit.2)[2]
}
names(slopes) = names(known)

# PIMLE/MLE LOO Estimates -------------------------------------------------


mle.loo.est = rep(NA, n.known)
for (k in 1:n.known) {
  known.ind.sub = (1:n.known)[-k]
  ## Step 1) Estimate the unknown subpopulation and degrees
  degrees.est = N * rowSums(y[,-k]) / sum(known[-k])
  
  ## Step 2) Get LOO estimates
  mle.loo.est[k] = N * sum(unlist(y[, k])) / sum(degrees.est)
}





# Calculate relative error ------------------------------------------------
mle.rel.err = known / mle.loo.est ## Either this or full cox transformation.
# mle.rel.err = mle.loo.est / known


# Plot against intercepts and slopes --------------------------------------

slope.df = data.frame(slope = slopes, name = names(known),
                      mle.err = mle.rel.err)

slope.melt = melt(slope.df, measure.vars = c("slope"))

gg.slope = ggplot(slope.melt, aes(x = value, y = mle.err, label = name)) +
  geom_text(size = 7, fontface = "bold") +
  ylab(expression(N[k] / hat(N)[k])) + xlab(expression(hat(beta)[1])) +
  theme_grey(base_size = 22) +
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
  ) +
  xlim(0.0008, 0.006)

gg.slope


ggsave(
  "./Figures/McCarty_slopes.pdf",
  gg.slope,
  width = 20,
  height = 10,
  units = "in"
)









gg.ratio = ggplot(subset(logratio.melt, variable %in%
                           c("MICHAEL", "GAVE BIRTH",
                             "CAR ACCIDENT", "SUICIDE")),
                  aes(y = LogRatio, x = degrees)) + geom_point() +
  facet_wrap(~variable, scales = "free_x") +
  stat_poly_line(se = FALSE, linewidth = 3) +
  stat_poly_eq(use_label(c("eq")), size = 8) +
  geom_point(size = 2) +
  xlab(expression(hat(d[i]))) + ylab(expression(y[i,k])) +
  theme_grey(base_size = 22) +
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

gg.ratio

ggsave(
  "./Figures/McCarty_degree_reg.pdf",
  gg.ratio,
  width = 20,
  height = 10,
  units = "in"
)







# Perform scaling using slope ---------------------------------------------

## Using slope.df
killworth.slope.scale = mle.err.pred = rep(NA, n.known)
for(k in 1:n.known){
  slope.sub = slope.df[-k,]
  
  # Raw error ---------------------------------------------------------------
  
  ## Fit linear regression
  lm.fit.sub = lm(mle.err ~ slope, data = slope.sub)
  
  ## Predict log relative error
  logmle.err.pred = predict(lm.fit.sub, newdata = data.frame(slope = slope.df$slope[k]))
  
  ## Predict relative error
  mle.err.pred[k] = logmle.err.pred
  
  ## Correct for killworth estimate
  killworth.slope.scale[k] = mle.loo.est[k] * mle.err.pred[k] ## Multiply because inverse error
}


# Numerical Results -------------------------------------------------------


## Slope.scale
(
  Metrics::mape(known, mle.loo.est) -
    Metrics::mape(known, killworth.slope.scale)
) /
  Metrics::mape(known, mle.loo.est) * 100

par(mfrow = c(1,2))

plot(known, mle.loo.est, pch = 16)
abline(0, 1, col = "red")

plot(known, killworth.slope.scale, pch = 16)
abline(0, 1, col = "red")


## With names

plot(known, mle.loo.est, col = "white")
text(known, mle.loo.est, labels = names(known))
abline(0, 1, col = "red")

plot(known, killworth.slope.scale, col = "white")
text(known, killworth.slope.scale, labels = names(known))
abline(0, 1, col = "red")






# Now do scaling without names --------------------------------------------

non.names.ind = c(1:N.k)[-c(1:12)]
known.non.names = known[non.names.ind]
n.non.names = length(non.names.ind)
slope.non.names.df = subset(slope.df, name %in% names(known)[non.names.ind])
mle.loo.non.names.est = mle.loo.est[non.names.ind]


## Using slope.df
killworth.slope.scale = mle.err.pred = rep(NA, n.non.names)
for(k in 1:n.non.names){
  slope.sub = slope.non.names.df[-k,]
  
  # Raw error ---------------------------------------------------------------
  
  ## Fit linear regression
  lm.fit.sub = lm(mle.err ~ slope, data = slope.sub)
  
  ## Predict log relative error
  logmle.err.pred = predict(lm.fit.sub, newdata = data.frame(slope = slope.df$slope[k]))
  
  ## Predict relative error
  mle.err.pred[k] = logmle.err.pred
  
  ## Correct for killworth estimate
  killworth.slope.scale[k] = mle.loo.non.names.est[k] * mle.err.pred[k] ## Multiply because inverse error
}


# Numerical Results -------------------------------------------------------


## Slope.scale
(
  Metrics::mape(known.non.names, mle.loo.non.names.est) -
    Metrics::mape(known.non.names, killworth.slope.scale)
) /
  Metrics::mape(known.non.names, mle.loo.non.names.est) * 100


par(mfrow = c(1,2))
plot(known.non.names, mle.loo.non.names.est, col = "white")
text(known.non.names, mle.loo.non.names.est, labels = names(known.non.names))
abline(0, 1, col = "red")

plot(known.non.names, killworth.slope.scale, col = "white")
text(known.non.names, killworth.slope.scale, labels = names(known.non.names))
abline(0, 1, col = "red")






rel.error.killworth = (known.non.names - mle.loo.non.names.est) / known.non.names
rel.error.scale = (known.non.names - killworth.slope.scale) / known.non.names

df.rel.error = data.frame(
  name = names(known.non.names),
  true = known.non.names,
  killworth = rel.error.killworth,
  scale = rel.error.scale
)
df.rel.error$name = toupper(df.rel.error$name)
df.rel.error$better = NA
df.rel.error$better[which(abs(df.rel.error$scale) <
                            abs(df.rel.error$killworth))] = "Yes"

df.melt.rel.error = reshape2::melt(df.rel.error,
                                   id.vars = c("name", "true", "better"))
df.melt.rel.error = subset(df.melt.rel.error,!is.na(true))
df.melt.rel.error$name = factor(df.melt.rel.error$name,
                                # levels = names(known))
                                levels = df.rel.error$name[order(known)])


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
    breaks = c("killworth", "scale"),
    labels = c("\nBasic Estimator \n", "\nAdjusted Estimator\n")
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
  "./Figures/McCarty_non_names_results.pdf",
  gg.res,
  width = 20,
  height = 10,
  units = "in"
)














# Now just names ----------------------------------------------------------

names.ind = c(1:12)
known.names = known[names.ind]
n.names = length(names.ind)
slope.names.df = subset(slope.df, name %in% names(known)[names.ind])
mle.loo.names.est = mle.loo.est[names.ind]


## Using slope.df
killworth.slope.scale = mle.err.pred = rep(NA, n.names)
for(k in 1:n.names){
  slope.sub = slope.names.df[-k,]
  
  # Raw error ---------------------------------------------------------------
  
  ## Fit linear regression
  lm.fit.sub = lm(mle.err ~ slope, data = slope.sub)
  
  ## Predict log relative error
  logmle.err.pred = predict(lm.fit.sub, newdata = data.frame(slope = slope.df$slope[k]))
  
  ## Predict relative error
  mle.err.pred[k] = logmle.err.pred
  
  ## Correct for killworth estimate
  killworth.slope.scale[k] = mle.loo.names.est[k] * mle.err.pred[k] ## Multiply because inverse error
}


# Numerical Results -------------------------------------------------------


## Slope.scale
(
  Metrics::mape(known.names, mle.loo.names.est) -
    Metrics::mape(known.names, killworth.slope.scale)
) /
  Metrics::mape(known.names, mle.loo.names.est) * 100


par(mfrow = c(1,2))
plot(known.names, mle.loo.names.est, col = "white")
text(known.names, mle.loo.names.est, labels = names(known.names))
abline(0, 1, col = "red")

plot(known.names, killworth.slope.scale, col = "white")
text(known.names, killworth.slope.scale, labels = names(known.names))
abline(0, 1, col = "red")








rel.error.killworth = (known.names - mle.loo.names.est) / known.names
rel.error.scale = (known.names - killworth.slope.scale) / known.names

df.rel.error = data.frame(
  name = names(known.names),
  true = known.names,
  killworth = rel.error.killworth,
  scale = rel.error.scale
)
df.rel.error$name = toupper(df.rel.error$name)
df.rel.error$better = NA
df.rel.error$better[which(abs(df.rel.error$scale) <
                            abs(df.rel.error$killworth))] = "Yes"

df.melt.rel.error = reshape2::melt(df.rel.error,
                                   id.vars = c("name", "true", "better"))
df.melt.rel.error = subset(df.melt.rel.error,!is.na(true))
df.melt.rel.error$name = factor(df.melt.rel.error$name,
                                # levels = names(known))
                                levels = df.rel.error$name[order(known)])


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
    breaks = c("killworth", "scale"),
    labels = c("\nBasic Estimator \n", "\nAdjusted Estimator\n")
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
  "./Figures/McCarty_names_results.pdf",
  gg.res,
  width = 20,
  height = 10,
  units = "in"
)
















# Now no twin or diabetes -------------------------------------------------


keep.ind = (1:N.k)[-c(21, 23)]
known.keep = known[keep.ind]
n.keep = length(keep.ind)
slope.keep.df = subset(slope.df, name %in% names(known)[keep.ind])
mle.loo.keep.est = mle.loo.est[keep.ind]


## Using slope.df
killworth.slope.scale = mle.err.pred = rep(NA, n.keep)
for(k in 1:n.keep){
  slope.sub = slope.keep.df[-k,]
  
  # Raw error ---------------------------------------------------------------
  
  ## Fit linear regression
  lm.fit.sub = lm(mle.err ~ slope, data = slope.sub)
  
  ## Predict log relative error
  logmle.err.pred = predict(lm.fit.sub, newdata = data.frame(slope = slope.df$slope[k]))
  
  ## Predict relative error
  mle.err.pred[k] = logmle.err.pred
  
  ## Correct for killworth estimate
  killworth.slope.scale[k] = mle.loo.keep.est[k] * mle.err.pred[k] ## Multiply because inverse error
}


# Numerical Results -------------------------------------------------------


## Slope.scale
(
  Metrics::mape(known.keep, mle.loo.keep.est) -
    Metrics::mape(known.keep, killworth.slope.scale)
) /
  Metrics::mape(known.keep, mle.loo.keep.est) * 100


par(mfrow = c(1,2))
plot(known.keep, mle.loo.keep.est, col = "white")
text(known.keep, mle.loo.keep.est, labels = names(known.keep))
abline(0, 1, col = "red")

plot(known.keep, killworth.slope.scale, col = "white")
text(known.keep, killworth.slope.scale, labels = names(known.keep))
abline(0, 1, col = "red")








rel.error.killworth = (known.keep - mle.loo.keep.est) / known.keep
rel.error.scale = (known.keep - killworth.slope.scale) / known.keep

df.rel.error = data.frame(
  name = names(known.keep),
  true = known.keep,
  killworth = rel.error.killworth,
  scale = rel.error.scale
)
df.rel.error$name = toupper(df.rel.error$name)
df.rel.error$better = NA
df.rel.error$better[which(abs(df.rel.error$scale) <
                            abs(df.rel.error$killworth))] = "Yes"

df.melt.rel.error = reshape2::melt(df.rel.error,
                                   id.vars = c("name", "true", "better"))
df.melt.rel.error = subset(df.melt.rel.error,!is.na(true))
df.melt.rel.error$name = factor(df.melt.rel.error$name,
                                # levels = names(known))
                                levels = df.rel.error$name[order(known)])


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
    breaks = c("killworth", "scale"),
    labels = c("\nBasic Estimator \n", "\nAdjusted Estimator\n")
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
  "./Figures/McCarty_twin_diab_results.pdf",
  gg.res,
  width = 20,
  height = 10,
  units = "in"
)




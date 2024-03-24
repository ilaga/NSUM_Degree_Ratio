library(ggplot2)
library(reshape2)
library(ggpmisc)

library(dplyr)





## Read in data
knownpop = read.csv("../Quality_vs_quantity_code/data/rwanda-varnames.csv")
nsum = readRDS("./Data/qvsq_meal.RDS") ## Meal

subpop.names = knownpop$varname
subpop.names = toupper(subpop.names)
subpop.names = gsub(".", " ", subpop.names, fixed = TRUE)


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
names(known) = subpop.names
n.known = ncol(y)
known.ind = 1:n.known
known.sizes = known

degrees.known = N * rowSums(y) / sum(known)


y = as.data.frame(y)



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
  xlim(0.004, 0.022)

gg.slope


ggsave(
  "./Figures/Rwanda_meal_slopes.pdf",
  gg.slope,
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







rel.error.killworth = (known - mle.loo.est) / known
rel.error.scale = (known - killworth.slope.scale) / known

df.rel.error = data.frame(
  name = names(known),
  true = known,
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
  "./Figures/Rwanda_all_results.pdf",
  gg.res,
  width = 20,
  height = 10,
  units = "in"
)

















# Now do scaling without names --------------------------------------------

non.names.ind = c(1:10)
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
















# Now also without priest -------------------------------------------------

non.names.ind = c(1:5, 7:10)
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
                                levels = df.rel.error$name[order(known.non.names)])


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
  "./Figures/Rwanda_meal_non_names_results.pdf",
  gg.res,
  width = 20,
  height = 10,
  units = "in"
)





















# Now just names ----------------------------------------------------------

names.ind = (1:N.k)[-c(1:10)]
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
















# Now just without priest -------------------------------------------------
non.priest.ind = c(1:N.k)[-6]
known.non.priest = known[non.priest.ind]
n.non.priest = length(non.priest.ind)
slope.non.priest.df = subset(slope.df, name %in% names(known)[non.priest.ind])
mle.loo.non.priest.est = mle.loo.est[non.priest.ind]


## Using slope.df
killworth.slope.scale = mle.err.pred = rep(NA, n.non.priest)
for(k in 1:n.non.priest){
  slope.sub = slope.non.priest.df[-k,]
  
  # Raw error ---------------------------------------------------------------
  
  ## Fit linear regression
  lm.fit.sub = lm(mle.err ~ slope, data = slope.sub)
  
  ## Predict log relative error
  logmle.err.pred = predict(lm.fit.sub, newdata = data.frame(slope = slope.df$slope[k]))
  
  ## Predict relative error
  mle.err.pred[k] = logmle.err.pred
  
  ## Correct for killworth estimate
  killworth.slope.scale[k] = mle.loo.non.priest.est[k] * mle.err.pred[k] ## Multiply because inverse error
}


# Numerical Results -------------------------------------------------------


## Slope.scale
(
  Metrics::mape(known.non.priest, mle.loo.non.priest.est) -
    Metrics::mape(known.non.priest, killworth.slope.scale)
) /
  Metrics::mape(known.non.priest, mle.loo.non.priest.est) * 100


par(mfrow = c(1,2))
plot(known.non.priest, mle.loo.non.priest.est, col = "white")
text(known.non.priest, mle.loo.non.priest.est, labels = names(known.non.priest))
abline(0, 1, col = "red")

plot(known.non.priest, killworth.slope.scale, col = "white")
text(known.non.priest, killworth.slope.scale, labels = names(known.non.priest))
abline(0, 1, col = "red")








rel.error.killworth = (known.non.priest - mle.loo.non.priest.est) / known.non.priest
rel.error.scale = (known.non.priest - killworth.slope.scale) / known.non.priest

df.rel.error = data.frame(
  name = names(known.non.priest),
  true = known.non.priest,
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
                                levels = df.rel.error$name[order(known.non.priest)])


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
  "./Figures/Rwanda_non_priest_results.pdf",
  gg.res,
  width = 20,
  height = 10,
  units = "in"
)


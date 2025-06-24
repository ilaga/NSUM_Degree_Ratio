library(ggplot2)
library(reshape2)
library(ggpmisc)
library(dplyr)
library(igraph)
library(grafify)

set.seed(52)

## Initialize Parameters
N.i <- 20000
N.k <- 20
means <- diag(seq(0.25, 0.5, length.out = N.k))
means <- means + 0.05

## Simulate SBM data
sbm.sim <- sample_sbm(N.i, means, rep(N.i / N.k, N.k), directed = F, loops = F)
node.groups <- rep(1:N.k, each = N.i / N.k)
node.groups <- model.matrix(~ factor(node.groups) - 1)
colnames(node.groups) <- 1:N.k
adj.mat <- as.matrix(get.adjacency(sbm.sim))

## Calculate subpopulation sizes
known <- known.sizes <- colSums(node.groups)
known / N.i

which.group <- apply(node.groups, 1, function(x) {
  which(x == 1)
})

known.degree <- rowSums(adj.mat)

## Find true degree ratio
degree.ratio <- rep(NA, N.k)
for (k in 1:N.k) {
  degree.ratio[k] <- mean(known.degree[which(which.group == k)]) / mean(known.degree)
}
degree.ratio


## Get ARD
ard <- adj.mat %*% node.groups

N.k <- ncol(ard)
N.i <- nrow(ard)
N <- nrow(ard)

n.known <- N.k
names(ard) <- names(known) <- 1:N.k
y <- y.known <- ard


## Define weights
degrees.loo.est <- matrix(NA, nrow = N.i, ncol = N.k)
ratio.mat <- matrix(NA, nrow = N.i, ncol = N.k)
for (k in 1:N.k) {
  degrees.loo.est[, k] <- N * rowSums(y[, -k]) / sum(known[-k])
  ratio.mat[, k] <- y[, k]
}


## Scale
ratio.mat <- scale(ratio.mat, center = FALSE, scale = apply(ratio.mat, 2, mean, na.rm = TRUE))

## Plot ratio against estimated loo degree
df.deg <- data.frame(degrees.loo.est)
names(df.deg) <- names(known)
df.ratio <- data.frame(ratio.mat)
names(df.ratio) <- names(known)

degrees.melt <- melt(df.deg)
logratio.melt <- melt(df.ratio)
logratio.melt$degrees <- degrees.melt$value
names(logratio.melt)[2] <- "LogRatio"


## Check slopes
slopes <- rep(NA, N.k)
for (k in 1:N.k) {
  lm.fit.2 <- lm(df.ratio[, k] ~ -1 + df.deg[, k])
  slopes[k] <- coef(lm.fit.2)[1]
}
names(slopes) <- names(known)

# PIMLE/MLE LOO Estimates -------------------------------------------------


mle.loo.est <- rep(NA, N.k)
for (k in 1:N.k) {
  known.ind.sub <- (1:N.k)[-k]
  ## Step 1) Estimate the unknown subpopulation and degrees
  degrees.est <- N * rowSums(y[, -k]) / sum(known[-k])

  ## Step 2) Get LOO estimates
  mle.loo.est[k] <- N * sum(unlist(y[, k])) / sum(degrees.est)
}





# Calculate relative error ------------------------------------------------
mle.rel.err <- known / mle.loo.est


# Plot against intercepts and slopes --------------------------------------

slope.df <- data.frame(
  slope = slopes, name = names(known),
  mle.err = mle.rel.err, logmle.err = log(mle.rel.err)
)

slope.melt <- melt(slope.df, measure.vars = c("slope"))

gg.slope <- ggplot(slope.melt, aes(x = value, y = mle.err, label = name)) +
  geom_point(size = 3) +
  ylab(expression(N[k] / hat(N)[k])) +
  xlab(expression(hat(beta)[1])) +
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
  geom_smooth(method = "lm")

gg.slope


ggsave(
  "./Figures/SBM_slopes.pdf",
  gg.slope,
  width = 20,
  height = 10,
  units = "in"
)



# Perform scaling using slope ---------------------------------------------

## Using slope.df

killworth.slope.scale <- mle.err.pred <- rep(NA, n.known)
for (k in 1:n.known) {
  slope.sub <- slope.df[-k, ]
  # Raw error ---------------------------------------------------------------

  ## Fit linear regression
  lm.fit.sub <- lm(mle.err ~ slope, data = slope.sub)

  ## Predict log relative error
  logmle.err.pred <- predict(lm.fit.sub, newdata = data.frame(slope = slope.df$slope[k]))

  ## Predict relative error
  mle.err.pred[k] <- logmle.err.pred

  ## Correct for killworth estimate
  killworth.slope.scale[k] <- mle.loo.est[k] * mle.err.pred[k]
}


# Numerical Results -------------------------------------------------------


known.known <- known
## Slope.scale
(
  Metrics::mape(known.known, mle.loo.est) -
    Metrics::mape(known.known, killworth.slope.scale)
) /
  Metrics::mape(known.known, mle.loo.est) * 100



par(mfrow = c(1, 2))

plot(known.known ~ mle.loo.est, pch = 16)
abline(0, 1, col = "red")

plot(known.known ~ killworth.slope.scale, pch = 16)
abline(0, 1, col = "red")









# Results plot ------------------------------------------------------------



rel.error.killworth <- (known.known - mle.loo.est) / known.known
rel.error.scale <- (known.known - killworth.slope.scale) / known.known

df.rel.error <- data.frame(
  name = names(known),
  true = known.known,
  killworth = rel.error.killworth,
  scale = rel.error.scale
)
df.rel.error$name <- toupper(df.rel.error$name)
df.rel.error$better <- NA
df.rel.error$better[which(abs(df.rel.error$scale) <
  abs(df.rel.error$killworth))] <- "Yes"

df.melt.rel.error <- reshape2::melt(df.rel.error,
  id.vars = c("name", "true", "better")
)
df.melt.rel.error <- subset(df.melt.rel.error, !is.na(true))
df.melt.rel.error$name <- factor(df.melt.rel.error$name,
  levels = names(known.known)[order(known.known)]
)


gg.res <- ggplot(
  df.melt.rel.error,
  aes(
    x = name,
    y = value,
    col = variable,
    group = variable,
    label = name
  )
) +
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
  scale_colour_grafify(
    na.value = "transparent",
    breaks = c("killworth", "scale"),
    labels = c("\nBasic Estimator \n", "\nAdjusted Estimator\n")
  ) +
  guides(color = guide_legend(override.aes = list(shape = 19))) +
  geom_abline(
    intercept = 0,
    slope = 0,
    linewidth = 1
  ) +
  xlab("Subpopulation") +
  ylab("Relative Error") +
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
  "./Figures/SBM_results.pdf",
  gg.res,
  width = 20,
  height = 10,
  units = "in"
)

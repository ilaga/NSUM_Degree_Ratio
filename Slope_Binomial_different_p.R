library(ggplot2)
library(reshape2)
library(ggpmisc)
library(grafify)

set.seed(27)
N <- 10^7
N.i <- 10000
N.k <- 20

known <- round(runif(N.k, 10^3, 10^6))
known <- known[sample(1:N.k)] ## Reorder subpopulation size
d <- round(runif(N.i, 10, 1000))
names(known) <- 1:N.k

c.k.1 <- seq(-1 / max(abs(d - mean(d))), 1 / max(abs(d - mean(d))), length.out = N.k)
c.k.2 <- seq(-1 / max(abs(d^(2) - mean(d^(2)))), 1 / max(abs(d^(2) - mean(d^(2)))), length.out = N.k)
c.k.3 <- seq(-1 / max(abs(d^(1 / 2) - mean(d^(1 / 2)))), 1 / max(abs(d^(1 / 2) - mean(d^(1 / 2)))), length.out = N.k)
c.k.4 <- seq(-1 / max(abs(d^(1 / 3) - mean(d^(1 / 3)))), 1 / max(abs(d^(1 / 3) - mean(d^(1 / 3)))), length.out = N.k)
c.k.5 <- seq(-1 / max(abs(d^(1 / 4) - mean(d^(1 / 4)))), 1 / max(abs(d^(1 / 4) - mean(d^(1 / 4)))), length.out = N.k)
c.k.6 <- seq(-1 / max(abs(d^(-1) - mean(d^(-1)))), 1 / max(abs(d^(-1) - mean(d^(-1)))), length.out = N.k)


## Simulate y
y <- f.mat <- matrix(NA, nrow = N.i, ncol = N.k)
for (k in 1:N.k) {
  p.sample <- sample(c(-2, -1, 1, 2), 1)
  if (p.sample == -2) {
    f.mat[, k] <- (1 + (d^(2) - mean(d^(2))) * c.k.2[k])
  } else if (p.sample == -1) {
    f.mat[, k] <- (1 + (d^(-1) - mean(d^(-1))) * c.k.6[k])
  } else if (p.sample == 1) {
    f.mat[, k] <- (1 + (d - mean(d)) * c.k.1[k])
  } else if (p.sample == 2) {
    f.mat[, k] <- (1 + (d^(1 / 2) - mean(d^(1 / 2))) * c.k.3[k])
  }

  y[, k] <- rbinom(N.i, d, known[k] / N * f.mat[, k])
}


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
  "./Figures/Killworth_dif_p_slopes.pdf",
  gg.slope,
  width = 20,
  height = 10,
  units = "in"
)



# Perform scaling using slope ---------------------------------------------

## Using slope.df
killworth.slope.scale <- mle.err.pred <- rep(NA, N.k)
for (k in 1:N.k) {
  slope.sub <- slope.df[-k, ]

  ## Fit linear regression
  lm.fit.sub <- lm(mle.err ~ slope, data = slope.sub)

  ## Predict log relative error
  logmle.err.pred <- predict(lm.fit.sub, newdata = data.frame(slope = slope.df$slope[k]))

  ## Predict relative error
  mle.err.pred[k] <- logmle.err.pred

  ## Correct for killworth estimate
  killworth.slope.scale[k] <- mle.loo.est[k] * mle.err.pred[k] ## Multiply because inverse error
}


# Numerical Results -------------------------------------------------------


## Slope.scale
(
  Metrics::mape(known, mle.loo.est) -
    Metrics::mape(known, killworth.slope.scale)
) /
  Metrics::mape(known, mle.loo.est) * 100

par(mfrow = c(1, 2))

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








# Results plot ------------------------------------------------------------



rel.error.killworth <- (known - mle.loo.est) / known
rel.error.scale <- (known - killworth.slope.scale) / known

df.rel.error <- data.frame(
  name = names(known),
  true = known,
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
  levels = names(known)
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
  "./Figures/Binomial_dif_p_results.pdf",
  gg.res,
  width = 20,
  height = 10,
  units = "in"
)

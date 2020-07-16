rm(list = ls())

library(modeest)  # Under trial
library("ggplot2")
library("plotly")

# Get gbdcd regression function
source("scripts/gbdcd.R")

# Load test dataset here
load(file = sprintf("data/four_144/c_four_144_50.RData"))

X <- cbind(1, dt[, 2])
y <- dt[, 1]

out <- gbdcd(
  y = y,
  X = X,
  viz = neighbors,
  n_iterations = 100000,
  burn_in = 50000,
  c = 0.001,
  coeffs_mu_prior = rep(0, 2),
  sigma_prior = 0.5,
  lambda = 1e-6,
  plot = T
)

k_hat_vector <- out$k.MCMC
k_hat_mode <- mfv(k_hat_vector)[1]
k_hat_summary <- summary(k_hat_vector)
hpd <- quantile(k_hat_vector, probs = c(0.025, 0.975))
hpd_size <- hpd[2] - hpd[1]

unname(c(k_hat_mode, k_hat_summary, hpd, hpd_size))

# Assess results ----------------------------------------------------------

# Proximity dendogram
n_cluster <- mfv(out$k.MCMC)[1]

plot(out$cluster.info,
  ylab = "d(i~j)", main = "Proximity Dendrogram", xlab = "",
  labels = grid$idx, cex = 0.6
)
groups <- cutree(out$cluster.info, k = n_cluster)
rect.hclust(out$cluster.info, k = n_cluster, border = "red")

# Grid plot
grid$result <- groups
plot(y ~ x, data = grid, pch = 15, col = result, cex = 2)
text(x = grid$x, y = grid$y, labels = grid$idx, pos = 1)


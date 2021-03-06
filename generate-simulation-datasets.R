rm(list=ls())

# import cluster structure
grid_size <- 64
load("data/cluster_64.RData")

# Dataset parameters
n_simulations <- 100
sigma <- 0.5
Z <- abs(qnorm(0.95))
x_cluster_a <- seq(0, 1, length = grid_size/2)
x_cluster_b <- seq(-1, 0, length = grid_size/2)
k <- sqrt(1 / sum(x_cluster_a^2))

# Coefficients for the linear regression equations
beta_cluster_a <- +Z * k * sigma * sqrt(2)
beta_cluster_b <- -Z * k * sigma * sqrt(2)

y_cluster_a <- 0 + beta_cluster_a * x_cluster_a + rnorm(length(x_cluster_a), sd = sigma)
y_cluster_b <- 0 + beta_cluster_b * x_cluster_b + rnorm(length(x_cluster_b), sd = sigma)

# Create datasets with one one cluster ------------------------------------
for (s_index in 1:n_simulations) {
  x <- seq(-1, 1, length = 64)
  y <- rnorm(64, mean = 0, sd = sigma)
  dt <- cbind(y, x)

  save(centers, neighbors, grid, dt, file = sprintf("data/one_%d/c_one_%d.RData", grid_size, s_index))
}

# Create datasets with two diagonal clusters ------------------------------
for (s_index in 1:n_simulations) {
  x_mat <- matrix(data = 0, 8, 8)
  y_mat <- matrix(data = 0, 8, 8)
  y_cluster_a <- 0 + beta_cluster_a * x_cluster_a + rnorm(length(x_cluster_a), sd = sigma)
  y_cluster_b <- 0 + beta_cluster_b * x_cluster_b + rnorm(length(x_cluster_b), sd = sigma)

  # Distribute cluster positions
  for (i in 1:8)
  {
    for (j in 1:8)
    {
      if (i > j) {
        x_mat[i, j] <- 1
        y_mat[i, j] <- 1
      } else {
        x_mat[i, j] <- 2
        y_mat[i, j] <- 2
      }
    }
  }

  # Correct some values on the diagonal
  x_mat[1, 1] <- 1
  y_mat[1, 1] <- 1
  x_mat[3, 3] <- 1
  y_mat[3, 3] <- 1
  x_mat[5, 5] <- 1
  y_mat[5, 5] <- 1
  x_mat[7, 7] <- 1
  y_mat[7, 7] <- 1

  # Insert variable values
  x_mat[which(x_mat == 1)] <- x_cluster_a
  x_mat[which(x_mat == 2)] <- x_cluster_b
  y_mat[which(y_mat == 1)] <- y_cluster_a
  y_mat[which(y_mat == 2)] <- y_cluster_b

  x <- as.vector(t(x_mat))
  y <- as.vector(t(y_mat))
  dt <- cbind(y, x)

  save(centers, neighbors, grid, dt, file = sprintf("data/diagonal_%d/c_diagonal_%d.RData", grid_size, s_index))
}

# Create datasets with four square clusters -------------------------------
for (s_index in 1:n_simulations) {
  x_mat <- matrix(data = 0, 8, 8)
  y_mat <- matrix(data = 0, 8, 8)
  y_cluster_a <- 0 + beta_cluster_a * x_cluster_a + rnorm(length(x_cluster_a), sd = sigma)
  y_cluster_b <- 0 + beta_cluster_b * x_cluster_b + rnorm(length(x_cluster_b), sd = sigma)

  # Distribute cluster positions
  for (i in 1:8)
  {
    for (j in 1:8)
    {
      if (i < 5 && j < 5) {
        x_mat[i, j] <- 1
        y_mat[i, j] <- 1
      } else if (i > 4 && j < 5) {
        x_mat[i, j] <- 2
        y_mat[i, j] <- 2
      } else if (i < 5 && j > 4) {
        x_mat[i, j] <- 2
        y_mat[i, j] <- 2
      } else if (i > 4 && j > 4) {
        x_mat[i, j] <- 1
        y_mat[i, j] <- 1
      }
    }
  }

  # Insert variable values
  x_mat[which(x_mat == 1)] <- x_cluster_a
  x_mat[which(x_mat == 2)] <- x_cluster_b
  y_mat[which(y_mat == 1)] <- y_cluster_a
  y_mat[which(y_mat == 2)] <- y_cluster_b

  x <- as.vector(t(x_mat))
  y <- as.vector(t(y_mat))
  dt <- cbind(y, x)

  save(centers, neighbors, grid, dt, file = sprintf("data/four_%d/c_four_%d.RData", grid_size, s_index))
}

# Create datasets with four square clusters - 144 regions -----------------------
grid_size <- 144
load("data/cluster_144.RData")

# Dataset parameters
n_simulations <- 100
sigma <- 0.5
Z <- abs(qnorm(0.95))
x_cluster_a <- seq(0, 1, length = grid_size/2)
x_cluster_b <- seq(-1, 0, length = grid_size/2)
k <- sqrt(1 / sum(x_cluster_a^2))

# Coefficients for the linear regression equations
beta_cluster_a <- +Z * k * sigma * sqrt(2)
beta_cluster_b <- -Z * k * sigma * sqrt(2)

y_cluster_a <- 0 + beta_cluster_a * x_cluster_a + rnorm(length(x_cluster_a), sd = sigma)
y_cluster_b <- 0 + beta_cluster_b * x_cluster_b + rnorm(length(x_cluster_b), sd = sigma)

for (s_index in 1:n_simulations) {
  x_mat <- matrix(data = 0, 12, 12)
  y_mat <- matrix(data = 0, 12, 12)
  y_cluster_a <- 0 + beta_cluster_a * x_cluster_a + rnorm(length(x_cluster_a), sd = sigma)
  y_cluster_b <- 0 + beta_cluster_b * x_cluster_b + rnorm(length(x_cluster_b), sd = sigma)
  
  # Distribute cluster positions
  for (i in 1:12)
  {
    for (j in 1:12)
    {
      if (i < 7 && j < 7) {
        x_mat[i, j] <- 1
        y_mat[i, j] <- 1
      } else if (i > 6 && j < 7) {
        x_mat[i, j] <- 2
        y_mat[i, j] <- 2
      } else if (i < 7 && j > 6) {
        x_mat[i, j] <- 2
        y_mat[i, j] <- 2
      } else if (i > 6 && j > 6) {
        x_mat[i, j] <- 1
        y_mat[i, j] <- 1
      }
    }
  }
  
  # Insert variable values
  x_mat[which(x_mat == 1)] <- x_cluster_a
  x_mat[which(x_mat == 2)] <- x_cluster_b
  y_mat[which(y_mat == 1)] <- y_cluster_a
  y_mat[which(y_mat == 2)] <- y_cluster_b
  
  x <- as.vector(t(x_mat))
  y <- as.vector(t(y_mat))
  dt <- cbind(y, x)
  
  save(centers, neighbors, grid, dt, file = sprintf("data/four_144/c_four_144_%d.RData", s_index))
}


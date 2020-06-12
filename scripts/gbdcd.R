# Define gbdcd regression function
gbdcd <- function(y,
                  X,
                  viz,
                  n_iterations = 1000000,
                  burn_in = 500000,
                  c = 0.35,
                  coeffs_mu_prior = rep(0, 2),
                  sigma_prior = sqrt(2),
                  lambda = 0.001, 
                  plot = TRUE) {

  # Test included because the partition function compilation takes time
  if (!exists("rcpp_partition", mode = "function")) source("scripts/partition.R")

  # Define important simulation parameters
  n_regions <- length(y)
  k_prob_prior <- ((1 - c)^(1:n_regions)) / sum((1 - c)^(1:n_regions))
  k_mean_prior <- round(sum((1:n_regions) * k_prob_prior))
  sigma2_prior <- sigma_prior^2
  sigma2 <- sigma2_prior
  a_0 <- 2.1 # Inverse Gamma shape parameter
  b_0 <- 1.1 # Inverse Gamma scale parameter

  n_coeffs <- dim(X)[2] # Number of estimated coefficients

  # Teste para a variância a priori dos betas
  sigma2_coeffs_prior <- sigma2_prior * diag(2)

  # Bookkeeping variables
  vec_k <- rep(NA, n_iterations)
  vec_sigma2 <- rep(NA, n_iterations)
  vec_steps <- rep(NA, n_iterations)
  vec_centers <- rep(NA, n_iterations)
  vec_accept <- rep(0, n_iterations)
  mat_freq <- matrix(0, n_regions, n_regions)
  mat_coeffs <- matrix(NA, n_regions, n_coeffs)
  mat_coeffs_hat <- array(NA, dim = c(n_regions, n_coeffs, n_iterations))
  cluster_centers <- sample.int(n_regions, size = 2, replace = FALSE)
  cluster_partitions <- rcpp_partition(viz, cluster_centers)

  # Initialize all cluster with coefficients equal to 0
  mat_coeffs[cluster_centers, ] <- 0

  # Start execution bar
  pb <- txtProgressBar(min = 0, max = n_iterations, char = "=", title = "progress bar")

  # Start the Markov Chain Monte Carlo simulation
  for (i in 1:n_iterations)
  {
    k <- length(cluster_centers) # Number of cluster centers on the current state

    # Select step type
    if (k == 1) {
      step_choice <- sample(c("birth", "death", "update"), size = 1, prob = c(0.8, 0, 0.2))
    } else if (k == n_regions) {
      step_choice <- sample(c("birth", "death", "update"), size = 1, prob = c(0, 0.8, 0.2))
    } else {
      step_choice <- sample(c("birth", "death", "update"), size = 1, prob = c(0.4, 0.4, 0.2))
    }

    # Define birth step
    if (step_choice == "birth") {

      # Select potential new cluster center and update cluster configuration
      new_center <- sample((1:n_regions)[!((1:n_regions) %in% cluster_centers)], size = 1)
      new_cluster_centers <- append(cluster_centers, new_center, after = sample(0:length(cluster_centers), size = 1))
      new_cluster_partitions <- rcpp_partition(viz, new_cluster_centers)
      y_cluster <- subset(y, new_cluster_partitions == new_center)
      X_cluster <- subset(X, new_cluster_partitions == new_center)

      # Propose regression coefficients for the new cluster
      auxiliar_proposed <- solve(crossprod(X_cluster, X_cluster) + lambda * diag(2))
      coeffs_mean_proposed <- auxiliar_proposed %*% crossprod(X_cluster, y_cluster)
      coeffs_cov_proposed <- sigma2 * auxiliar_proposed

      mat_coeffs[new_center, ] <- mvtnorm::rmvnorm(1, mean = coeffs_mean_proposed, sigma = coeffs_cov_proposed)

      phi <- mvtnorm::dmvnorm(mat_coeffs[new_center, ], mean = coeffs_mean_proposed, sigma = coeffs_cov_proposed)
      phi_k_plus_1 <- mvtnorm::dmvnorm(mat_coeffs[new_center, ], mean = coeffs_mu_prior, sigma = sigma2_coeffs_prior)

      # Define likelihood ratio and step acceptance probability
      llk <- sum(dnorm(y, mean = rowSums(mat_coeffs[cluster_partitions, ] * X), sd = sqrt(sigma2), log = TRUE))
      llk_plus_1 <- sum(dnorm(y, mean = rowSums(mat_coeffs[new_cluster_partitions, ] * X), sd = sqrt(sigma2), log = TRUE))
      llk_ratio <- exp(llk_plus_1 - llk)
      accept_prob <- llk_ratio * (1 - c) * 1 * (phi_k_plus_1 / phi)

      if (is.na(accept_prob)) {
        accept_prob <- 0
      } # Handle errors on small probabilities

      alpha <- min(1, accept_prob)

      if (runif(1) < alpha) {
        cluster_centers <- new_cluster_centers
        cluster_partitions <- new_cluster_partitions
        vec_accept[i] <- 1
      }
    }

    # Define death step
    if (step_choice == "death") {

      # Select potential cluster center to remove and update cluster configuration
      center <- sample(cluster_centers, 1)
      new_cluster_centers <- cluster_centers[-which(cluster_centers == center)]
      new_cluster_partitions <- rcpp_partition(viz, new_cluster_centers)
      y_cluster <- subset(y, cluster_partitions == center)
      X_cluster <- subset(X, cluster_partitions == center)

      # Calculate coefficient distribution parameters for the proposed configuration
      auxiliar_proposed <- solve(crossprod(X_cluster, X_cluster) + lambda * diag(2))
      coeffs_mean_proposed <- auxiliar_proposed %*% crossprod(X_cluster, y_cluster)
      coeffs_cov_proposed <- sigma2 * auxiliar_proposed

      phi <- mvtnorm::dmvnorm(mat_coeffs[center, ], mean = coeffs_mean_proposed, sigma = coeffs_cov_proposed)
      phi_k_minus_1 <- mvtnorm::dmvnorm(mat_coeffs[center, ], mean = coeffs_mu_prior, sigma = sigma2_coeffs_prior)

      # Define likelihood ratio and step acceptance probability
      llk <- sum(dnorm(y, mean = rowSums(mat_coeffs[cluster_partitions, ] * X), sd = sqrt(sigma2), log = TRUE))
      llk_minus_1 <- sum(dnorm(y, mean = rowSums(mat_coeffs[new_cluster_partitions, ] * X), sd = sqrt(sigma2), log = TRUE))
      llk_ratio <- exp(llk_minus_1 - llk)
      accept_prob <- llk_ratio * (1 / (1 - c)) * 1 * (phi / phi_k_minus_1)

      if (is.na(accept_prob)) {
        accept_prob <- 0
      } # Handle errors on small probabilities

      alpha <- min(1, accept_prob)

      if (runif(1) < alpha) {
        cluster_centers <- new_cluster_centers
        cluster_partitions <- new_cluster_partitions
        vec_accept[i] <- 1
      }
    }

    # Define update step
    if (step_choice == "update") {

      # Update regression coefficients
      for (center_update in cluster_centers) {

        # Select the data points of an specific cluster
        y_cluster <- subset(y, cluster_partitions == center_update)
        X_cluster <- subset(X, cluster_partitions == center_update)

        # Update regression coefficients on an specific cluster
        auxiliar_proposed <- solve(crossprod(X_cluster, X_cluster) + lambda * diag(2))
        coeffs_mean_proposed <- auxiliar_proposed %*% crossprod(X_cluster, y_cluster)
        coeffs_cov_proposed <- sigma2 * auxiliar_proposed

        mat_coeffs[center_update, ] <- mvtnorm::rmvnorm(1, mean = coeffs_mean_proposed, sigma = coeffs_cov_proposed)
      }

      # Update variance parameter
      a_n <- a_0 + n_regions / 2
      b_n <- b_0 + sum((y - rowSums(mat_coeffs[cluster_partitions, ] * X))^2) / 2
      sigma2 <- 1 / rgamma(1, shape = a_n, rate = b_n)

      vec_accept[i] <- 1
    }

    # Record chain state
    vec_k[i] <- length(cluster_centers)
    vec_sigma2[i] <- sigma2
    vec_steps[i] <- step_choice
    mat_coeffs_hat[, , i] <- mat_coeffs[cluster_partitions, ]
    vec_centers[i] <- paste(cluster_centers, collapse = ";")

    # Set progress bar
    setTxtProgressBar(pb, i)

    # Fill 'connection' frequency matrix
    if (i > burn_in) mat_freq <- mat_freq + rcpp_freqmatrix(cluster_partitions)
  }

  # End of MCMC iterations
  Sys.sleep(1)
  close(pb)

  # Chain results processing - remove burn-in
  seq.burn <- -(1:burn_in)
  vec_k <- vec_k[seq.burn]
  vec_sigma2 <- vec_sigma2[seq.burn]
  vec_steps <- vec_steps[seq.burn]
  mat_coeffs_hat <- mat_coeffs_hat[, , seq.burn]
  vec_accept <- vec_accept[seq.burn]
  vec_centers <- vec_centers[seq.burn]

  # Posteriors for the coefficients and credible intervals
  coeff_hat <- apply(mat_coeffs_hat, MARGIN = c(1, 2), FUN = median)
  coeff_lwr <- apply(mat_coeffs_hat, MARGIN = c(1, 2), FUN = function(x) quantile(x, probs = 0.05))
  coeff_upr <- apply(mat_coeffs_hat, MARGIN = c(1, 2), FUN = function(x) quantile(x, probs = 0.95))

  # Processing 'connection' frequency matrix
  mat_connections <- mat_freq
  mat_connections[lower.tri(mat_connections)] <- t(mat_connections)[lower.tri(mat_connections)]
  maximum <- max(as.vector(mat_connections), na.rm = TRUE)
  mat_connections <- maximum - mat_connections
  diag(mat_connections) <- maximum + 1

  # Hierarchical clustering on 'connection' frequency data
  clusters <- hclust(as.dist(mat_connections), c("single", "complete")[1])
  
  # Plot results
  if (plot == TRUE) {
    
  plot(1:length(vec_k), vec_k, type = "l", xlab = "step", main = "k-MCMC")
  barplot(table(vec_k) / sum(table(vec_k)), col = "light blue", main = "k-Posteriori")
  
  cat("\n Estimates for k: \n")
  print(summary(vec_k))
  print(quantile(vec_k, probs = c(0.05, 0.95)))
  
  cat("\n Estimates for sigma2: \n")
  print(summary(vec_sigma2))
  print(quantile(vec_sigma2, probs = c(0.05, 0.95)))
  
  cat("Step acceptance frequency: \n \n")
  print(aggregate(vec_accept ~ vec_steps, FUN = mean))
  
  plot(clusters, ylab = "proximity", main = "Proximity Dendogram", xlab = "", cex = 0.7)
  }

  output <- list(
    coeff.info = list(coeff_lwr, coeff_hat, coeff_upr),
    cluster.info = clusters,
    matConexoes = mat_connections,
    k.MCMC = vec_k,
    vec.centros = vec_centers
  )

  return(output)
}

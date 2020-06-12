rm(list = ls())

library(doParallel)

# Register cluster for parallel execution
cl <- makeCluster(3)  # Number of nodes
registerDoParallel(cl)

# Read gbdcd regression function
source("scripts/gbdcd.R")

# Simulation parameters
number_of_tests <- 10
cluster_types <- c("one", "diagonal", "four")
c_prioris <- c(0.001, 0.333)
n_iterations <- 100000
burn_in <- 50000
coeffs_prior <- rep(0, 2)
lambda <- 1e-6

start_time <- Sys.time()
for (cluster_index in 1:length(cluster_types)) {
  for (priori_index in 1:length(c_prioris)) {
    
    result <- foreach(i = 1:number_of_tests, .combine = rbind) %dopar% {

      ## Load data
      load(file = sprintf(
        "data/%s/c_%s_%d.RData",
        cluster_types[cluster_index],
        cluster_types[cluster_index],
        i
      ))

      # Prepare input data
      X <- cbind(1, dt[, 2])
      y <- dt[, 1]

      # Run gbdcd
      out <- gbdcd(
        y = y,
        X = X,
        viz = neighbors,
        n_iterations = n_iterations,
        burn_in = burn_in,
        c = c_prioris[priori_index],
        coeffs_mu_prior = rep(0, 2),
        sigma_prior = 0.5,
        lambda = lambda,
        plot = F
      )

      k_hat_vector <- out$k.MCMC
      k_hat_mode <- modeest::mfv(k_hat_vector)[1]
      k_hat_summary <- summary(k_hat_vector)
      hpd <- quantile(k_hat_vector, probs = c(0.025, 0.975))
      hpd_size <- hpd[2] - hpd[1]

      unname(c(k_hat_mode, k_hat_summary, hpd, hpd_size)) # Result
    }
    
    # Set column names
    colnames(result) <- c(
      "mode_k",
      "min_k",
      "qt_1",
      "median_k",
      "mean_k",
      "qt_3",
      "max_k",
      "hpd_1",
      "hpd_2",
      "hpd_size"
    )

    save(result, file = sprintf(
      "results/result_cluster_%s_t_%d_c_%d_l_%s.RData",
      cluster_types[cluster_index],
      number_of_tests,
      priori_index,
      as.character(lambda)
    ))
  }
}
end_time <- Sys.time()
end_time - start_time

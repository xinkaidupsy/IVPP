#' Simulate data for a (multi-group) panelGVAR model
#'
#' This function generates a (multi-group) panel GVAR dataset
#'
#' @param n_node number of nodes
#' @param n_person an integer denoting the sample size of each group, default to 500
#' @param n_time number of waves
#' @param n_group number of groups
#'
#' @return A list of temporal and contemporaneous networks
#' %>%
#' @importFrom dplyr %>%

# To-dos
# add PDC list
sim_panelGVAR <- function(temp_ls,
                          beta_ls,
                          omega_ls,
                          n_node,
                          n_person = 500,
                          n_time,
                          n_group){



  # setup default -----------------------------------------------------------

  # * denotes Hadamard product
  beta = sqrt(t(PDC)^2 * (diag(sigma) %o% diag(kappa)) / (1-t(PDC)^2))

  # Generate the means (same for all variables in a wave):
  means <- seq_len(n_time) * mean_trend

  # the identity matrix for the networks
  I <- diag(n_node)

  # Variable names:
  varnames <- paste0("V",seq_len(n_node))

  # Generate a random between person correlation matrix for each group:
  between <- cov2cor(clusterGeneration::genPositiveDefMat(n_node)$Sigma)

  # generate data -----------------------------------------------------------


  for (g in 1:n_group) {

    # ----- sigma_zeta and beta -----

    # sigma_zeta is the sigma (cov matrix) of innovations zeta (contemporaneous)
    # P4 of https://www.tandfonline.com/doi/full/10.1080/00273171.2018.1454823
    sigma_zeta_per_wave[[paste0("g", g)]] <-
      lapply(seq_len(n_time),function(t){
        cov2cor(solve(I-cont_ls[[g]][[t]]))
      }) %>% setNames(paste0("t",seq_len(n_time)))

    # beta is the transpose of temporaneous
    beta_per_wave[[paste0("g", g)]] <-
      lapply(temp_ls[[g]], t) %>%
      setNames(paste0("t",seq_len(n_time-1)))

    # stationary distribution for t1
    stationary_wave_1[[paste0("g", g)]] <- matrix(solve(kronecker(I, I) - kronecker(beta_per_wave[[g]][[1]], beta_per_wave[[g]][[1]])) %*% c(sigma_zeta_per_wave[[g]][[1]]), n_node, n_node)


    # ----- data generation -----

    # we obtain this by calculating sigma_eta (latent variance) from the sigma_zeta (innovation variance) and beta matrix (cross-wave cov matrix)
    # Page 211 of https://link.springer.com/article/10.1007/s11336-020-09697-3

    stationary_wave_1[[paste0("g", g)]] <- matrix(solve(kronecker(I, I) - kronecker(beta_per_wave[[g]][[1]], beta_per_wave[[g]][[1]])) %*% c(sigma_zeta_per_wave[[g]][[1]]), n_node, n_node)

    # Simulate within-person data per wave. First wave 1:
    within_person_data_per_wave[[paste0("g", g)]] <- list(
      rmvnorm(n_person, sigma = stationary_wave_1[[g]])
    )

    # Simulate the other waves:
    for (t in 2:n_time){
      within_person_data_per_wave[[g]][[t]] <- t(beta_per_wave[[g]][[t-1]] %*% t(within_person_data_per_wave[[g]][[t-1]])) + rmvnorm(n_person, sigma = sigma_zeta_per_wave[[g]][[t]])
    }

    # Generate the means per person:
    person_means_base[[paste0("g", g)]] <- rmvnorm(n_person,sigma = between)
    person_means_per_wave[[paste0("g", g)]] <- lapply(seq_len(n_time),function(t)person_means_base[[g]] + means[t])

    # Add the data together:
    data[[paste0("g", g)]] <- lapply(seq_len(n_time), function(t) {
      data <- data.table(person_means_per_wave[[g]][[t]] + within_person_data_per_wave[[g]][[t]])
      setnames(data, varnames)
      data[, `:=` (subject = (g-1)*(n_person) + seq_len(.N),
                   time = t,
                   group = g)]
    })

    # Merge data across waves (long format):
    data_long[[paste0("g", g)]] <- rbindlist(data[[g]])[order(subject, time)]



  } # end: for(g in 1:n_group)

  # convert list to data.table
  data_merged <- rbindlist(data_long)

  # return data
  return(data_merged)
}

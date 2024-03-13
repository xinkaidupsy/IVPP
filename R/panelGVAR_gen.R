#' This function generates a (multi-group) panel GVAR dataset
#'
#' @param n_node number of nodes
#' @param n_person sample size
#' @param n_time number of waves
#' @param n_group number of groups
#' @param p_rewire extent of group difference
#' @param mean_trend mean trends in data
#' @param p_rewire_temp probability rewiring temporal networks across groups
#' @param p_rewire_cont probability rewiring temporal networks across groups
#'
#' @return A dataframe / datatable that has n_node variables, with subject id, wave number, and group membership specified (if multiple groups)
#' @export
#'
#' @importFrom bootnet genGGM
#' @importFrom clusterGeneration genPositiveDefMat
#' @importFrom dplyr %>%

panelGVAR_gen <- function(n_node = 6,
                          n_person = 500,
                          n_time = 4,
                          p_rewire = 0.5,
                          n_group = 1,
                          mean_trend = 0,
                          p_rewire_temp = 0,
                          p_rewire_cont = 0){

  # Initialization --------------------------------------------------------


  # Generate the means (same for all variables in a wave):
  means <- seq_len(n_time) * mean_trend

  # variances are all 1

  # Generate a random between person correlation matrix for each group:
  between <- cov2cor(clusterGeneration::genPositiveDefMat(n_node)$Sigma)

  # the identity matrix for the networks
  I <- diag(n_node)

  # Variable names:
  varnames <- paste0("V",seq_len(n_node))

  # base networks for each group
  cont_base_ls <-
    temp_base_ls <-
    # list for all networks
    cont_ls <-
    temp_ls <-

    # sigma_zeta, beta, and wave 1 distribution
    sigma_zeta_per_wave <-
    beta_per_wave <-
    stationary_wave_1 <-
    within_person_data_per_wave <-
    person_means_base <-
    person_means_per_wave <-
    data <-
    data_long <- list()


  # generate networks -------------------------------------------------------


  # the contemporaneous network for g1
  cont_base_ls[["g1"]] <- bootnet::genGGM(n_node, propPositive = 0.8)

  # the temporal network for g1
  temp_base_ls[["g1"]] <- diag(0.5, n_node)
  temp_base_ls[["g1"]][cbind(1:n_node, ((1:n_node) + 1) %% n_node + 1)] <-
    sample(c(0.25, -0.25), n_node, replace = TRUE, prob = c(0.8, 0.2))


  # If there are multiple groups, generate the rest
  if (n_group > 1) {

    # Use a previous network as the base to rewire instead of always using the first network
    # Start with second network
    for (g in 2:n_group) {
      # cont_base for g2 onwards
      cont_base_ls[[paste0("g",g)]] <- rewire(cont_base_ls[[g-1]], p = p_rewire, directed = FALSE)
      # temp_base for g2 onwards
      temp_base_ls[[paste0("g",g)]] <- rewire(temp_base_ls[[g-1]], p = p_rewire, directed = TRUE)
    }

  } # end: if(n_group > 1)


  # generate data -----------------------------------------------------------


  for (g in 1:n_group) {


    # ----- all true temp and cont networks -----


    # cont per wave for group g
    cont_ls[[paste0("g",g)]] <-
      # # rewire with probability p_rewire_cont if there is non-stationarity
      lapply(seq_len(n_time),function(x){
        rewire(cont_base_ls[[g]], p = p_rewire_cont, directed = FALSE)
      }) %>%
      # set dimension names
      setNames(paste0("t",seq_len(n_time)))

    # does not rewire first wave of each group
    cont_ls[[g]][[1]] <- cont_base_ls[[g]]

    # temp per wave for group g
    temp_ls[[paste0("g",g)]] <-
      # rewire with probability p_rewire_temp if there is non-stationarity
      lapply(seq_len(n_time-1),function(x){
        rewire(temp_base_ls[[g]], p = p_rewire_temp, directed = TRUE)
      }) %>%
      # set dimension names
      setNames(paste0("t",seq_len(n_time-1)))

    # does not rewire first wave of each group
    temp_ls[[g]][[1]] <- temp_base_ls[[g]]

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



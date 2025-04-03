#' Simulate data for a (multi-group) panelGVAR model
#'
#' This function generates data for the input (multi-group) panelGVAR model
#'
#' @param temp_base_ls a list of temporal networks of all groups
#' @param beta_base_ls a list of beta matrices of all groups
#' @param cont_base_ls a list of contemporaneous networks of all groups
#' @param n_node number of nodes
#' @param n_person an integer denoting the sample size of each group, default to 500
#' @param n_time number of waves, default to 3
#' @param n_group number of groups
#' @param mean_trend a numeric value indicating the extent of mean trends in data, default to 0
#' @param p_rewire_temp a numeric value between 0 and 1 indicating the extent of non-stationarity in temporal networks, default to 0
#' @param p_rewire_cont a numeric value between 0 and 1 indicating the extent of non-stationarity in contemporaneous networks, default to 0
#' @param save_nets a logical value indicating whether to save the data-generating networks, default to FALSE
#'
#' @author Xinkai Du
#' Maintainer: Xinkai Du <xinkai.du.xd@gmail.com>
#'
#' @returns A list of temporal and contemporaneous networks
#' @import mvtnorm dplyr
#' @importFrom stats cov2cor setNames
#' @export sim_panelGVAR
#' @examples
#' library(IVPP)
#' # Generate the network
#' net_ls <- gen_panelGVAR(n_node = 6,
#'                         p_rewire_temp = 0.5,
#'                         p_rewire_cont = 0,
#'                         n_group = 3)
#'
#' # Generate the data
#' data <- sim_panelGVAR(temp_base_ls = net_ls$temporal,
#'                       cont_base_ls = net_ls$omega_zeta_within,
#'                       n_person = 500,
#'                       n_time = 4,
#'                       n_group = 3,
#'                       n_node = 6)


sim_panelGVAR <- function(temp_base_ls,
                          beta_base_ls,
                          cont_base_ls,
                          # kappa_base_ls,
                          n_node,
                          n_person = 500,
                          n_time = 3,
                          n_group,
                          mean_trend = 0,
                          p_rewire_temp = 0,
                          p_rewire_cont = 0,
                          save_nets = FALSE){



  # set defaults ------------------------------------------------------------

  # # need data-generating temporal network to proceed
  # if(missing(temp_base_ls) &
  #    missing(beta_base_ls) &
  #    missing(PDC_base_ls)){
  #   stop("please input the list of temporal networks of all groups
  #        to one of temp_base_ls, beta_base_ls, and PDC_base_ls")
  # }

  # need data-generating contemporaneous network to proceed
  if(missing(cont_base_ls)){
    stop("please input the list of contemporaneous networks of all groups to cont_base_ls")
  }


  # obtain temp_base_ls if missing
  if(missing(temp_base_ls)){

    # can obtain directly from beta
    if(!missing(beta_base_ls)){

      temp_base_ls <- lapply(beta_base_ls, t)

    } else {
      stop("please input either temp_base_ls or beta_base_ls to proceed")
    }

  } # end: if(missing(temp_base_ls))


  # n_node
  if(missing(n_node)){
    n_node <- cont_base_ls[[1]] %>% nrow
  }

  n_n <- cont_base_ls[[1]] %>% nrow
  if(n_node != n_n){
    stop("please specify the correct number of nodes")
  }

  # n_groups
  if(missing(n_group)){
    n_group <- length(cont_base_ls)
  }

  n_g <- length(cont_base_ls)
  if(n_group != n_g){
    stop("please specify the correct number of groups")
  }




  # ----- Initialize means and tool vectors -----


  # Generate the means (same for all variables in a wave):
  means <- seq_len(n_time) * mean_trend

  # the identity matrix for the networks
  I <- diag(n_node)

  # Variable names:
  varnames <- paste0("V",seq_len(n_node))

  # Generate a random between person correlation matrix for each group:
  between <- cov2cor(clusterGeneration::genPositiveDefMat(n_node)$Sigma)


  # ----- network for each wave -----

  # initialize lists
  cont_ls <-
    beta_ls <-
    temp_ls <-
    # distributions
    # sigma_zeta, beta, and wave 1 distribution
    sigma_zeta_per_wave <-
    beta_per_wave <-
    stationary_wave_1 <-
    within_person_data_per_wave <-
    person_means_base <-
    person_means_per_wave <-
    data <-
    data_long <- list()

  # ----- all true temp, cont networks, and beta -----

  # cont list for all waves and all groups
  cont_ls <- lapply(seq_len(n_group), function(g){

    # create cont for each wave in group g
    cont_ls <- lapply(seq_len(n_time),function(x){
      # rewire across wave with probability p_rewire_cont
      rewire(cont_base_ls[[g]], p = p_rewire_cont, directed = FALSE)

    }) %>% setNames(paste0("t",seq_len(n_time))) # set dimension names

    # does not rewire first wave
    cont_ls[[1]] <- cont_base_ls[[g]]

    return(cont_ls)

  }) %>% setNames(paste0("g",seq_len(n_group)))

  # temp list of all waves for all groups
  temp_ls <- lapply(seq_len(n_group), function(g){
    temp_ls <- lapply(seq_len(n_time-1),function(t){
      # rewire across wave with probability p_rewire_temp
      rewire(temp_base_ls[[g]], p = p_rewire_temp, directed = TRUE)
    }) %>% setNames(paste0("t",seq_len(n_time-1))) # set dimension names

    # does not rewire first wave of each group
    temp_ls[[1]] <- temp_base_ls[[g]]

    return(temp_ls)

  }) %>% setNames(paste0("g",seq_len(n_group)))


  # beta is the transpose of temp
  beta_ls <- lapply(seq_len(n_group), function(g){
    beta_ls <- lapply(temp_ls[[g]], t) %>% setNames(paste0("t", seq_len(n_time-1)))
  }) %>% setNames(paste0("g",seq_len(n_group)))


  # ----- generate data -----


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
      mvtnorm::rmvnorm(n_person, sigma = stationary_wave_1[[g]])
    )

    # Simulate the other waves:
    for (t in 2:n_time){
      within_person_data_per_wave[[g]][[t]] <- t(beta_per_wave[[g]][[t-1]] %*% t(within_person_data_per_wave[[g]][[t-1]])) + mvtnorm::rmvnorm(n_person, sigma = sigma_zeta_per_wave[[g]][[t]])
    }

    # Generate the means per person:
    person_means_base[[paste0("g", g)]] <- mvtnorm::rmvnorm(n_person,sigma = between)
    person_means_per_wave[[paste0("g", g)]] <- lapply(seq_len(n_time),function(t)person_means_base[[g]] + means[t])

    # Add the data together:
    data[[paste0("g", g)]] <- lapply(seq_len(n_time),function(t){
      data <- as.data.frame(person_means_per_wave[[g]][[t]] + within_person_data_per_wave[[g]][[t]])
      names(data) <- varnames
      data$subject <- (g-1)*(n_person) + seq_len(nrow(data))
      data$time <- t
      data$group <- g
      data
    })

    # Merge data (long format):
    data_long[[paste0("g", g)]] <- bind_rows(data[[g]]) %>% arrange(.data$subject,.data$time,.data$group)

  } # end: for(g in 1:n_group)

  # convert list to data.table
  data_merged <- bind_rows(data_long)

  if (save_nets) {

    # compute kappa for all waves in each group
    kappa_ls <- lapply(cont_ls, function(g) {
      lapply(g, function(omega){

        # implicit delta, assume variance = 1 for all variables
        kappa <- solve(cov2cor(solve(diag(n_node) - omega)))
        # correct floating error
        kappa[abs(omega) < sqrt(.Machine$double.eps) & !diag(n_node)==1] <- 0

        return(kappa)

      })
    })

    # compute PDC
    PDC_ls <- lapply(seq_len(n_group), function(g){

      lapply(seq_len(n_time - 1), function(t){
        kappa <- kappa_ls[[g]][[t]] # here kappa of wave t is not used
        sigma <- solve(kappa)
        beta <- beta_ls[[g]][[t]]
        PDC <- t(beta / sqrt(diag(sigma) %o% diag(kappa) + beta^2))

        return(PDC)
      }) %>% setNames(paste0("t", seq_len(n_time-1)))

    }) %>% setNames(paste0("g", seq_len(n_group))) # set names g1-3

    nets <- list(
      beta = beta_ls,
      temporal = temp_ls,
      PDC = PDC_ls,
      kappa = kappa_ls,
      omega_zeta_within = cont_ls
    )
    nets_n_data <- list(
      nets = nets,
      data = data_merged
    )

    return(nets_n_data)

  } else {
    # return data
    return(data_merged)
  } # end: if(save_nets)

} # end: sim_panelGVAR

#' Simulate data for a (multi-group) N = 1 GVAR model
#'
#' This function generates data for the input (multi-group) N = 1 GVAR model
#'
#' @param temp_base_ls a list of temporal networks of all individuals
#' @param beta_base_ls a list of beta matrices of all individuals
#' @param cont_base_ls a list of contemporaneous networks of all individuals
#' @param kappa_base_ls a list of precision matricies of all individuals
#' @param n_node number of nodes
#' @param n_time number of measurements per person, default to 50
#' @param n_person number of individuals to generate data for
#' @param save_nets a logical value indicating whether to save the data-generating networks, default to FALSE
#'
#' @author Xinkai Du
#' Maintainer: Xinkai Du <xinkai.du.xd@gmail.com>
#'
#' @return A list of temporal and contemporaneous networks
#' @import mvtnorm dplyr
#' @importFrom graphicalVAR graphicalVARsim
#' @export sim_tsGVAR
#' @examples
#' library(IVPP)
#' # Generate the network
#' net_ls <- gen_tsGVAR(n_node = 6,
#'                      p_rewire_temp = 0.5,
#'                      p_rewire_cont = 0,
#'                      n_persons = 3)
#'
#' # Generate the data
#' data <- sim_tsGVAR(beta_base_ls = net_ls$beta,
#'                    kappa_base_ls = net_ls$kappa,
#'                    # n_person = 3,
#'                    n_time = 50)

sim_tsGVAR <- function(temp_base_ls,
                       beta_base_ls,
                       cont_base_ls,
                       kappa_base_ls,
                       n_node,
                       n_time = 50,
                       n_person,
                       save_nets = FALSE){

  # ----- set defaults -----

  # obtain kappa
  if(missing(kappa_base_ls)){

    # if contemporaneous networks are available, obtain directly
    if (!missing(cont_base_ls)) {

      # n_node
      if(missing(n_node)){
        n_node <- cont_base_ls[[1]] %>% nrow
      }

      n_n <- cont_base_ls[[1]] %>% nrow
      if(n_node != n_n){
        stop("please specify the correct number of nodes")
      }

      # obtain kappa from omega
      kappa_base_ls <- lapply(cont_base_ls, function(omega) {

        # implicit delta, assume variance = 1 for all variables
        kappa <- solve(cov2cor(solve(diag(n_node) - omega)))
        # correct floating error
        kappa[abs(omega) < sqrt(.Machine$double.eps) & !diag(n_node)==1] <- 0

        return(kappa)

      })

      # end: if(!missing(cont_base_ls))
    } else {
      stop("please input either kappa_base_ls or cont_base_ls to proceed")
    }

  } # end: if(missing(kappa_base_ls))

  # obtain beta
  if(missing(beta_base_ls)){

    # if temporal networks are available, obtain directly
    if (!missing(temp_base_ls)) {

      beta_base_ls <- lapply(temp_base_ls, t)

    } else {
      stop("please input either temp_base_ls, or beta_base_ls to proceed")
    }

  } # end: if(missing(beta_base_ls))

  # n_node
  if(missing(n_node)){
    n_node <- kappa_base_ls[[1]] %>% nrow
  }

  n_n <- kappa_base_ls[[1]] %>% nrow
  if(n_node != n_n){
    stop("please specify the correct number of nodes")
  }

  # n_person
  if(missing(n_person)){
    n_person <- length(kappa_base_ls)
  }

  n_p <- length(kappa_base_ls)
  if(n_person != n_p){
    stop("please specify the correct number of individuals")
  }

  # ----- generate data-----

  # initialize lists
  data <- list()

  data <- lapply(seq_len(n_person), function(p){

    # set.seed(1)

    # Simulate the data
    data <- graphicalVAR::graphicalVARsim(
      nTime = n_time,
      beta = beta_base_ls[[p]],
      kappa = kappa_base_ls[[p]]
      ) %>% as.data.frame

    data$id <- p

    return(data)

  })

  data_merged <- bind_rows(data)


    # return data
    return(data_merged)

} # end: sim_panelGVAR


#' Generate a (multi-group) panelGVAR model
#'
#' This function generates a (multi-group) panel GVAR model. Currently it's focused on generating
#'
#' @param n_node an integer denoting the number of nodes
#' @param n_time an integer denoting the number of waves
#' @param n_group an integer denoting the number of groups
#' @param p_rewire a numeric value between 0-1 denoting the extent of group difference
#' @param mean_trend a numeric value denoting the extent of mean trends in data
#' @param p_rewire_temp a numeric value between 0-1 denoting the probability rewiring temporal networks across groups
#' @param p_rewire_cont a numeric value between 0-1 denoting the probability rewiring temporal networks across groups
#'
#' @return A list of beta, PDC, kappa and contemporaneous networks
#'
#' @details
#' \code{beta} can be transposed to obtain the temporal network;
#' \code{PDC} is the partial directed correlation matrix, which is a standardized version of temporal network;
#' \code{kappa} is the precision matrix denoting conditional (in)dependence,
#' which is a inverse of covariance matrix denoting the (dependence) among variables;
#' kappa can be further standardized to the contemporaneous networks (\code{omega_zeta_within})
#'
#'
#' @importFrom bootnet genGGM
#' @importFrom clusterGeneration genPositiveDefMat
#' @importFrom dplyr %>%


gen_panelGVAR <- function(n_node = 6,
                          n_time = 4,
                          p_rewire = 0.5,
                          n_group = 1){


  # check inputs ------------------------------------------------------------

  if(p_rewire < 0 | p_rewire > 1){
    stop("probability of rewiring should be between 0 and 1")
  }

  # Initialization --------------------------------------------------------

  # base networks for each group
  cont_base_ls <- temp_base_ls <- list()

  # generate networks -------------------------------------------------------

  # the temporal network for g1
  temp_base_ls[["g1"]] <- diag(0.5, n_node)
  temp_base_ls[["g1"]][cbind(1:n_node, ((1:n_node) + 1) %% n_node + 1)] <-
    sample(c(0.25, -0.25), n_node, replace = TRUE, prob = c(0.8, 0.2))

  # the contemporaneous network for g1
  cont_base_ls[["g1"]] <- bootnet::genGGM(n_node, propPositive = 0.8)



  # If there are multiple groups, generate the rest
  if (n_group > 1) {

    # Use a previous network as the base to rewire instead of always using the first network
    # Start with second network
    for (g in 2:n_group) {
      # temp_base for g2 onwards
      temp_base_ls[[paste0("g",g)]] <- rewire(temp_base_ls[[g-1]], p = p_rewire, directed = TRUE)
      # cont_base for g2 onwards
      cont_base_ls[[paste0("g",g)]] <- rewire(cont_base_ls[[g-1]], p = p_rewire, directed = FALSE)
    }

  } # end: if(n_group > 1)

  # ----- all true temp and cont networks -----

  # precision matrix
  kappa_ls <- lapply(cont_base_ls, function(omega) {

      # implicit delta, assume variance = 1 for all variables
      kappa <- solve(cov2cor(solve(diag(n_node) - omega)))
      # correct floating error
      kappa[abs(omega) < sqrt(.Machine$double.eps) & !diag(n_node)==1] <- 0

      return(kappa)

  })

  # beta is the transpose of temp
  beta_base_ls <- lapply(temp_base_ls, t)

  # compute PDC
  PDC_ls <- lapply(seq_len(n_group), function(g){

      kappa <- kappa_ls[[g]]
      sigma <- solve(kappa)
      beta <- beta_base_ls[[g]]
      PDC <- t(beta / sqrt(diag(sigma) %o% diag(kappa) + beta^2))

      return(PDC)

  }) %>% setNames(paste0("g", seq_len(n_group))) # set names g1-3

  nets <- list(
    beta = beta_base_ls,
    temporal = temp_base_ls,
    PDC = PDC_ls,
    kappa = kappa_ls,
    omega_zeta_within = cont_base_ls
  )

  return(nets)

}

#' Rewire Networks
#'
#' This function takes in a network and rewire it with probability p
#'
#' @param net the network matrix that you want to rewire
#' @param directed logicals. Write TRUE if the input net is a directed network and FALSE otherwise
#' @param p probability of rewiring
#'
#' @return A network rewired from the input
#'

rewire <- function(net, p = 0, directed){
  net <- as.matrix(net)
  rownames(net) <- colnames(net) <- NULL
  if (missing(directed)){
    directed <- !all(net == t(net))
  }

  if (directed){
    edges <- c(net)
  } else {
    edges <- net[lower.tri(net,diag=FALSE)]
  }

  n_edge <- length(edges)

  for (i in which(edges!=0 & runif(n_edge) < p)){
    rewire_to <- sample(seq_len(n_edge)[-i],1)
    edge1 <- edges[i]
    edge2 <- edges[rewire_to]
    edges[i] <- edge2
    edges[rewire_to] <- edge1
  }

  if (directed){
    net[] <- edges
  } else {
    net[lower.tri(net,diag=FALSE)] <- edges
    net[upper.tri(net,diag=FALSE)] <-  t(net)[upper.tri(net,diag=FALSE)]
  }

  return(net)
}

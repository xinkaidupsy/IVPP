#' Generate a (multi-group) panelGVAR model
#'
#' This function generates a (multi-group) panel GVAR model. Currently generating temporal and contemporaneous networks
#'
#' @param n_node an integer denoting the number of nodes
#' @param n_group an integer denoting the number of groups
#' @param p_rewire_temp a numeric value between 0-1 denoting the extent of group difference in the temporal network
#' @param p_rewire_cont a numeric value between 0-1 denoting the extent of group difference in the contemporaneous network
#'
#' @author Xinkai Du
#' Maintainer: Xinkai Du <xinkai.du.xd@gmail.com>
#'
#' @return A list of beta, PDC, kappa and contemporaneous networks
#'
#' @details
#' \code{beta} can be transposed to obtain the temporal network;
#' \code{PDC} is the partial directed correlation matrix, which is a standardized version of temporal network;
#' \code{kappa} is the precision matrix denoting conditional (in)dependence,
#' which is a inverse of covariance matrix denoting the (dependence) among variables;
#' kappa can be further standardized to the contemporaneous networks (\code{omega_zeta_within})
#' @import bootnet dplyr clusterGeneration
#' @importFrom stats setNames cov2cor
#' @export gen_panelGVAR
#' @examples
#' library(IVPP)
#' # Generate the network
#' net_ls <- gen_panelGVAR(n_node = 6,
#'                         p_rewire_temp = 0.5,
#'                         p_rewire_cont = 0,
#'                         n_group = 2)



gen_panelGVAR <- function(n_node = 6,
                          p_rewire_temp = 0.5,
                          p_rewire_cont = 0.5,
                          n_group = 1){


  # check inputs ------------------------------------------------------------

  if(p_rewire_temp < 0 | p_rewire_temp > 1 | p_rewire_cont < 0 | p_rewire_cont > 1){
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
      temp_base_ls[[paste0("g",g)]] <- rewire(temp_base_ls[[g-1]], p = p_rewire_temp, directed = TRUE)
      # cont_base for g2 onwards
      cont_base_ls[[paste0("g",g)]] <- rewire(cont_base_ls[[g-1]], p = p_rewire_cont, directed = FALSE)
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

#' Generate time-series GVAR model for multiple (heterogeneous) individuals
#'
#' This function generates time-series GVAR model for multiple individuals that demonstrates difference or simularity.
#' Currently generating temporal and contemporaneous networks
#'
#' @param n_node an integer denoting the number of nodes
#' @param n_persons an integer denoting the number of individuals to generate tsGVAR for
#' @param p_rewire_temp a numeric value between 0-1 denoting the extent of individual difference in the temporal network
#' @param p_rewire_cont a numeric value between 0-1 denoting the extent of individual difference in the contemporaneous network
#'
#' @author Xinkai Du
#' Maintainer: Xinkai Du <xinkai.du.xd@gmail.com>
#'
#' @return A list of beta, PDC, kappa and contemporaneous networks
#'
#' @details
#' \code{beta} can be transposed to obtain the temporal network;
#' \code{PDC} is the partial directed correlation matrix, which is a standardized version of temporal network;
#' \code{kappa} is the precision matrix denoting conditional (in)dependence,
#' which is a inverse of covariance matrix denoting the (dependence) among variables;
#' kappa can be further standardized to the contemporaneous networks (\code{omega_zeta_within})
#' @import bootnet dplyr clusterGeneration
#' @importFrom stats setNames cov2cor
#' @export gen_tsGVAR
#' @examples
#' library(IVPP)
#'
#' # Generate the network
#' net_ls <- gen_tsGVAR(n_node = 6,
#'                      p_rewire_temp = 0.5,
#'                      p_rewire_cont = 0,
#'                      n_persons = 2)


gen_tsGVAR <- function(n_node = 6,
                       p_rewire_temp = 0.5,
                       p_rewire_cont = 0.5,
                       n_persons = 1){


  # check inputs ------------------------------------------------------------

  if(p_rewire_temp < 0 | p_rewire_temp > 1 | p_rewire_cont < 0 | p_rewire_cont > 1){
    stop("probability of rewiring should be between 0 and 1")
  }

  # Initialization --------------------------------------------------------

  # base networks for each group
  cont_base_ls <- temp_base_ls <- list()

  # generate networks -------------------------------------------------------

  # the temporal network for g1
  temp_base_ls[["p1"]] <- diag(0.5, n_node)
  temp_base_ls[["p1"]][cbind(1:n_node, ((1:n_node) + 1) %% n_node + 1)] <-
    sample(c(0.25, -0.25), n_node, replace = TRUE, prob = c(0.8, 0.2))

  # the contemporaneous network for g1
  cont_base_ls[["p1"]] <- bootnet::genGGM(n_node, propPositive = 0.8)



  # If there are multiple groups, generate the rest
  if (n_persons > 1) {

    # Use a previous network as the base to rewire instead of always using the first network
    # Start with second network
    for (g in 2:n_persons) {
      # temp_base for g2 onwards
      temp_base_ls[[paste0("p",g)]] <- rewire(temp_base_ls[[g-1]], p = p_rewire_temp, directed = TRUE)
      # cont_base for g2 onwards
      cont_base_ls[[paste0("p",g)]] <- rewire(cont_base_ls[[g-1]], p = p_rewire_cont, directed = FALSE)
    }

  } # end: if(n_persons > 1)

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
  PDC_ls <- lapply(seq_len(n_persons), function(g){

    kappa <- kappa_ls[[g]]
    sigma <- solve(kappa)
    beta <- beta_base_ls[[g]]
    PDC <- t(beta / sqrt(diag(sigma) %o% diag(kappa) + beta^2))

    return(PDC)

  }) %>% setNames(paste0("p", seq_len(n_persons))) # set names g1-3

  nets <- list(
    beta = beta_base_ls,
    temporal = temp_base_ls,
    PDC = PDC_ls,
    kappa = kappa_ls,
    omega_zeta = cont_base_ls
  )

  return(nets)

}


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

  for (i in which(edges!=0 & stats::runif(n_edge) < p)){
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

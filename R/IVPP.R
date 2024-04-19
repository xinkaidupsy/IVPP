
#' The omnibus test of the invariance partial pruning (IVPP) algorithm for multi-group panel GVAR models
#'
#' This function implements the IVPP algorithm to compare networks in the multi-group panel GVAR models.
#' The IVPP algorithm is a two-step procedure that first conducts an omnibus test of network difference and then performs partial pruning for the specific edge-level differences.
#' This function returns the results of the omnibus test.
#' Currently only supports the comparison of temporal and contemporaneous networks.
#'
#' @param data A data frame containing the long-formated panel data
#' @param vars A character vector of variable names
#' @param idvar A character string specifying subject IDs
#' @param beepvar A character string specifying the name of wave (time) variable
#' @param groups A character string specifying the name of group variable
#' @param constrain A character vector specifying the network you want to constrain equal in the omnibus test.
#' Specify "both" if you want to constrain on both temporal or contemporaneous networks to be equal across groups.
#' Specify "temporal" if you want to constrain the temporal network to be equal.
#' Specify "contemporaneous" if you want to constrain the contemporaneous network to be equal.
#' See \link{Details} for more information.
#' @param net_type A character vector specifying the type of networks to be compared.
#' Specify "saturated" if you want to estimate and compare the saturated networks.
#' Specify "sparse" if you want to estimate and compare the pruned networks.
#' @param prune_alpha A numeric value specifying the alpha level for the pruning (if net_type is sparse).
#' @param estimator A character string specifying the estimator to be used. Must be "FIML"
#' @param standardize A character string specifying the type of standardization to be used.
#' "none" (default) for no standardization, "z" for z-scores,
#' and "quantile" for a non-parametric transformation to the quantiles of the marginal standard normal distribution.
#' @param ... Additional arguments to be passed to the \code{\link[psychonetrics]{dlvm1}} function
#'
#' @details
#' The comparison between the fully unconstrained (free) model and tempEq model is a test for group equality in temporal networks.
#' The comparison between fully constrained model (bothEq) and tempEq is a test for group equality in contemporaneous networks.
#' Similarly, the comparison between the free model and contEq model is a test for group equality in contemporaneous networks,
#' and the comparison between bothEq and contEq is a test for group equality in temporal networks.
#'
#' @return The results of comparison between models of different levels of constraints
#'
#' @import dplyr
#' @import tidyr
#' @import psychonetrics

# # Generate the network
# net_ls <- gen_panelGVAR(n_node = 6,
#                         n_time = 4,
#                         p_rewire = 0.5,
#                         n_group = 3)
#
# # Generate the data
# data_merged <- sim_panelGVAR(temp_base_ls = net_ls$temporal,
#                              cont_base_ls = net_ls$omega_zeta_within,
#                              n_person = 500,
#                              n_time = 4,
#                              n_group = 3,
#                              n_node = 6)
#
# vars <- paste0("V",1:6)
# idvar = "subject"
# beepvar = "time"
# groups = "group"
# standardize = "z"
# estimator = "FIML"
#
# panel_omni(data = data_merged,
#            vars = paste0("V",1:6),
#            idvar = "subject",
#            beepvar = "time",
#            groups = "group",
#            net = "temporal",
#            net_type = "saturated",
#            estimator = "FIML",
#            standardize = "z")

# logic:
# do the omnibus 1
# if not significant, stop and return no difference
# if signigicant, omibus 2
# Do partial pruning for the significant

panel_omni <- function(data,
                       vars,
                       idvar,
                       beepvar,
                       groups,
                       # test = c("omnibus", "partial_prune"),
                       constrain = c("both", "temporal", "contemporaneous"),
                       # vsModel = c("bothEq", "free"),
                       net_type = c("saturated", "sparse"),
                       prune_alpha = 0.01,
                       estimator = "FIML",
                       standardize = c("none", "z","quantile"),
                       ...){

  # ----- check missing arguments -----
  if(missing(data)){
    stop("data is missing")
  }

  # vars
  if(missing(vars)){
    stop("specify the variable names")
  }

  # id
  if(missing(idvar)){
    stop("specify the id variable")
  }

  # beep
  if(missing(beepvar)){
    stop("specify the beep (wave) variable")
  }

  # group
  if(missing(group)){
    stop("specify the group variable")
  }

  # net_type
  if(missing(net_type)){
    stop("specify the type of networks to be compared")
  }

  # estimator
  if (estimator != "FIML"){
    stop("Only 'FIML' supported currently.")
  }

  # ----- saturated -----

  # estimate the saturated free model
  mod_saturated <- ml_gvar(data_merged,
                           vars = vars,
                           idvar = idvar,
                           beepvar = beepvar,
                           groups = groups,
                           standardize = standardize,
                           estimator = estimator,
                           between = "chol",
                           ...) %>% runmodel %>% suppressWarnings



  if(net_type == "saturated"){

    # model comparisons
    if(constrain == "both"){

      # estimate the fully-constrained model
      mod_saturated_bothEq <- mod_saturated %>%
        groupequal(matrix = "beta") %>%
        groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

      # compare with the free model
      comp <- psychonetrics::compare(free = mod_saturated,
                                     bothEq = mod_saturated_bothEq)

      mod_ls <- list(free = mod_saturated,
                     bothEq = mod_saturated_bothEq)

      # end: if(net == "both")
    } else if (constrain == "temporal"){

      # estimate the model that constrains temporal networks to be equal
      mod_saturated_tempEq <- mod_saturated %>%
        groupequal(matrix = "beta") %>% runmodel %>% suppressWarnings

      # compare with the free model and fully-constrained model
      comp <- psychonetrics::compare(free = mod_saturated,
                                     tempEq = mod_saturated_tempEq,
                                     bothEq = mod_saturated_bothEq)

      mod_ls <- list(free = mod_saturated,
                     tempEq = mod_saturated_tempEq,
                     bothEq = mod_saturated_bothEq)

    # end: if(net == "temporal")
    } else if (constrain == "contemporaneous"){

      mod_saturated_contEq <- mod_saturated %>%
        groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

      comp <- psychonetrics::compare(free = mod_saturated,
                                     contEq = mod_saturated_contEq,
                                     bothEq = mod_saturated_bothEq)

      mod_ls <- list(free = mod_saturated,
                     contEq = mod_saturated_contEq,
                     bothEq = mod_saturated_bothEq)

    # end: if(net == "contemporaneous")
    } else {
      stop("constrain either 'both', 'temporal', or 'contemporaneous'")
    } # end: if(net == xxx)


  # end: if(net_type == "saturated")
  } else { # if (net_type == "pruned")

    # The free union model
    mod_union <- mod_saturated %>% prune(alpha = prune_alpha) %>% unionmodel %>% runmodel %>% suppressWarnings

    # the fully-constrained union model
    mod_union_bothEq <- mod_union %>%
      groupequal(matrix = "beta") %>%
      groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

    # model comparisons
    if(constrain == "both"){

      # compare with the free model
      comp <- psychonetrics::compare(free = mod_union,
                                     bothEq = mod_union_bothEq)

      mod_ls <- list(free = mod_union,
                     bothEq = mod_union_bothEq)

      # end: if(net == "both")
    } else if (constrain == "temporal"){

      # estimate the model that constrains temporal networks to be equal
      mod_union_tempEq <- mod_union %>%
        groupequal(matrix = "beta") %>% runmodel %>% suppressWarnings

      # compare with the free model and fully-constrained model
      comp <- psychonetrics::compare(free = mod_union,
                                     tempEq = mod_union_tempEq,
                                     bothEq = mod_union_bothEq)

      mod_ls <- list(free = mod_union,
                     tempEq = mod_union_tempEq,
                     bothEq = mod_union_bothEq)

      # end: if(net == "temporal")
    } else if (constrain == "contemporaneous"){

      # estimate the model that constrains contemporaneous networks to be equal
      mod_union_contEq <- mod_union %>%
        groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

      # compare with the free model and fully-constrained model
      comp <- psychonetrics::compare(free = mod_union,
                                     contEq = mod_union_contEq,
                                     bothEq = mod_union_bothEq)

      mod_ls <- list(free = mod_union,
                     contEq = mod_union_contEq,
                     bothEq = mod_union_bothEq)

    # end: if(net == "contemporaneous")
    } else {
      stop("constrain either 'both', 'temporal', or 'contemporaneous'")
    } # end: if(net == xxx)


  } # end: else (net_type == "pruned")

  # return the comparison results

  return(list(comp_res = comp, mod = mod_ls))

} # end: panel_omni


#' The partial_pruning test of the invariance partial pruning (IVPP) algorithm for multi-group panel GVAR models
#'
#' This function implements the IVPP algorithm to compare networks in the multi-group panel GVAR models.
#' The IVPP algorithm is a two-step procedure that first conducts an omnibus test of network difference and then performs partial pruning for the specific edge-level differences.
#' This function returns the results of the omnibus test.
#' Currently only supports the comparison of temporal and contemporaneous networks.
#'
#' @param data A data frame containing the long-formated panel data
#' @param vars A character vector of variable names
#' @param idvar A character string specifying subject IDs
#' @param beepvar A character string specifying the name of wave (time) variable
#' @param groups A character string specifying the name of group variable
#' @param prune_net A character vector specifying the network you want to partial prune on.
#' Specify "both" if you want to partial prune both temporal or contemporaneous networks.
#' Specify "temporal" if you want to partial prune the temporal network.
#' Specify "contemporaneous" if you want to partial prune the contemporaneous network.
#' @param estimator A character string specifying the estimator to be used. Must be "FIML"
#' @param standardize A character string specifying the type of standardization to be used.
#' "none" (default) for no standardization, "z" for z-scores,
#' and "quantile" for a non-parametric transformation to the quantiles of the marginal standard normal distribution.
#' @param ... Additional arguments to be passed to the \code{\link[psychonetrics]{dlvm1}} function
#'
#' @details
#' The comparison between the fully unconstrained (free) model and tempEq model is a test for group equality in temporal networks.
#' The comparison between fully constrained model (bothEq) and tempEq is a test for group equality in contemporaneous networks.
#' Similarly, the comparison between the free model and contEq model is a test for group equality in contemporaneous networks,
#' and the comparison between bothEq and contEq is a test for group equality in temporal networks.
#'
#' @return The results of comparison between models of different levels of constraints
#'
#' @import dplyr
#' @import tidyr
#' @import psychonetrics
#'
#'

panel_pp <- function(data,
                     vars,
                     idvar,
                     beepvar,
                     groups,
                     # test = c("omnibus", "partial_prune"),
                     prune_net = c("both", "temporal", "contemporaneous"),
                     # vsModel = c("bothEq", "free"),
                     prune_alpha = 0.01,
                     estimator = "FIML",
                     standardize = c("none", "z","quantile"),
                     ...){
  # ----- check missing arguments -----
  if(missing(data)){
    stop("data is missing")
  }

  # vars
  if(missing(vars)){
    stop("specify the variable names")
  }

  # id
  if(missing(idvar)){
    stop("specify the id variable")
  }

  # beep
  if(missing(beepvar)){
    stop("specify the beep (wave) variable")
  }

  # group
  if(missing(group)){
    stop("specify the group variable")
  }

  # net_type
  if(missing(net_type)){
    stop("specify the type of networks to be compared")
  }

  # estimator
  if (estimator != "FIML"){
    stop("Only 'FIML' supported currently.")
  }

  # ----- saturated -----

  # estimate the saturated free model
  mod_saturated <- ml_gvar(data_merged,
                           vars = vars,
                           idvar = idvar,
                           beepvar = beepvar,
                           groups = groups,
                           standardize = standardize,
                           estimator = estimator,
                           between = "chol",
                           ...) %>% runmodel %>% suppressWarnings



  if(net_type == "saturated"){

    # model comparisons
    if(constrain == "both"){

      # estimate the fully-constrained model
      mod_saturated_bothEq <- mod_saturated %>%
        groupequal(matrix = "beta") %>%
        groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

      # compare with the free model
      comp <- psychonetrics::compare(free = mod_saturated,
                                     bothEq = mod_saturated_bothEq)

      # end: if(net == "both")
    } else if (constrain == "temporal"){

      # estimate the model that constrains temporal networks to be equal
      mod_saturated_tempEq <- mod_saturated %>%
        groupequal(matrix = "beta") %>% runmodel %>% suppressWarnings

      # compare with the free model and fully-constrained model
      comp <- psychonetrics::compare(free = mod_saturated,
                                     tempEq = mod_saturated_tempEq,
                                     bothEq = mod_saturated_bothEq)

      # end: if(net == "temporal")
    } else if (constrain == "contemporaneous"){

      mod_saturated_contEq <- mod_saturated %>%
        groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

      comp <- psychonetrics::compare(free = mod_saturated,
                                     contEq = mod_saturated_contEq,
                                     bothEq = mod_saturated_bothEq)

      # end: if(net == "contemporaneous")
    } else {
      stop("constrain either 'both', 'temporal', or 'contemporaneous'")
    } # end: if(net == xxx)


    # end: if(net_type == "saturated")
  } else { # if (net_type == "pruned")

    # The free union model
    mod_union <- mod_saturated %>% prune(alpha = prune_alpha) %>% unionmodel %>% runmodel %>% suppressWarnings

    # the fully-constrained union model
    mod_union_bothEq <- mod_union %>%
      groupequal(matrix = "beta") %>%
      groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

    # model comparisons
    if(constrain == "both"){

      # compare with the free model
      comp <- psychonetrics::compare(free = mod_union,
                                     bothEq = mod_union_bothEq)

      # end: if(net == "both")
    } else if (constrain == "temporal"){

      # estimate the model that constrains temporal networks to be equal
      mod_union_tempEq <- mod_union %>%
        groupequal(matrix = "beta") %>% runmodel %>% suppressWarnings

      # compare with the free model and fully-constrained model
      comp <- psychonetrics::compare(free = mod_union,
                                     tempEq = mod_union_tempEq,
                                     bothEq = mod_union_bothEq)

      # end: if(net == "temporal")
    } else if (constrain == "contemporaneous"){

      # estimate the model that constrains contemporaneous networks to be equal
      mod_union_contEq <- mod_union %>%
        groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

      # compare with the free model and fully-constrained model
      comp <- psychonetrics::compare(free = mod_union,
                                     contEq = mod_union_contEq,
                                     bothEq = mod_union_bothEq)

      # end: if(net == "contemporaneous")
    } else {
      stop("constrain either 'both', 'temporal', or 'contemporaneous'")
    } # end: if(net == xxx)


  } # end: else (net_type == "pruned")

  # return the comparison results
  return(comp)

} # end: panel_pp

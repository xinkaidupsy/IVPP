
#' The omnibus test of the invariance partial pruning (IVPP) algorithm for multi-group panel GVAR models
#'
#' This function implements the IVPP algorithm to compare networks in the multi-group panel GVAR models.
#' The IVPP algorithm is a two-step procedure that first conducts an omnibus test of network difference and then performs partial pruning for the specific edge-level differences.
#' Currently supports the comparison of temporal and contemporaneous networks.
#'
#' @param data A data frame containing the long-formated panel data
#' @param vars A character vector of variable names
#' @param idvar A character string specifying subject IDs
#' @param beepvar A character string specifying the name of wave (time) variable
#' @param groups A character string specifying the name of group variable
#' @param test A character vector specifying the network you want to test group-equality on in the omnibus test.
#' Specify "both" if you want to test on both temporal or contemporaneous networks.
#' Specify "temporal" if you want to test only on the temporal network.
#' Specify "contemporaneous" if you want to test only on the contemporaneous network.
#' See \link{Details} for more information.
#' @param net_type A character vector specifying the type of networks to be compared.
#' Specify "saturated" if you want to estimate and compare the saturated networks.
#' Specify "sparse" if you want to estimate and compare the pruned networks.
#' @param prune_alpha A numeric value specifying the alpha level for the pruning (if net_type = "sparse").
#' @param partial_prune A logical value specifying whether to conduct partial pruning test or not.
#' @param p_prune_alpha A numeric value specifying the alpha level for the partial pruning (if partial_prune = TRUE).
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
#' @return A list containing the results of IVPP and networks of all groups.
#'
#' @import dplyr
#' @import tidyr
#' @import psychonetrics

IVPP_panel <- function(data,
                       vars,
                       idvar,
                       beepvar,
                       groups,
                       # test = c("omnibus", "partial_prune"),
                       test = c("both", "temporal", "contemporaneous"),
                       # vsModel = c("bothEq", "free"),
                       net_type = c("saturated", "sparse"),
                       partial_prune,
                       prune_alpha = 0.01,
                       p_prune_alpha = 0.01,
                       estimator = "FIML",
                       standardize = c("none", "z","quantile"),
                       ...){

  # ----- argument check -----

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
  if(missing(groups)){
    stop("specify the group variable")
  }

  # test
  if(!(test %in% c("both", "temporal", "contemporaneous"))){
    stop("test should be either on 'both', 'temporal', or 'contemporaneous'")
  }

  # net_type
  if(missing(net_type)){
    stop("specify the type of networks to be compared to net_type")
  }

  if(!(net_type %in% c("saturated", "sparse"))){
    stop("network_type is either 'saturated or 'sparse'")
  }

  # partial_prune
  if(missing(partial_prune)){
    stop("specify whether to conduct partial pruning or not")
  }

  if(!is.logical(partial_prune)){
    stop("partial_prune should be a logical value")
  }

  # estimator
  if (estimator != "FIML"){
    stop("Only 'FIML' supported currently.")
  }

  # ----- omnibus test -----

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


  # omnibus test for saturated & sparse networks
  if(net_type == "saturated"){

    # estimate the fully-constrained model
    mod_saturated_bothEq <- mod_saturated %>%
      groupequal(matrix = "beta") %>%
      groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

    # model comparisons
    if(test == "both"){

      # compare with the free model
      comp_vs_free <- comp_vs_bothEq <- psychonetrics::compare(free = mod_saturated,
                                                               bothEq = mod_saturated_bothEq)

      # generate the comparison table
      tab_vs_free <- data.frame(
        AIC = comp_vs_free$AIC,
        delta_AIC = c(NA, diff(comp_vs_free$AIC)),
        BIC = comp_vs_free$BIC,
        delta_BIC = c(NA, diff(comp_vs_free$BIC)),
        row.names = comp_vs_free$model
        # chi_sq = comp_vs_free$Chisq,
        # delta_chisq = comp_vs_free$Chisq_diff,
        # p_value = comp_vs_free$p_value
      )

      tab_vs_bothEq <- data.frame(
        AIC = comp_vs_bothEq$AIC,
        delta_AIC = c(NA, diff(comp_vs_bothEq$AIC)),
        BIC = comp_vs_bothEq$BIC,
        delta_BIC = c(NA, diff(comp_vs_bothEq$BIC)),
        row.names = comp_vs_bothEq$model
        # chi_sq = comp_vs_bothEq$Chisq,
        # delta_chisq = comp_vs_bothEq$Chisq_diff,
        # p_value = comp_vs_bothEq$p_value
      )



      # end: if(test == "both")
    } else if (test == "temporal"){

      # estimate the model that constrains temporal networks to be equal
      mod_saturated_tempEq <- mod_saturated %>%
        groupequal(matrix = "beta") %>% runmodel %>% suppressWarnings

      # estimate the model that constrains contemporaneous networks to be equal
      mod_saturated_contEq <- mod_saturated %>%
        groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

      # compare with the free model
      comp_vs_free <- psychonetrics::compare(free = mod_saturated,
                                     tempEq = mod_saturated_tempEq)

      tab_vs_free <- data.frame(
        AIC = comp_vs_free$AIC,
        delta_AIC = c(NA, diff(comp_vs_free$AIC)),
        BIC = comp_vs_free$BIC,
        delta_BIC = c(NA, diff(comp_vs_free$BIC)),
        row.names = comp_vs_free$model
        # chi_sq = comp_vs_free$Chisq,
        # delta_chisq = comp_vs_free$Chisq_diff,
        # p_value = comp_vs_free$p_value
      )

      # compare with the full-constrained model
      comp_vs_bothEq <- psychonetrics::compare(contEq = mod_saturated_contEq,
                                               bothEq = mod_saturated_bothEq)

      tab_vs_bothEq <- data.frame(
        AIC = comp_vs_bothEq$AIC,
        delta_AIC = c(NA, diff(comp_vs_bothEq$AIC)),
        BIC = comp_vs_bothEq$BIC,
        delta_BIC = c(NA, diff(comp_vs_bothEq$BIC)),
        row.names = comp_vs_bothEq$model
        # chi_sq = comp_vs_bothEq$Chisq,
        # delta_chisq = comp_vs_bothEq$Chisq_diff,
        # p_value = comp_vs_bothEq$p_value
      )



    # end: if(test == "temporal")
    } else if (test == "contemporaneous"){

      # estimate the model that constrains temporal networks to be equal
      mod_saturated_tempEq <- mod_saturated %>%
        groupequal(matrix = "beta") %>% runmodel %>% suppressWarnings

      # estimate the model that constrains contemporaneous networks to be equal
      mod_saturated_contEq <- mod_saturated %>%
        groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

      # compare with the free model
      comp_vs_free <- psychonetrics::compare(free = mod_saturated,
                                             contEq = mod_saturated_contEq)

      tab_vs_free <- data.frame(
        AIC = comp_vs_free$AIC,
        delta_AIC = c(NA, diff(comp_vs_free$AIC)),
        BIC = comp_vs_free$BIC,
        delta_BIC = c(NA, diff(comp_vs_free$BIC)),
        row.names = comp_vs_free$model
        # chi_sq = comp_vs_free$Chisq,
        # delta_chisq = comp_vs_free$Chisq_diff,
        # p_value = comp_vs_free$p_value
      )

      # compare with the full-constrained model
      comp_vs_bothEq <- psychonetrics::compare(tempEq = mod_saturated_tempEq,
                                               bothEq = mod_saturated_bothEq)

      tab_vs_bothEq <- data.frame(
        AIC = comp_vs_bothEq$AIC,
        delta_AIC = c(NA, diff(comp_vs_bothEq$AIC)),
        BIC = comp_vs_bothEq$BIC,
        delta_BIC = c(NA, diff(comp_vs_bothEq$BIC)),
        row.names = comp_vs_bothEq$model
        # chi_sq = comp_vs_bothEq$Chisq,
        # delta_chisq = comp_vs_bothEq$Chisq_diff,
        # p_value = comp_vs_bothEq$p_value
      )

    # end: if(test == "contemporaneous")
    }


  # end: if(net_type == "saturated")
  } else { # if (net_type == "pruned")

    # The free union model
    mod_union <- mod_saturated %>% prune(alpha = prune_alpha) %>% unionmodel %>% runmodel %>% suppressWarnings

    # the fully-constrained union model
    mod_union_bothEq <- mod_union %>%
      groupequal(matrix = "beta") %>%
      groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

    # model comparisons
    if(test == "both"){

      # compare with the free model
      comp_vs_free <- comp_vs_bothEq <- psychonetrics::compare(free = mod_union,
                                                               bothEq = mod_union_bothEq)

      # generate the comparison table
      tab_vs_free <- data.frame(
        AIC = comp_vs_free$AIC,
        delta_AIC = c(NA, diff(comp_vs_free$AIC)),
        BIC = comp_vs_free$BIC,
        delta_BIC = c(NA, diff(comp_vs_free$BIC)),
        row.names = comp_vs_free$model
        # chi_sq = comp_vs_free$Chisq,
        # delta_chisq = comp_vs_free$Chisq_diff,
        # p_value = comp_vs_free$p_value
      )

      tab_vs_bothEq <- data.frame(
        AIC = comp_vs_bothEq$AIC,
        delta_AIC = c(NA, diff(comp_vs_bothEq$AIC)),
        BIC = comp_vs_bothEq$BIC,
        delta_BIC = c(NA, diff(comp_vs_bothEq$BIC)),
        row.names = comp_vs_bothEq$model
        # chi_sq = comp_vs_bothEq$Chisq,
        # delta_chisq = comp_vs_bothEq$Chisq_diff,
        # p_value = comp_vs_bothEq$p_value
      )



      # end: if(test == "both")
    } else if (test == "temporal"){

      # estimate the model that constrains temporal networks to be equal
      mod_union_tempEq <- mod_union %>%
        groupequal(matrix = "beta") %>% runmodel %>% suppressWarnings

      # estimate the model that constrains contemporaneous networks to be equal
      mod_union_contEq <- mod_union %>%
        groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

      # compare with the free model
      comp_vs_free <- psychonetrics::compare(free = mod_union,
                                             tempEq = mod_union_tempEq)

      tab_vs_free <- data.frame(
        AIC = comp_vs_free$AIC,
        delta_AIC = c(NA, diff(comp_vs_free$AIC)),
        BIC = comp_vs_free$BIC,
        delta_BIC = c(NA, diff(comp_vs_free$BIC)),
        row.names = comp_vs_free$model
        # chi_sq = comp_vs_free$Chisq,
        # delta_chisq = comp_vs_free$Chisq_diff,
        # p_value = comp_vs_free$p_value
      )

      # compare with the full-constrained model
      comp_vs_bothEq <- psychonetrics::compare(contEq = mod_union_contEq,
                                               bothEq = mod_union_bothEq)

      tab_vs_bothEq <- data.frame(
        AIC = comp_vs_bothEq$AIC,
        delta_AIC = c(NA, diff(comp_vs_bothEq$AIC)),
        BIC = comp_vs_bothEq$BIC,
        delta_BIC = c(NA, diff(comp_vs_bothEq$BIC)),
        row.names = comp_vs_bothEq$model
        # chi_sq = comp_vs_bothEq$Chisq,
        # delta_chisq = comp_vs_bothEq$Chisq_diff,
        # p_value = comp_vs_bothEq$p_value
      )



      # end: if(test == "temporal")
    } else if (test == "contemporaneous"){

      # estimate the model that constrains temporal networks to be equal
      mod_union_tempEq <- mod_union %>%
        groupequal(matrix = "beta") %>% runmodel %>% suppressWarnings

      # estimate the model that constrains contemporaneous networks to be equal
      mod_union_contEq <- mod_union %>%
        groupequal(matrix = "omega_zeta_within") %>% runmodel %>% suppressWarnings

      # compare with the free model
      comp_vs_free <- psychonetrics::compare(free = mod_union,
                                             contEq = mod_union_contEq)

      tab_vs_free <- data.frame(
        AIC = comp_vs_free$AIC,
        delta_AIC = c(NA, diff(comp_vs_free$AIC)),
        BIC = comp_vs_free$BIC,
        delta_BIC = c(NA, diff(comp_vs_free$BIC)),
        row.names = comp_vs_free$model
        # chi_sq = comp_vs_free$Chisq,
        # delta_chisq = comp_vs_free$Chisq_diff,
        # p_value = comp_vs_free$p_value
      )

      # compare with the full-constrained model
      comp_vs_bothEq <- psychonetrics::compare(tempEq = mod_union_tempEq,
                                               bothEq = mod_union_bothEq)

      tab_vs_bothEq <- data.frame(
        AIC = comp_vs_bothEq$AIC,
        delta_AIC = c(NA, diff(comp_vs_bothEq$AIC)),
        BIC = comp_vs_bothEq$BIC,
        delta_BIC = c(NA, diff(comp_vs_bothEq$BIC)),
        row.names = comp_vs_bothEq$model
        # chi_sq = comp_vs_bothEq$Chisq,
        # delta_chisq = comp_vs_bothEq$Chisq_diff,
        # p_value = comp_vs_bothEq$p_value
      )



    # end: if(test == "contemporaneous")
    }  # end: if(test == xxx)


  } # end: else (net_type == "pruned")

  # return the comparison results
  tab <- list(
    vs_free = tab_vs_free,
    vs_bothEq = tab_vs_bothEq
  )
#
#   class(tab_vs_free) <- class(tab_vs_bothEq) <- c("panel_omni", "data.frame")

  # ----- partial pruning -----

  if (partial_prune) {

    # model comparisons
    if(test == "both"){

      mod_pp <- mod_saturated %>%
        partialprune(matrices = "beta", alpha = p_prune_alpha) %>%
        partialprune(matrices = "omega_zeta_within", alpha = p_prune_alpha) %>% runmodel %>% suppressWarnings

      # save networks
      save_matrix = c("PDC", "beta", "omega_zeta_within")

      mat <- lapply(save_matrix, function(m){
        m <- getmatrix(mod_pp, m)
        return(m)
      }) %>% setNames(save_matrix)

      names(mat) <- gsub("PDC", "temporal", names(mat))
      names(mat) <- gsub("omega_zeta_within", "contemporaneous", names(mat))



      # end: if(test == "both")
    } else if (test == "temporal"){

      browser()

      mod_pp <- mod_saturated_contEq %>%
        partialprune(matrices = "beta", alpha = p_prune_alpha) %>% runmodel %>% suppressWarnings

      # save networks
      save_matrix = c("PDC", "beta", "omega_zeta_within")

      mat <- lapply(save_matrix, function(m){
        m <- getmatrix(mod_pp, m)
        return(m)
      }) %>% setNames(save_matrix)

      names(mat) <- gsub("PDC", "temporal", names(mat))
      names(mat) <- gsub("omega_zeta_within", "contemporaneous", names(mat))

      # for the convenience of interpretation, correct the contemp of groups other than g1
      for(g in 2:length(mat$contemporaneous)){
        mat$contemporaneous[[g]] <- mat$contemporaneous[[1]]
      }



      # end: if(test == "temporal")
    } else if (test == "contemporaneous"){

      mod_pp <- mod_saturated_tempEq %>%
        partialprune(matrices = "omega_zeta_within", alpha = p_prune_alpha) %>% runmodel %>% suppressWarnings

      # save networks
      save_matrix = c("PDC", "beta", "omega_zeta_within")

      mat <- lapply(save_matrix, function(m){
        m <- getmatrix(mod_pp, m)
        return(m)
      }) %>% setNames(save_matrix)

      names(mat) <- gsub("PDC", "temporal", names(mat))
      names(mat) <- gsub("omega_zeta_within", "contemporaneous", names(mat))

      # for the convenience of interpretation, correct the temp and beta of groups other than g1
      for(g in 2:length(mat$temporal)){
        mat$temporal[[g]] <- mat$temporal[[1]]
        mat$beta[[g]] <- mat$beta[[1]]
      }

    } # end: if(test == xxx)

  } # end if (partial_prune)
  # warn exploratory pruning
  cat("\nResults of partial pruning are explratory. Be careful to interpret if group-equality constraints decreased AIC or BIC")

  return(list(
    omnibus = tab,
    partial_prune = mat
  ))

} # end: IVPP_panel





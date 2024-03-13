#
# # logic:
# # do the omnibus 1
# # if not significant, stop and return no difference
# # if signigicant, omibus 2
# # Do partial pruning for the significant
#
#
#
# ##=================================
# ##  Model 1: Saturated panel GVAR
# ##=================================
#
# mod_1_saturated <- ml_gvar(data_merged,
#                            vars = varnames,
#                            idvar = "subject",
#                            beepvar = "time",
#                            groups = "group",
#                            standardize = "z",
#                            between = "chol") %>% runmodel
#
#
# # model that constrain temporal networks to be equal across two groups
# mod_1_saturated_tempEq <- ml_gvar(data_merged,
#                                   vars = varnames,
#                                   idvar = "subject",
#                                   beepvar = "time",
#                                   groups = "group",
#                                   standardize = "z",
#                                   between = "chol") %>% groupequal(matrix = "beta") %>% runmodel
#
# # model that constrain contemporaneous networks to be equal across two groups
# mod_1_saturated_contEq <- ml_gvar(data_merged,
#                                   vars = varnames,
#                                   idvar = "subject",
#                                   beepvar = "time",
#                                   groups = "group",
#                                   standardize = "z",
#                                   between = "chol") %>% groupequal(matrix = "omega_zeta_within") %>% runmodel
#
# # model that constrain both networks to be equal across two groups
# mod_1_saturated_bothEq <- ml_gvar(data_merged,
#                                   vars = varnames,
#                                   idvar = "subject",
#                                   beepvar = "time",
#                                   groups = "group",
#                                   standardize = "z",
#                                   between = "chol") %>% groupequal(matrix = "beta") %>% groupequal(matrix = "omega_zeta_within") %>% runmodel
#
#
# ##======================================
# ##  Model 2: Union model (saturated_free)
# ##======================================
#
# # The pruned (sparse) union model
# mod_2_union <- mod_1_saturated %>% prune(alpha = 0.01) %>% unionmodel %>% runmodel
#
# ##==============================================
# ##  Model 3: Constrained models (sparse_equal)
# ##==============================================
#
# mod_3_tempEq = mod_2_union %>% groupequal(matrix = "beta") %>% runmodel
# mod_3_contEq = mod_2_union %>% groupequal(matrix = "omega_zeta_within") %>% runmodel
# mod_3_bothEq = mod_2_union %>% groupequal(matrix = "beta") %>% groupequal(matrix = "omega_zeta_within") %>% runmodel
#
# ##---------------------------------------------------------------
# ##  Model comparison as an omnibus test of network difference  --
# ##---------------------------------------------------------------
#
# ##.............
# ##  saturated
# ##.............
#
# # Check whether group equal constraints in temporal & contemporaneous network lowered model fit
# # these tables always order based on df
# modComp_sat = psychonetrics::compare(saturated_free = mod_1_saturated,
#                                      saturated_tempEq = mod_1_saturated_tempEq,
#                                      saturated_contEq = mod_1_saturated_contEq,
#                                      saturated_bothEq = mod_1_saturated_bothEq) %>% data.table
#
# # Reshape the data.table to have model names as rows
# melted_modComp_sat <- melt(modComp_sat, measure.vars = names(modComp_sat)[-1], id.vars = "model", variable.name = "metric")
#
# # create a new column specifying model & metrics
# melted_modComp_sat[, model_metric := paste(model, metric, sep = "_")]
#
# # cast the table to a single row; remove the first column with only .
# modComp_sat = dcast(melted_modComp_sat, 1 ~ model_metric, value.var = "value")[, -1]
#
# ##..........
# ##  sparse
# ##..........
#
#
# # Check whether group equal constraints in temporal & contemporaneous network lowered model fit
#
# modComp_spa = psychonetrics::compare(sparse_free = mod_2_union,
#                                      sparse_tempEq = mod_3_tempEq,
#                                      sparse_contEq = mod_3_contEq,
#                                      sparse_bothEq = mod_3_bothEq) %>% data.table
#
# # Reshape the data.table to have model names as rows
# melted_modComp_spa <- melt(modComp_spa, measure.vars = names(modComp_spa)[-1], id.vars = "model", variable.name = "metric")
#
# # create a new column specifying model & metrics
# melted_modComp_spa[, model_metric := paste(model, metric, sep = "_")]
#
# # cast the table to a single row; remove the first column with only .
# modComp_spa = dcast(melted_modComp_spa, 1 ~ model_metric, value.var = "value")[, -1]
#
# # merge results of saturated & sparse model
# results = c(modComp_sat, modComp_spa)

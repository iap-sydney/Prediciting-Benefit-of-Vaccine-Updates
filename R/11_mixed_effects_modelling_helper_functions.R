get_mixed_effects_model = function(dataset = regression_data_old_updated, 
                                   outcome = "fold_rise",
                                   include_interaction = F,
                                   include_exposure = T,
                                   include_conflicted_variable = F,
                                   exclude_conflicted_data = F,
                                   use_only_papers_with_old_updated_data = F,
                                   include_interaction_cat = F,
                                   include_exposure_cat = F, 
                                   fixed_comparison_effect_for_post_neuts = T,
                                   include_valency = F) {
  
  model_equation_string = paste0(outcome," ~ 1 +(1 |old_updated_pairing_index) + (1|paper) ",
                                 ifelse(outcome=="fold_rise" | fixed_comparison_effect_for_post_neuts==F," + (1+immunogen_chronology||comparison_no)"," + (-1+immunogen_chronology||comparison_no)"),
                                 case_when(include_interaction&include_exposure&include_exposure_cat ~ " + immunogen_chronology + immunogen_chronology*(prior_doses+has_more_than_2_exposures)",
                                           include_interaction&include_exposure ~ " + immunogen_chronology + immunogen_chronology*prior_doses",
                                           include_interaction_cat&include_exposure_cat ~ " + immunogen_chronology + immunogen_chronology*has_more_than_2_exposures",
                                           T~" + immunogen_chronology"),
                                 ifelse(include_exposure," + prior_doses",""),
                                 if_else(include_conflicted_variable,"+ is_conflicted",""),
                                 if_else(include_exposure_cat,"+ has_more_than_2_exposures",""),
                                 if_else(outcome=="post_neut" & fixed_comparison_effect_for_post_neuts==T,"+ comparison_no",""),
                                 if_else(include_valency,"+startsWith(valency,'bivalent')","")
                                 #if_else(include_valency,"+valency","")
  )
  data_str = paste0("filter(dataset, !is.na(",outcome,")",
                    ifelse(exclude_conflicted_data, ", is_conflicted==FALSE",""),
                    ifelse(use_only_papers_with_old_updated_data, ", old_updated_has_pair==TRUE",""),
                    ")")
  
  
  glmmTMB_str = paste0("model = glmmTMB(",model_equation_string,",",
                       #"control = glmmTMBControl(optimizer = optim, optArgs = list(method=\"BFGS\")),",
                       "data =",data_str,",weights = sqrt(n))")
  
  eval(parse(text = glmmTMB_str))
  model
}


get_mixed_effects_params = function(dataset = regression_data_old_updated, 
                                    outcome = "fold_rise",
                                    include_interaction = F,
                                    include_exposure = T,
                                    include_conflicted_variable = F,
                                    exclude_conflicted_data = F,
                                    use_only_papers_with_old_updated_data = F,
                                    include_interaction_cat = F,
                                    include_exposure_cat = F,
                                    fixed_comparison_effect_for_post_neuts = T, 
                                    include_valency = F) {
  model = get_mixed_effects_model(dataset, 
                                  outcome,
                                  include_interaction,
                                  include_exposure,
                                  include_conflicted_variable,
                                  exclude_conflicted_data,
                                  use_only_papers_with_old_updated_data,
                                  include_interaction_cat,
                                  include_exposure_cat,
                                  fixed_comparison_effect_for_post_neuts, 
                                  include_valency)
  
  coeffs = (summary(model)$coefficients$cond)
  
  model_name = paste0("summary_old_updated_model_",outcome,
                      case_when(include_exposure_cat&include_exposure~"with_exposure_AND_has_more_than_2_exposures",
                                include_exposure ~"_with_exposure",
                                include_exposure_cat~"_with_has_more_than_2_exposures",
                                T~"_no_exposure"),
                      ifelse(include_interaction,"_interaction",""),
                      ifelse(include_conflicted_variable,"_with_conflicted",""),
                      ifelse(exclude_conflicted_data,"_exclude_conflicted_data",""),
                      ifelse(use_only_papers_with_old_updated_data,"_use_only_paired_papers",""),
                      ifelse(include_valency,"_with_valency",""))
  
  
  mixed_eff_params = data.frame(matrix(nrow=1,ncol=length(mixed_eff_params_colnames)))
  model_params = c(round(summary(model)$AIC[1],2),
                   round(c(coeffs["(Intercept)",c(1,4)],
                           coeffs["immunogen_chronologyupdated_immunogen",c(1,4)],
                           if(include_exposure) coeffs["prior_doses",c(1,4)] else c(NA,NA),
                           if(include_interaction&include_exposure){coeffs["immunogen_chronologyupdated_immunogen:prior_doses",c(1,4)]} else {c(NA,NA)},
                           if(include_exposure_cat) coeffs["has_more_than_2_exposuresTRUE",c(1,4)]else c(NA,NA),
                           if(include_interaction_cat) coeffs["immunogen_chronologyupdated_immunogen:has_more_than_2_exposuresTRUE",c(1,4)]else c(NA,NA),
                           if(include_conflicted_variable) coeffs["is_conflictedTRUE",c(1,4)] else c(NA,NA),
                           if(include_valency) coeffs["startsWith(valency, \"bivalent\")TRUE",c(1,4)] else c(NA,NA),
                           coeffs["(Intercept)",2], coeffs["immunogen_chronologyupdated_immunogen",2], 
                           ifelse(include_exposure, coeffs["prior_doses",2],NA),
                           ifelse(include_interaction, coeffs["immunogen_chronologyupdated_immunogen:prior_doses",2],NA), 
                           ifelse(include_interaction_cat, coeffs["immunogen_chronologyupdated_immunogen:has_more_than_2_exposuresTRUE",2],NA),
                           ifelse(include_valency, coeffs["startsWith(valency, \"bivalent\")TRUE",2],NA)),max_dec_pl))
  mixed_eff_params[1,]=c(model_name,outcome,include_exposure,include_exposure_cat,include_interaction,include_conflicted_variable,include_valency,
                         str_detect(model_name,"summary_old"),str_detect(model_name,"updated"),str_detect(model_name,"matched"),
                         model_params)
  colnames(mixed_eff_params) = mixed_eff_params_colnames
  mixed_eff_params
}
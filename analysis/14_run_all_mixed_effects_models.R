outcomes = c("fold_rise","post_neut")

mixed_eff_params = NULL
for(use_outcome in outcomes){
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_interaction = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_exposure_cat = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_exposure_cat = T, include_interaction_cat = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_exposure_cat = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_exposure_cat = T, include_interaction = T))
  
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_interaction = T, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_exposure_cat = T, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_exposure_cat = T, include_interaction_cat = T, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_exposure_cat = T, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_exposure_cat = T, include_interaction = T, include_conflicted_variable = T))
  
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_interaction = T, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_exposure_cat = T, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_exposure_cat = T, include_interaction_cat = T, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_exposure_cat = T, include_conflicted_variable = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_exposure_cat = T, include_interaction = T, include_conflicted_variable = T))
  
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, exclude_conflicted_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, exclude_conflicted_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_interaction = T, exclude_conflicted_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_exposure_cat = T, exclude_conflicted_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_exposure_cat = T, include_interaction_cat = T, exclude_conflicted_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_exposure_cat = T, exclude_conflicted_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_exposure_cat = T, include_interaction = T, exclude_conflicted_data = T))
                           
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, use_only_papers_with_old_updated_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, use_only_papers_with_old_updated_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_interaction = T, use_only_papers_with_old_updated_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_exposure_cat = T, use_only_papers_with_old_updated_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=F, include_exposure_cat = T, include_interaction_cat = T, use_only_papers_with_old_updated_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_exposure_cat = T, use_only_papers_with_old_updated_data = T))
  mixed_eff_params = rbind(mixed_eff_params, get_mixed_effects_params(outcome = use_outcome, include_exposure=T, include_exposure_cat = T, include_interaction = T, use_only_papers_with_old_updated_data = T))
                            
}
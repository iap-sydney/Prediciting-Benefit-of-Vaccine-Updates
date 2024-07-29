generate_dummy_data = function() {
  dummy_data = data.table()
  for (use_interaction in c(TRUE,FALSE)){
    #eval(parse(text=paste0("fold_rise_model = ",filter(mixed_eff_params,dependant_var=="fold_rise", includes_matched==FALSE,with_conflicted==FALSE,with_exposure_cts==TRUE,with_exposure_2_cat==FALSE, includes_interaction==use_interaction)$model_name)))
    fold_rise_model = summary(get_mixed_effects_model(outcome = "fold_rise",include_exposure = T, include_interaction = use_interaction))
    #eval(parse(text=paste0("post_neut_model = ",filter(mixed_eff_params,dependant_var=="post_neut", includes_matched==FALSE,with_conflicted==FALSE,with_exposure_cts==TRUE,with_exposure_2_cat==FALSE, includes_interaction==use_interaction)$model_name)))
    post_neut_model = summary(get_mixed_effects_model(outcome = "post_neut",include_exposure = T, include_interaction = use_interaction))
    fold_coeffs = fold_rise_model$coefficients$cond
    neut_coeffs = post_neut_model$coefficients$cond
    comparisons = rownames(neut_coeffs)[startsWith(rownames(neut_coeffs),"comparison")] 
    names(comparisons) = paste("Comparison",str_remove(comparisons, "comparison_no"))
    if(nrow(neut_coeffs)<7){
      neut_coeffs = rbind(neut_coeffs,rep(0,ncol(neut_coeffs)))
      rownames(neut_coeffs)[7]="immunogen_chronologyupdated_immunogen:prior_doses"
    } 
    neut_coeffs = neut_coeffs[c("(Intercept)","immunogen_chronologyupdated_immunogen","prior_doses",comparisons,"immunogen_chronologyupdated_immunogen:prior_doses"),]
    rownames(neut_coeffs) = c("intercept","updated","prior_exposures",names(comparisons),"interaction")
    
    if(nrow(fold_coeffs)<4){
      fold_coeffs = rbind(fold_coeffs,rep(0,ncol(fold_coeffs)))
      rownames(fold_coeffs)[4]="immunogen_chronologyupdated_immunogen:prior_doses"
    }
    fold_coeffs = fold_coeffs[c("(Intercept)","immunogen_chronologyupdated_immunogen","prior_doses","immunogen_chronologyupdated_immunogen:prior_doses"),]
    rownames(fold_coeffs) = c("intercept","updated","prior_exposures","interaction")
    
    
    dummy_data_temp = as.data.table(expand.grid(
      prior_exposures = c(2,3,4,4.3),
      comparison_no = c(c(1:4),NA), 
      immunogen_chronology = factor(c("old_immunogen","updated_immunogen"), levels = c("old_immunogen","updated_immunogen"))
    ))  %>% mutate(
      post_neut = neut_coeffs["intercept",1]+
        ifelse(immunogen_chronology=="updated_immunogen", neut_coeffs["updated",1]+(prior_exposures-2)*neut_coeffs["interaction",1],0) +
        (prior_exposures-2)*neut_coeffs["prior_exposures",1]+
        case_match(comparison_no, 2~neut_coeffs["Comparison 2",1],
                   3~neut_coeffs["Comparison 3",1],
                   4~neut_coeffs["Comparison 4",1],
                   T~0),
      fold_rise = fold_coeffs["intercept",1]+
        ifelse(immunogen_chronology=="updated_immunogen", fold_coeffs["updated",1]+(prior_exposures-2)*fold_coeffs["interaction",1],0) + 
        (prior_exposures-2)*fold_coeffs["prior_exposures",1],
      comparison_no = factor(ifelse(is.na(comparison_no),"any",comparison_no)), 
      interaction = factor(ifelse(use_interaction, "Yes","No"), levels = c("Yes","No")),
      jitter=0
    ) %>% 
      melt(measure.vars = c("post_neut","fold_rise")) %>%
      mutate(
        variable_label = ifelse(variable=="fold_rise","Fold Rise", "Neutralisation Titres"),
        imm_chron_label = ifelse(immunogen_chronology=="old_immunogen","Old Immunogen","Updated Immunogen")
      )
    dummy_data = rbind(dummy_data,dummy_data_temp)
  }
  dummy_data
}

generate_plot_data_points = function(dataset = all_old_updated_matched_data){
  data_points = dataset %>% 
    melt(measure.vars = c("post_neut","fold_rise")) %>% 
    filter(immunogen_chronology != "matched_immunogen") %>% 
    mutate(interaction = "No",
           variable_label = ifelse(variable=="fold_rise","Fold Rise", "Neutralisation Titres"),
           imm_chron_label = ifelse(immunogen_chronology=="old_immunogen","Old Immunogen","Updated Immunogen")
    )
  data_points$jitter = rnorm(nrow(data_points), sd=small_jitter_width)
  data_points
}

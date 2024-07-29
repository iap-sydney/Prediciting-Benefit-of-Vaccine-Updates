make_study_table = function(dataset = regression_data_old_updated){
  
  # Make the Table with the number of papers etc ----------------------------
  value_table = get_study_table_values(dataset)
  study_table = data.frame(matrix(nrow = 3, ncol = 0))
  for (j in c(1:nrow(value_table))){
    col_vals = NULL
    col_vals = c(paste0(value_table[j,"n_either"]," (",value_table[j,"n_either_fr"],")"),
                 paste0(value_table[j,"n_both"]," (",value_table[j,"n_both_fr"],")"),
                 paste0(value_table[j,"fr_model"]," [",value_table[j,"lo"],"-",value_table[j,"hi"],"]")
    )
    study_table = cbind(study_table,(col_vals))
  }
  colnames(study_table) = rownames(value_table)
  study_table
}

get_study_table_values = function(dataset){
  old_updated_papers =  dataset %>%
    filter(immunogen_chronology %in% c("old_immunogen","updated_immunogen")) %>%
    group_by(comparison_no) %>%
    summarise(n = length(unique(paper2)))
  old_updated_fr_papers = dataset %>%
    filter(immunogen_chronology %in% c("old_immunogen","updated_immunogen"), 
           !is.na(fold_rise)) %>%
    group_by(comparison_no) %>%
    summarise(n = length(unique(paper2)))
  
  old_and_updated = data.frame(matrix(nrow=0,ncol=3))
  for(comp in c(1:4)){
    old_papers = filter(dataset,comparison_no==comp,
                        immunogen_chronology=="old_immunogen")$paper2
    old_papers_fr = filter(dataset,comparison_no==comp,
                           immunogen_chronology=="old_immunogen", !is.na(fold_rise))$paper2
    updated_papers = filter(dataset,comparison_no==comp,
                            immunogen_chronology=="updated_immunogen")$paper2
    updated_papers_fr = filter(dataset,comparison_no==comp,
                               immunogen_chronology=="updated_immunogen", !is.na(fold_rise))$paper2
    old_and_updated_papers = sum(unique(old_papers) %in% unique(updated_papers))
    old_and_updated_papers_fr = sum(unique(old_papers_fr) %in% unique(updated_papers_fr))
    old_and_updated = rbind(old_and_updated,c(comp,old_and_updated_papers,old_and_updated_papers_fr))
  }
  colnames(old_and_updated) = c("comparison_no","n_both","n_both_fr")
  combined_paper_count = merge(old_updated_papers,old_updated_fr_papers, by = "comparison_no")
  colnames(combined_paper_count) = c("comparison_no","n_either","n_either_fr")
  table_vals = (cbind(combined_paper_count,old_and_updated[,c(2,3)]))
  old_papers = filter(dataset,
                      immunogen_chronology=="old_immunogen")$paper2
  old_papers_fr = filter(dataset,
                         immunogen_chronology=="old_immunogen", !is.na(fold_rise))$paper2
  updated_papers = filter(dataset,
                          immunogen_chronology=="updated_immunogen")$paper2
  updated_papers_fr = filter(dataset,
                             immunogen_chronology=="updated_immunogen", !is.na(fold_rise))$paper2
  old_and_updated_papers = sum(unique(old_papers) %in% unique(updated_papers))
  old_and_updated_papers_fr = sum(unique(old_papers_fr) %in% unique(updated_papers_fr))
  
  old_or_updated_papers = unique(filter(dataset, immunogen_chronology %in% c("old_immunogen","updated_immunogen"))$paper2)
  old_or_updated_papers_fr = unique(filter(dataset, immunogen_chronology %in% c("old_immunogen","updated_immunogen"),!is.na(fold_rise))$paper2)
  combined_row = c("combined",
                   num_papers(old_or_updated_papers),
                   num_papers(old_or_updated_papers_fr),
                   old_and_updated_papers,old_and_updated_papers_fr
  )
  # Make the first column not a factor anymore
  table_vals$comparison_no = as.character(table_vals$comparison_no)
  table_vals = rbind(table_vals,combined_row)
  rownames(table_vals) = c(c(1:4),"combined")
  table_vals = table_vals[,c(2:5)]
  
  fold_rise_cols = matrix(nrow=5,ncol = 2)
  for (i in c(1:5)){
    if (i==5){
      compvals = c(1:4)
    } else compvals = i
    this_test = t.test(log10(filter(all_old_updated_matched_data, !is.na(fold_rise), immunogen_chronology == "updated_immunogen", comparison_no %in% compvals )$fold_rise),
                                log10(filter(all_old_updated_matched_data, !is.na(fold_rise), immunogen_chronology == "old_immunogen", comparison_no %in% compvals )$fold_rise))
    
    fold_rise_cols[i,] = c(round(10^(-diff(this_test$estimate)),2),round(this_test$p.value,2))
  }
  fold_rise_model_vals = matrix(nrow=5,ncol = 3)
  
  post_neut_cols = matrix(nrow=5,ncol = 2)
  for (i in c(1:5)){
    if (i==5){
      compvals = c(1:4)
    } else compvals = i
    this_test = t.test(log10(filter(all_old_updated_matched_data, !is.na(post_neut), immunogen_chronology == "updated_immunogen", comparison_no %in% compvals )$post_neut),
                       log10(filter(all_old_updated_matched_data, !is.na(post_neut), immunogen_chronology == "old_immunogen", comparison_no %in% compvals )$post_neut))
    
    post_neut_cols[i,] = c(round(10^(-diff(this_test$estimate)),2),round(this_test$p.value,2))
  }
  colnames(post_neut_cols) = c("pn_data","pval")
  
  simple_model_cols = NULL
  for (i in c(1:5)){
    if (i==5){
      compvals = c(1:4)
    } else compvals = i
    simple_model = glmmTMB(fold_rise ~ 1  + (1|paper)  + immunogen_chronology,
                           data =filter(dataset,comparison_no %in% compvals, !is.na(fold_rise)),weights = sqrt(n))
    val = summary(simple_model)$coefficients$cond[2,1]
    se = summary(simple_model)$coefficients$cond[2,2]
    simple_model_cols = rbind(simple_model_cols,
                              c(round(10^val,2),round(10^(val-1.96*se),2),round(10^(val+1.96*se),2)))
  }
  colnames(simple_model_cols) = c("fr_model","lo","hi")
  table_vals = cbind(table_vals,fold_rise_cols,simple_model_cols,post_neut_cols)
}

num_papers = function(dataset){
  if(is.data.frame(dataset)){
    papers = dataset$paper2
  } else {
    papers = dataset
  }
  len = length(unique(papers))+ifelse("Collier+MillerNEJM" %in% papers,1,0)
  len
}

make_mixed_effects_results_table_for_manuscript = function(dataset = regression_data_old_updated, 
                                 outcome = "fold_rise",
                                 include_interaction = F,
                                 include_exposure = T,
                                 include_conflicted_variable = F,
                                 exclude_conflicted_data = F,
                                 use_only_papers_with_old_updated_data = F,
                                 include_interaction_cat = F,
                                 include_exposure_cat = F,
                                 fixed_comparison_effect_for_post_neuts = T, 
                                 include_valency=F){
  model_params = get_mixed_effects_params(dataset, 
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
  results_table = data.frame(matrix(signif(as.numeric(model_params[c("intercept","intercept_stderr","intercept_pval","updated_benefit_param","updated_benefit_stderr", "updated_benefit_pval",
                                                        "exp_cts_param","exp_cts_stderr","exp_cts_pval")]),3),nrow=3,ncol=3,byrow = T))
  rownames(results_table) = c("Intercept","Immunogen\n(updated vs older)","Exposure number")
  colnames(results_table) = c("Param_value","Std_err","p-value")
  results_table
}


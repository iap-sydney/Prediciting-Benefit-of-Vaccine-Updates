
table1 = make_study_table(regression_data_old_updated)
print(table1)


study_numbers_text = NULL
study_numbers_text = c(study_numbers_text, paste0("We identified ",
                            length(unique(combined_dataset_to_use$paper2)),
                            " studies[",ref_start,"-",ref_start+length(unique(combined_dataset_to_use$paper2))-1,
                            "] in which neutralising antibody titres were measured"))
num_immunogens_per_paper = summarise(combined_dataset_to_use, n_imms=length(unique(immunogen_details)), .by = paper2)
study_numbers_text = c(study_numbers_text,paste0("Of these, ",
                            length(unique(filter(combined_dataset_to_use,!is.na(fold_rise))$paper2)),
                            " studies reported pre-neutralisation titres (though in ",
                            length(unique(filter(combined_dataset_to_use,!is.na(fold_rise))$paper2)[unique(filter(combined_dataset_to_use,!is.na(fold_rise))$paper2) %in% unique(filter(combined_dataset_to_use,is.na(fold_rise))$paper2)]),
                            " of these studies pre-boost titres were not reported for all immunogen/variant combinations) and ",
                            length(unique(combined_dataset_to_use$paper2))-length(unique(filter(combined_dataset_to_use,!is.na(fold_rise))$paper2)),
                            " studies only reported post-boost titres. ",
                            length(filter(num_immunogens_per_paper,n_imms==1)$paper2),
                            " studies analysed responses to a single relevant immunogen-containing booster, ",
                            length(filter(num_immunogens_per_paper,n_imms==2)$paper2)," compared neutralising antibodies after two different immunogen-containing boosters, and ",
                            length(filter(num_immunogens_per_paper,n_imms>=3)$paper2)," compared responses after 3 or more immunogen-containing boosters.")
)
num_comparisons_per_paper = summarise(regression_data_old_updated, n_comps=length(unique(comparison_no)), .by = paper2)
num_comparisons_per_paper_per_immnunogen_chron = regression_data_old_updated %>%
  summarise(n_comps=length(unique(comparison_no)), .by = c(paper2,immunogen_chronology)) %>%
  summarise(.by = paper2, n_comps = min(n_comps))
study_numbers_text = c(study_numbers_text,paste0("Of the ",
                                                 length(unique(combined_dataset_to_use$paper2))," studies we identified, only ",
                                                 length(filter(num_comparisons_per_paper, n_comps==4)$paper2)," studies[X, Y] contained data that could be used to contribute to all four comparisons, however ",
                                                 length(filter(num_comparisons_per_paper_per_immnunogen_chron, n_comps==4)$paper2)," study contained data for both the old and updated immunogens across all four comparisons. In addition, ",
                                                 length(unique(filter(combined_dataset_to_use, !(paper2 %in% regression_data_old_updated$paper2))$paper2))," study could not be used in our analysis as it did not include titres after one of the relevant immunogens[Z], leaving ",
                                                 length(unique(combined_dataset_to_use$paper2)) - length(unique(filter(combined_dataset_to_use, !(paper2 %in% regression_data_old_updated$paper2))$paper2)),
                                                 " studies that could contribute to one or more of the comparisons identified in Figure 2A.")
)
print(study_numbers_text)

model_numbers_text = NULL

simple_model = glmmTMB(fold_rise ~ 1  + (1|paper)  + immunogen_chronology,
                       data =filter(regression_data_old_updated, comparison_no %in% c(1:4), !is.na(fold_rise)),weights = sqrt(n))
val = summary(simple_model)$coefficients$cond["immunogen_chronologyupdated_immunogen","Estimate"]
se = summary(simple_model)$coefficients$cond["immunogen_chronologyupdated_immunogen","Std. Error"]
se_mults = c(0,-1.96,1.96)
simple_benefit_ests = round(10^(val+se_mults*se),2)
full_model_results = get_mixed_effects_params(regression_data_old_updated,"fold_rise",F,T)
full_model_exposure_vals = as.numeric(c(full_model_results$exp_cts_param,full_model_results$exp_cts_stderr))
full_model_exposure_ests = round(10^(-full_model_exposure_vals[1]+full_model_exposure_vals[2]*se_mults),1)
full_model_benefit_vals = as.numeric(c(full_model_results$updated_benefit_param,full_model_results$updated_benefit_stderr, full_model_results$updated_benefit_pval))
full_model_benefit_ests = round(10^(full_model_benefit_vals[1]+full_model_benefit_vals[2]*se_mults),2)

regression_data_old_updated = mutate(regression_data_old_updated,
                                     cohort = paste0(paper,group, prior_doses,prior_status, immunogen_group)
                                     )
model_numbers_text = c(model_numbers_text, 
                       paste0("Considering all the data identified from studies that could contribute to the comparisons (",
                              length(unique(filter(regression_data_old_updated, !is.na(fold_rise))$cohort))," cohorts from ",
                              length(unique(filter(regression_data_old_updated,!is.na(fold_rise))$paper2)),
                              " studies), we found that the neutralisation titres against the specified future variant were, on average, ",
                              simple_benefit_ests[1],"-fold [95% CI ",simple_benefit_ests[2],"-",simple_benefit_ests[3],
                              "] greater after boosting with the updated immunogen compared to the older immunogen."),
                       paste0("the overall level of boosting declined with an increasing number of prior exposures (by ",
                              full_model_exposure_ests[1],"-fold per additional exposure [95%CI ",
                              full_model_exposure_ests[2],"-",
                              full_model_exposure_ests[3],") see Table 2). "),
                       paste0("However, even after accounting for the number of prior exposures, boosting with an updated immunogen was estimated to give a ",
                              full_model_benefit_ests[1],"-fold better boost than boosting with an older immunogen [95%CI = ",
                              full_model_benefit_ests[2],"-",full_model_benefit_ests[3],
                              ", p=",round(full_model_benefit_vals[3],2),
                              "). That is, on average, neutralisation titres to the future variant were boosted ",
                              100*(full_model_benefit_ests[1]-1),"% higher after an updated immunogen boost, than after an older immunogen boost. ")

)
print(model_numbers_text)

sensitivity_text = NULL
sensitivity_model_uninfected = get_mixed_effects_params(filter(regression_data_old_updated,prior_status == "uninfected"),"fold_rise",F,T)
sensitivity_model_uninfected_benefit = as.numeric(c(sensitivity_model_uninfected$updated_benefit_param,sensitivity_model_uninfected$updated_benefit_stderr))
sensitivity_model_uninfected_benefit_ests = round(10^(sensitivity_model_uninfected_benefit[1]+sensitivity_model_uninfected_benefit[2]*se_mults),2)

sensitivity_model_homogeneous = get_mixed_effects_params(filter(regression_data_old_updated,prior_status != "mixed"),"fold_rise",F,T)
sensitivity_model_homogeneous_benefit = as.numeric(c(sensitivity_model_homogeneous$updated_benefit_param,sensitivity_model_homogeneous$updated_benefit_stderr))
sensitivity_model_homogeneous_benefit_ests = round(10^(sensitivity_model_homogeneous_benefit[1]+sensitivity_model_homogeneous_benefit[2]*se_mults),2)

sensitivity_model_noVoC = get_mixed_effects_params(filter(regression_data_old_updated,comparison_no %in% c(2,4)),"fold_rise",F,T)
sensitivity_model_noVoC_benefit = as.numeric(c(sensitivity_model_noVoC$updated_benefit_param,sensitivity_model_noVoC$updated_benefit_stderr))
sensitivity_model_noVoC_benefit_ests = round(10^(sensitivity_model_noVoC_benefit[1]+sensitivity_model_noVoC_benefit[2]*se_mults),2)

sensitivity_model_paired = get_mixed_effects_params(filter(regression_data_old_updated,old_updated_has_pair==TRUE),"fold_rise",F,T)
sensitivity_model_paired_benefit = as.numeric(c(sensitivity_model_paired$updated_benefit_param,sensitivity_model_paired$updated_benefit_stderr))
sensitivity_model_paired_benefit_ests = round(10^(sensitivity_model_paired_benefit[1]+sensitivity_model_paired_benefit[2]*se_mults),2)

sensitivity_model_conflicted = get_mixed_effects_params(filter(regression_data_old_updated,is_conflicted==FALSE),"fold_rise",F,T)
sensitivity_model_conflicted_benefit = as.numeric(c(sensitivity_model_conflicted$updated_benefit_param,sensitivity_model_conflicted$updated_benefit_stderr))
sensitivity_model_conflicted_benefit_ests = round(10^(sensitivity_model_conflicted_benefit[1]+sensitivity_model_conflicted_benefit[2]*se_mults),2)

sensitivity_model_deployed = get_mixed_effects_params(filter(regression_data_old_updated,comparison_no %in% c(2,4), immunogen_details %in% c("Ancestral","Ancestral_BA.5","Ancestral_BA.1","Ancestral_BA1","Ancestral_BA5")),"fold_rise",F,T)
sensitivity_model_deployed_benefit = as.numeric(c(sensitivity_model_deployed$updated_benefit_param,sensitivity_model_deployed$updated_benefit_stderr))
sensitivity_model_deployed_benefit_ests = round(10^(sensitivity_model_deployed_benefit[1]+sensitivity_model_deployed_benefit[2]*se_mults),2)


sensitivity_text = c(sensitivity_text,
                     paste0("We found that in both cases an updated immunogen also provided a greater boost to neutralisation titres (",
                            sensitivity_model_uninfected_benefit_ests[1]," [95%CI ",
                            sensitivity_model_uninfected_benefit_ests[2],"-",sensitivity_model_uninfected_benefit_ests[3],"] and ",
                            sensitivity_model_homogeneous_benefit_ests[1]," [95%CI ",
                            sensitivity_model_homogeneous_benefit_ests[2],"-",sensitivity_model_homogeneous_benefit_ests[3],"] respectively). "),
                     paste0("We found an almost identical benefit from the use of updated booster immunogen when excluding any comparisons that involved the early VoC immunogens (",
                            sensitivity_model_noVoC_benefit_ests[1]," [95%CI ",
                            sensitivity_model_noVoC_benefit_ests[2],"-",sensitivity_model_noVoC_benefit_ests[3],
                            "] benefit of the updated immunogen"),
                     paste0("Including only these data, we found a similar benefit from the use of an updated booster immunogen (",
                            sensitivity_model_deployed_benefit_ests[1],"-fold [95%CI ",
                            sensitivity_model_deployed_benefit_ests[2],"-",sensitivity_model_deployed_benefit_ests[3],
                            "] benefit of the updated immunogen,"),
                     paste0("For example, ",
                            num_papers(filter(regression_data_old_updated, is_conflicted==T))," of the ",
                            num_papers(regression_data_old_updated)," studies (",
                            round(100*num_papers(filter(regression_data_old_updated, is_conflicted==T))/num_papers(regression_data_old_updated),0),"%), that we used in our analysis were sponsored by the companies producing the vaccines. We therefore performed analysis of the subset of ",
                            num_papers(filter(regression_data_old_updated,is_conflicted==FALSE)$paper)," studies with no pharmaceutical company involvement and found a similar benefit of updating the booster immunogen to that described above (update advantage of ",
                            sensitivity_model_conflicted_benefit_ests[1]," [95%CI ",
                            sensitivity_model_conflicted_benefit_ests[2],"-",sensitivity_model_conflicted_benefit_ests[3],
                            "], "))
print(sensitivity_text)

absolute_text = NULL
absolute_model = get_mixed_effects_params(regression_data_old_updated,"post_neut",F,T, fixed_comparison_effect_for_post_neuts=F)
absolute_model_benefit = as.numeric(c(absolute_model$updated_benefit_param,absolute_model$updated_benefit_stderr))
absolute_model_benefit_ests = round(10^(absolute_model_benefit[1]+absolute_model_benefit[2]*se_mults),2)

absolute_text = c(absolute_text,
                     paste0("However, under ",
                            round(100*length(unique(filter(regression_data_old_updated, !is.na(fold_rise))$paper2))/length(unique(filter(regression_data_old_updated)$paper2)),1),
                            "% (",
                            length(unique(filter(regression_data_old_updated, !is.na(fold_rise))$paper2)),"/",
                            length(unique(filter(regression_data_old_updated)$paper2)),") of the studies we identified reported both pre- and post-boost titres "),
                     paste0("We found that the absolute neutralisation titres after boosting with an updated immunogen were estimated to be ",
                            absolute_model_benefit_ests[1]," [95%CI ",
                            absolute_model_benefit_ests[2],"-",absolute_model_benefit_ests[3],"] higher than after boosting with an older immunogen"))
print(absolute_text)

table2 = make_mixed_effects_results_table_for_manuscript(regression_data_old_updated, outcome = "fold_rise")
print(table2)
table2_text = NULL
for (i in c(1:ncol(table2))){
  table2_text[i] = paste0(round(table2$Param_value[i],3),
                          " (",round(table2$Param_value[i]+table2$Std_err[i]*se_mults[2],3),"-",
                          round(table2$Param_value[i]+table2$Std_err[i]*se_mults[3],3),"), p=",
                          table2$`p-value`[i])
}
print(table2_text)

or_text = NULL
ORbenefit = round(getOR("severe",benefit = 1.4, returnCIs = T),2)
ORbenefit2 = round(getOR("severe",benefit = 1.4, returnCIs = T)^2,2)
or_text = c(or_text,
                  paste0("If an updated booster provides on average a 40% boost in titre compared to an older immunogen, then the odds ratio of protection from COVID-19 in a homogeneous population is ",
                         ORbenefit[1]," [95% CI ",ORbenefit[2],"-",ORbenefit[3],"];"),
                  paste0("compared to one receiving a booster that is two immunogens behind, is predicted to be ",
                         ORbenefit[1],"^2=",ORbenefit2[1]," [95% CI ",ORbenefit2[2],"-",ORbenefit2[3],"]. ")
                  )
print(or_text)

clinical_text = NULL
natmedparams = get_natmed2021_params()
benefit = 1.4
base_protect = .6
min_protect = .3
max_protect = .8
base_neut = get_neut_from_efficacy(base_protect)
base_severe = LogisticModel_PercentUninfected(log10(base_neut),natmedparams$sig, natmedparams$logk_Severe, natmedparams$C50_Severe)
updated_protect = LogisticModel_PercentUninfected(log10(base_neut*benefit),natmedparams$sig, natmedparams$logk, natmedparams$C50)
updated_severe = LogisticModel_PercentUninfected(log10(base_neut*benefit),natmedparams$sig, natmedparams$logk_Severe, natmedparams$C50_Severe)

min_neut = get_neut_from_efficacy(min_protect)
min_severe = LogisticModel_PercentUninfected(log10(min_neut),natmedparams$sig, natmedparams$logk_Severe, natmedparams$C50_Severe)
updated_min_protect = LogisticModel_PercentUninfected(log10(min_neut*benefit),natmedparams$sig, natmedparams$logk, natmedparams$C50)
updated_min_severe = LogisticModel_PercentUninfected(log10(min_neut*benefit),natmedparams$sig, natmedparams$logk_Severe, natmedparams$C50_Severe)

max_neut = get_neut_from_efficacy(max_protect)
max_severe = LogisticModel_PercentUninfected(log10(max_neut),natmedparams$sig, natmedparams$logk_Severe, natmedparams$C50_Severe)
updated_max_protect = LogisticModel_PercentUninfected(log10(max_neut*benefit),natmedparams$sig, natmedparams$logk, natmedparams$C50)
updated_max_severe = LogisticModel_PercentUninfected(log10(max_neut*benefit),natmedparams$sig, natmedparams$logk_Severe, natmedparams$C50_Severe)


clinical_text = c(clinical_text,
                  paste0("If an older booster vaccine gave rise to ",
                     round(100*base_protect,0),"% protection in the population against symptomatic disease and ",
                     round(100*base_severe,0),"% protection against severe disease, then the updated booster would be predicted to increase protection to ",
                     round(100*updated_protect,0),"% against symptomatic and ",
                     round(100*updated_severe,0),"% against severe disease. If a clinical trial were to be run, we would expect to observe a relative efficacy of ",
                     round(100*(1-(1-updated_protect)/(1-base_protect)),0),"% and ",
                     round(100*(1-(1-updated_severe)/(1-base_severe)),0),"% for the updated immunogen compared to the historical immunogen against symptomatic and severe disease respectively."),
                paste0("More generally, if we assume that an older booster offers between ",
                     round(100*min_protect,0),"% and ",
                     round(100*max_protect,0),"% protection against symptomatic disease. We can then compare this to the protection provided by an updated boosted immunogen that gives a ",
                     benefit,"-fold improvement in neutralisation titres. In this case, we would observe an ",
                     round(100*(1-(1-updated_min_protect)/(1-min_protect)),0),"-",round(100*(1-(1-updated_max_protect)/(1-max_protect)),0),
                     "% effectiveness of the updated booster vaccine against symptomatic disease and a ",
                     round(100*(1-(1-updated_min_severe)/(1-min_severe)),0),"-",round(100*(1-(1-updated_max_severe)/(1-max_severe)),0),
                     "% effectiveness against severe disease compared to the older vaccine immunogen"))
print(clinical_text)

discussion_text = NULL
conflicted_model = get_mixed_effects_params(regression_data_old_updated,"fold_rise",
                                            include_interaction=F,include_exposure=T, exclude_conflicted_data = T)
conflicted_model_benefit = as.numeric(c(conflicted_model$updated_benefit_param,conflicted_model$updated_benefit_stderr))
conflicted_model_benefit_ests = round(10^(conflicted_model_benefit[1]+conflicted_model_benefit[2]*se_mults),2)
conflicted_model2 = get_mixed_effects_params(regression_data_old_updated,"fold_rise",
                                            include_interaction=F,include_exposure=T, include_conflicted_variable = T)

discussion_text = paste0(discussion_text, 
                         paste0("(i.e. paired comparison in the same study), we found very similar results (update advantage of ",
                                sensitivity_model_paired_benefit_ests[1]," [95%CI ",
                                sensitivity_model_paired_benefit_ests[2],"-",sensitivity_model_paired_benefit_ests[3],
                                "], "),
                         paste0("We therefore performed analysis of the subset of ",
                                num_papers(filter(regression_data_old_updated, is_conflicted=="FALSE", !is.na(fold_rise))),
                                " studies with pre-boost neutralisation titres and no pharmaceutical company involvement and found a similar benefit of updating the booster immunogen to that described above (update advantage of ",
                                  conflicted_model_benefit_ests[1],"-fold [95%CI ",
                                conflicted_model_benefit_ests[2],"-",conflicted_model_benefit_ests[3],"], "),
                         paste0("). Additionally, in our multiple regression analysis when we included pharmaceutical involvement as a potential factor, we found this variable was ",
                                ifelse(as.numeric(conflicted_model2$conflicted_pval)>.05,"not ",""),"significant.")
                         )
print(discussion_text)


supp_table_2 = make_mixed_effects_results_table_for_manuscript(dataset = filter(regression_data_old_updated, prior_status=="uninfected"),
                                                               outcome = "fold_rise",
                                                               include_interaction = F,
                                                               include_exposure = T)
print(supp_table_2)

supp_table_3 = make_mixed_effects_results_table_for_manuscript(dataset = filter(regression_data_old_updated, prior_status!="mixed"),
                                                               outcome = "fold_rise",
                                                               include_interaction = F,
                                                               include_exposure = T)
print(supp_table_3)

supp_table_4 = make_mixed_effects_results_table_for_manuscript(dataset = filter(regression_data_old_updated, comparison_no %in% c(2,4)),
                                                               outcome = "fold_rise",
                                                               include_interaction = F,
                                                               include_exposure = T)
print(supp_table_4)

supp_table_5 = make_mixed_effects_results_table_for_manuscript(dataset = regression_data_old_updated,
                                                               outcome = "post_neut",
                                                               include_interaction = F,
                                                               include_exposure = T,
                                                               fixed_comparison_effect_for_post_neuts = F)
print(supp_table_5)

supp_table_6 = make_mixed_effects_results_table_for_manuscript(dataset = filter(regression_data_old_updated, old_updated_has_pair == TRUE),
                                                               outcome = "fold_rise",
                                                               include_interaction = F,
                                                               include_exposure = T)
print(supp_table_6)

supp_table_7 = make_mixed_effects_results_table_for_manuscript(dataset = filter(regression_data_old_updated, is_conflicted == "FALSE"),
                                                               outcome = "fold_rise",
                                                               include_interaction = F,
                                                               include_exposure = T)
print(supp_table_7)
supp_table_7a = make_mixed_effects_results_table_for_manuscript(dataset = regression_data_old_updated,
                                                               outcome = "fold_rise",
                                                               include_interaction = F,
                                                               include_exposure = T,
                                                               exclude_conflicted_data=T)
print(supp_table_7a==supp_table_7)

supp_table_8 = make_mixed_effects_results_table_for_manuscript(dataset = filter(regression_data_old_updated, comparison_no %in% c(2,4),immunogen_details %in% c("Ancestral","Ancestral_BA.5","Ancestral_BA.1","Ancestral_BA1","Ancestral_BA5")),
                                                               outcome = "fold_rise",
                                                               include_interaction = F,
                                                               include_exposure = T)
print(supp_table_8)

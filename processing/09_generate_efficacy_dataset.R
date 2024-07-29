## Efficacy dataset
s_small = .0001
s=.01
baseline_eff = c(seq(s_small,s-s_small,s_small),seq(s,1-s,s))

rises_by_dose = c(13.401453,6.782245,3.432378)
names(rises_by_dose) = paste0(c(2:4),"doses")
updated_boost = 1.465228

eff_cols = c("baseline_eff","baseline_neut","sympt_baseline",paste0("sympt_",names(rises_by_dose)),paste0("sympt_updated_",names(rises_by_dose)),
             "severe_baseline",paste0("severe_",names(rises_by_dose)),paste0("severe_updated_",names(rises_by_dose)),
             paste0("sympt_",names(rises_by_dose),"_increase_over_baseline"),
             paste0("severe_",names(rises_by_dose),"_increase_over_baseline"),
             paste0("sympt_",names(rises_by_dose),"_increase_over_baseline_with_updated"),
             paste0("severe_",names(rises_by_dose),"_increase_over_baseline_with_updated"),
             paste0("sympt_",names(rises_by_dose),"_increase_with_updated"),
             paste0("severe_",names(rises_by_dose),"_increase_with_updated")
             )
effs = data.frame(matrix(nrow = 0, ncol = length(eff_cols))) 
# Assign column names
colnames(effs) = eff_cols
baseline_neut = c()
for (i in c(1:length(baseline_eff))){
  this_neut = get_neut_from_efficacy(baseline_eff[i])
  sympt_eff = baseline_eff[i]
  sympt_eff_by_dose = get_efficacy_from_neut(this_neut*rises_by_dose)
  sympt_eff_by_dose_updated = get_efficacy_from_neut(this_neut*updated_boost*rises_by_dose)
  severe_eff = get_efficacy_from_neut(this_neut,4)
  severe_eff_by_dose = get_efficacy_from_neut(this_neut*rises_by_dose,4)
  severe_eff_by_dose_updated = get_efficacy_from_neut(this_neut*updated_boost*rises_by_dose,4)
  this_row = c(baseline_eff[i],this_neut, sympt_eff,
              sympt_eff_by_dose,sympt_eff_by_dose_updated,
              severe_eff,
              severe_eff_by_dose,severe_eff_by_dose_updated,
              sympt_eff_by_dose-sympt_eff,
              severe_eff_by_dose-severe_eff,
              sympt_eff_by_dose_updated-sympt_eff,
              severe_eff_by_dose_updated-severe_eff,
              sympt_eff_by_dose_updated-sympt_eff_by_dose,
              severe_eff_by_dose_updated - severe_eff_by_dose)
  effs[i,]=this_row
}

eff_df  = effs %>%
  as.data.table() %>%
  melt(id.vars = c("baseline_eff","baseline_neut")) %>%
  mutate(outcome = factor(case_when(str_detect(variable,"severe")~"Severe",
                                    str_detect(variable,"sympt")~"Symptomatic",
                                    T~"NA"),
                          levels = c("Symptomatic","Severe","NA")),
         data_type = factor(case_when(
           str_detect(variable,"increase_with_updated")~"Difference from Older Immunogen",
           str_detect(variable,"increase")~"Difference from Baseline",
          T~"Efficacy"), levels = c("Efficacy","Difference from Older Immunogen","Difference from Baseline")),
         immunogen_chronology = factor(case_when(str_detect(variable,"updated")~"updated_immunogen",
                                                 T~"old_immunogen"), levels = immunogen_chronology_levels),
         #comparison = factor(case_when(str_detect(variable,"increase_with_updated")~"cf older_booster",
         #                              data_type=="Difference"~"",
         #                                T~""), levels = c("cf older_booster","")),
         dosesStart = str_locate(variable,"doses")[,1]-1,
         doses = as.numeric(str_sub(variable,dosesStart,dosesStart)),
         doses = ifelse(is.na(doses),0,doses)
)

# updated vs old efficacy -------------------------------------------------

updated_boost = 1.4

eff_old_updated_cols = c("old_boost_eff","old_boost_neut","sympt_old_boost","sympt_updated_boost",
             "severe_old_boost","severe_updated_boost",
             "sympt_increase_over_old",
             "severe_increase_over_old","rel_eff_sympt",
             "rel_eff_severe")

effs_old_updated = data.frame(matrix(nrow = 0, ncol = length(eff_old_updated_cols))) 
# Assign column names
colnames(effs_old_updated) = eff_old_updated_cols
baseline_neut = c()
for (i in c(1:length(baseline_eff))){
  this_neut = get_neut_from_efficacy(baseline_eff[i])
  sympt_eff_old = baseline_eff[i]
  sympt_eff_updated = get_efficacy_from_neut(this_neut*updated_boost)
  severe_eff_old = get_efficacy_from_neut(this_neut,4)
  severe_eff_updated = get_efficacy_from_neut(this_neut*updated_boost,4)
  rel_eff_sympt = (1-(1-sympt_eff_updated)/(1-sympt_eff_old))
  rel_eff_severe = (1-(1-severe_eff_updated)/(1-severe_eff_old))
  this_row = c(baseline_eff[i],this_neut, sympt_eff_old,
               sympt_eff_updated,
               severe_eff_old,
               severe_eff_updated,
               sympt_eff_updated-sympt_eff_old,
               severe_eff_updated-severe_eff_old,
               rel_eff_sympt,
               rel_eff_severe
               )
  effs_old_updated[i,]=this_row
}

effs_old_updated_df  = effs_old_updated %>%
  as.data.table() %>%
  melt(id.vars = c("old_boost_eff","old_boost_neut")) %>%
  mutate(outcome = factor(case_when(str_detect(variable,"severe")~"Severe",
                                    str_detect(variable,"sympt")~"Symptomatic",
                                    T~"NA"),
                          levels = c("Symptomatic","Severe","NA")),
         data_type = factor(case_when(
           str_detect(variable,"increase_over_old")~eff_plot_factors[2],
           str_detect(variable,"rel_eff")~eff_plot_factors[3],
           T~eff_plot_factors[1]), levels = eff_plot_factors),
         immunogen_chronology = factor(case_when(str_detect(variable,"updated")~"updated_immunogen",
                                                 str_detect(variable,"increase_over_old")~"updated_immunogen",
                                                 str_detect(variable,"rel_eff")~"updated_immunogen",
                                                 T~"old_immunogen"), levels = immunogen_chronology_levels),
         letter_label = "A"
         #comparison = factor(case_when(str_detect(variable,"increase_with_updated")~"cf older_booster",
         #                              data_type=="Difference"~"",
         #                                T~""), levels = c("cf older_booster","")),

  )
effs_old_updated_df$letter_label
for (i in c(1:nrow(effs_old_updated_df))){
  effs_old_updated_df$letter_label[i] = intToUtf8(2*as.numeric(effs_old_updated_df$data_type[i])+(as.numeric(effs_old_updated_df$outcome[i])-1)+63)
}

rm(effs, baseline_eff, baseline_neut, s, effs_old_updated)

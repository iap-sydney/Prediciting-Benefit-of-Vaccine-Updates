
# combine_data_by_immunogen = function(immunogen_data_A, immunogen_data_B, matched=TRUE, suffixes = c(".old",".new"), comparison_table, r, combined_immunogen_data = NULL){
old_updated_data_paired_temp = data.frame()
old_updated_data_all_temp = data.frame()
updated_matched_data_paired_temp = data.frame()
updated_matched_data_all_temp = data.frame()

# Decide which dataset to use, and if we are calling functions in their orginal format ("orig") or in the cleaned data format ("cleaned")
#combined_dataset_to_use = combined_reduced_human_data
#called_from = "orig" 
combined_dataset_to_use = filter(updated_immunogen_neut_data, use_in_analysis)
called_from = "cleaned"

for (r in c(1:nrow(comparison_table))){
  # Get the old, updated and matched data for this comparison
  this_old_immunogen_data = filter(combined_dataset_to_use,
                                   immunogen_group == comparison_table[r,"old_immunogen"],
                                   variant_group ==comparison_table[r,"future_variant"])
  this_updated_immunogen_data = filter(combined_dataset_to_use,
                                       immunogen_group == comparison_table[r,"updated_immunogen"],
                                       variant_group ==comparison_table[r,"future_variant"])
  this_matched_immunogen_data = filter(combined_dataset_to_use,
                                       immunogen_group == comparison_table[r,"future_variant"],
                                       variant_group ==comparison_table[r,"future_variant"])
  #print(paste("Comparison",r,"old-updated paired"))
  old_updated_data_paired_temp = combine_data_by_immunogen(this_old_immunogen_data, this_updated_immunogen_data,suffixes = c(".old",".new"),
                                                           comparison_table, r, paired=TRUE, old_updated_data_paired_temp, called_from)
  #print(paste("Comparison",r,"old-updated all"))
  old_updated_data_all_temp = combine_data_by_immunogen(this_old_immunogen_data, this_updated_immunogen_data,suffixes = c(".old",".updated"),
                                                        comparison_table, r, paired=FALSE, old_updated_data_all_temp, called_from)
  #print(paste("Comparison",r,"updated-matched paired"))
  updated_matched_data_paired_temp = combine_data_by_immunogen(this_updated_immunogen_data, this_matched_immunogen_data,suffixes = c(".updated",".matched"),
                                                               comparison_table, r, paired=TRUE, updated_matched_data_paired_temp, called_from)
  #print(paste("Comparison",r,"updated-matched all"))
  updated_matched_data_all_temp = combine_data_by_immunogen(this_updated_immunogen_data, this_matched_immunogen_data,suffixes = c(".updated",".matched"),
                                                            comparison_table, r, paired=FALSE, updated_matched_data_all_temp, called_from)
}
old_updated_matched_data_all_temp = distinct(rbind(old_updated_data_all_temp,updated_matched_data_all_temp))

old_updated_data_all_temp2 = combine_with_ancestral_data(old_updated_data_all_temp, called_from = called_from)
old_updated_data_paired_temp2 = combine_with_ancestral_data(old_updated_data_paired_temp, called_from = called_from)
updated_matched_data_all_temp2 = combine_with_ancestral_data(updated_matched_data_all_temp, called_from = called_from)
updated_matched_data_paired_temp2 = combine_with_ancestral_data(updated_matched_data_paired_temp, called_from = called_from)
old_updated_matched_data_all_temp2 = combine_with_ancestral_data(old_updated_matched_data_all_temp, called_from = called_from)

select_cols_join = c("paper","paper2","pairing_index","dataset", "group","another_group","prior_status", "prior_doses","old_immunogen","updated_immunogen","future_variant","comparison_name","comparison_no","variant","variant_group")
select_cols_vars = c("crhd_index","crhd_and_comparison_index","immunogen_chronology","n","fold_rise","post_neut", "prior_doses","propn_prior_infected","post_neut.ancestral", "post_neut_normalised","immunogen_group","immunogen_details","valency","is_conflicted", "use_in_analysis")

old_updated_data_paired = old_updated_data_paired_temp2 %>%
  select(any_of(c(select_cols_join,select_cols_vars))) %>%
  mutate(
    post_neut_normalised = post_neut/post_neut.ancestral,
    old_updated_has_pair = T,
    old_updated_pairing_index = pairing_index) %>%
  select(-pairing_index) %>%
  distinct()

old_updated_data_unpaired = old_updated_data_all_temp2 %>%
  select(any_of(c(select_cols_join,select_cols_vars))) %>%
  mutate(
    post_neut_normalised = post_neut/post_neut.ancestral,
    old_updated_has_pair = F,
    old_updated_pairing_index = 100+pairing_index) %>%
  select(-pairing_index) %>%
  filter(!(crhd_and_comparison_index %in% old_updated_data_paired$crhd_and_comparison_index)) %>%
  distinct() 

old_updated_data_all = bind_rows(old_updated_data_paired,old_updated_data_unpaired) %>%
  distinct()

updated_matched_data_paired = updated_matched_data_paired_temp2 %>%
  select(any_of(c(select_cols_join,select_cols_vars))) %>%
  mutate(
    post_neut_normalised = post_neut/post_neut.ancestral,
    updated_matched_has_pair = T,
    updated_matched_pairing_index = pairing_index) %>%
  select(-pairing_index) %>%
  distinct()

updated_matched_data_unpaired = updated_matched_data_all_temp2 %>%
  select(any_of(c(select_cols_join,select_cols_vars))) %>%
  mutate(
    post_neut_normalised = post_neut/post_neut.ancestral,
    updated_matched_has_pair = F,
    updated_matched_pairing_index = 100+pairing_index) %>%
  select(-pairing_index) %>%
  filter(!(crhd_and_comparison_index %in% updated_matched_data_paired$crhd_and_comparison_index)) %>%
  distinct()

updated_matched_data_all = bind_rows(updated_matched_data_paired,updated_matched_data_unpaired) %>%
  distinct()


all_old_updated_matched_data = bind_rows(old_updated_data_all,updated_matched_data_all) %>%
  arrange(crhd_and_comparison_index) %>%
  distinct(pick(!(matches("old_updated") | matches("updated_matched")))) %>% # make sure the updated immunogens are distinct rows
  left_join(old_updated_data_all, by = head(colnames(old_updated_data_all),-2)) %>%
  left_join(updated_matched_data_all, by = head(colnames(updated_matched_data_all),-2)) %>%
  rename(variant_details = variant) %>%
  relocate(comparison_no, .before=old_immunogen) %>%
  relocate(variant_details, .after=variant_group) %>%
  relocate(starts_with("crhd"), .before=comparison_no) %>%
  relocate(starts_with("immunogen"), .before=variant_group) %>%
  relocate(use_in_analysis, .after=last_col()) %>%
  relocate(propn_prior_infected, .after=last_col()) %>%
  filter(use_in_analysis == TRUE) %>%
  mutate(paper = factor(paper, levels = paper_list_cleaned),
         #prior_exposures = prior_doses+ case_match(prior_status, "infected"~1,"uninfected"~0, "mixed"~.5)
         prior_exposures = prior_doses+propn_prior_infected
  )

rm(list=grep("select_cols",ls(),value=TRUE))
rm(list=grep("_temp",ls(),value=TRUE))

write.csv(all_old_updated_matched_data, file = paste0("output/all_old_updated_matched_data_v",version_no,"_",Sys.Date(),".csv"), row.names = F)




old_updated_paired_dataset_temp = all_old_updated_matched_data %>%
  filter((old_updated_has_pair), !is.na(fold_rise), immunogen_chronology != "matched_immunogen") %>%
  select(c("old_updated_pairing_index","fold_rise","post_neut","immunogen_chronology","paper","paper2","comparison_no","n"))
old_updated_paired_dataset_temp2 = data.frame()

for(i in unique(old_updated_paired_dataset_temp$old_updated_pairing_index)){
  updateds = filter(old_updated_paired_dataset_temp, old_updated_pairing_index==i, immunogen_chronology=="updated_immunogen")
  olds = filter(old_updated_paired_dataset_temp, old_updated_pairing_index==i, immunogen_chronology=="old_immunogen")
  old_updated_paired_dataset_temp2 = rbind(old_updated_paired_dataset_temp2, full_join(olds, updateds, by = c("old_updated_pairing_index","paper","paper2","comparison_no"),suffix = c(".old",".updated"),relationship = "many-to-many"))
}

old_updated_paired_dataset = old_updated_paired_dataset_temp2 %>% 
  select(-c(immunogen_chronology.old,immunogen_chronology.updated))%>%
  mutate(old_updated_pairing_index_2 = row_number()) %>%
  melt(id.vars = c("old_updated_pairing_index","old_updated_pairing_index_2","paper","paper2","comparison_no")) %>%
  mutate(immunogen_chronology = factor(case_match(variable,
                                                  "fold_rise.old"~"old_immunogen",
                                                  "fold_rise.updated"~"updated_immunogen",
                                                  "n.old"~"old_immunogen",
                                                  "n.updated"~"updated_immunogen",
                                                  "post_neut.old"~"old_immunogen",
                                                  "post_neut.updated"~"updated_immunogen"),
                                       levels = immunogen_chronology_levels),
         paper = factor(paper, levels = paper_list_cleaned),
         variable = case_match(as.character(variable),
                               c("fold_rise.old","fold_rise.updated")~"fold_rise",
                               c("n.old","n.updated")~"n",
                               c("post_neut.old","post_neut.updated")~"post_neut")) %>%
  dcast(old_updated_pairing_index+old_updated_pairing_index_2+paper+paper2+comparison_no+immunogen_chronology~...,value.var="value")

rm(old_updated_paired_dataset_temp,old_updated_paired_dataset_temp2)



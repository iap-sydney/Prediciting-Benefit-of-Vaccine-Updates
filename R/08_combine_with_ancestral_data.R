# This function combines a dataframe with matching ancestralk data and removes anythign with only 1 prior dose
combine_with_ancestral_data = function(df, df_ancestral = combined_dataset_to_use, called_from="orig"){
  data_cols = c("crhd_index","crhd_and_comparison_index","pre_neut","post_neut","fold_rise")
  if(called_from=="orig") {
    select_cols = c("paper","paper2","group","prior_doses","prior_status","another_group","immunogen")
  } else if (called_from == "cleaned"){
    select_cols = c("paper","paper2","group","prior_doses_group","prior_doses","prior_status","immunogen_details","immunogen_group", "is_conflicted", "reference","use_in_analysis")
  }
  
  df = df %>%
    left_join(filter(df_ancestral,variant_group=="Ancestral") %>% select(any_of(c(select_cols,data_cols))),
              by = select_cols, suffix = c("",".ancestral")) %>%
    filter(prior_doses>1) %>%
    mutate(comparison_name = paste(old_immunogen,"vs",updated_immunogen,"for",future_variant)) %>%
    distinct()
  if(called_from == "orig"){
    df = df %>% rename(immunogen_details = immunogen) %>%
      mutate(use_in_analysis = TRUE)
  }
  df
}

read_cleaned_neut_data = function(data_file, data_sheet="extracted pre post", data_range="A1:N1000", dataset = NA, data_summary_file = NA){
  neut_data = read_xlsx(path=data_file, sheet = data_sheet, range = data_range) %>%
    as.data.table() %>%
    filter(!is.na(paper))
  paper_list_new = unique(neut_data$paper)
  neut_data = neut_data %>%
    mutate(
      paper = factor(paper, levels = paper_list_new),
      immunogen_group_check = case_when(str_detect(immunogen_details,"XBB")~"XBB",
                                        str_detect(immunogen_details,"BA5") | str_detect(immunogen_details,"BA.5") ~"BA.5",
                                        str_detect(immunogen_details,"BA1") | str_detect(immunogen_details,"BA.1") ~ "BA.1",
                                        str_detect(immunogen_details,"Beta") | str_detect(immunogen_details,"Delta")~ "early VOC",
                                        T~immunogen_group),
      immunogen_group = factor(immunogen_group, levels = immunogen_group_levels_reduced),
      variant_group_check = case_when(str_detect(variant,"XBB")~"XBB",
                                      variant_group %in% c( "BA.2","BA.4.6","CH.1.1") | str_detect(variant,"BQ") | str_detect(variant,"BF") ~"other Omicron",
                                      variant_group %in% c("Beta","Delta")~ "early VOC",
                                      T~"variant_group"), # This can just be used to check we have correctly grouped the variants
      variant_group = factor(variant_group,levels = variant_group_levels_reduced),
      n = N,
      fold_rise = post_neut/pre_neut,
      censored_above = ifelse(((!is.na(UpperLOD)&post_neut>=UpperLOD) | pre_neut<= LOD),TRUE,FALSE),
      censored_completely = if_else((!is.na(UpperLOD)&pre_neut>=UpperLOD & post_neut>=UpperLOD) | (pre_neut<=LOD & post_neut<=LOD),TRUE,FALSE),
      dataset = "cleaned",
      prior_status = case_match(prior_status,
                                c("Naïve","Uninfected","uninfected")~"uninfected",
                                c("Infected","Prior infection", "infected")~"infected",
                                .default = "mixed"),
      propn_prior_infected = ifelse(is.na(propn_prior_infected),case_match(prior_status, "infected"~1,"uninfected"~0, "mixed"~.5),propn_prior_infected),
      crhd_index = row_number()
      
    ) %>% # reorder the columns
    relocate(crhd_index) %>%
    relocate(fold_rise, .after=post_neut) %>%
    relocate(data_from, .after=last_col())
  
  incorrect_variant_groups = with(neut_data, which(as.character(variant_group) != variant_group_check))
  incorrect_immunogen_groups = with(neut_data, which(as.character(immunogen_group) != immunogen_group_check))
  if(length(incorrect_variant_groups)>0){
    #print(paste("Incorrectly grouped variants in rows",incorrect_variant_groups))
  }
  if(length(incorrect_immunogen_groups)>0){
    #print(paste("Incorrectly grouped immunogens in rows",incorrect_immunogen_groups))
  }
  
  neut_data = neut_data %>%
    select(-c(variant_group_check, immunogen_group_check))
  
  if (!is.na(data_summary_file)){
    summary_data_old = read_xlsx(path=data_summary_file) %>%
      as.data.table() %>%
      filter(!is.na(paper)) %>%
      select(c("paper","is_conflicted","reference")) %>%
      mutate(is_conflicted = case_match(is_conflicted,"Yes"~TRUE,.default = FALSE))
    summary_data = summary_data_old
    for(r in c(1:nrow(summary_data))) {
      this_paper = summary_data$paper[r]
      #print(paste("r=",r,this_paper))
      if (!is.na(this_paper) & sum(summary_data$paper[c(r:nrow(summary_data))]==this_paper)>1){
        paper_set = filter(summary_data, paper == this_paper)
        refs = paste(paper_set$reference, collapse = " & ")
        #print(refs)
        summary_data[r,"reference"] = refs
        summary_data = rbind(summary_data[c(1:r),],filter(summary_data[c(r+1:nrow(summary_data)),],paper != this_paper))
      }
    }
    neut_data = left_join(neut_data,summary_data,by = "paper")
  }
}

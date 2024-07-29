
combine_data_by_immunogen = function(immunogen_data_A, immunogen_data_B, paired=TRUE, suffixes = c(".old",".new"), comparison_table, r, combined_immunogen_data = NULL, called_from = "orig"){
  if(is.null(combined_immunogen_data)){
    combined_immunogen_data = data.frame()
  }
  if(called_from =="orig"){
    select_cols = c("paper","paper2","dataset","group","prior_doses","prior_status","variant","variant_group","another_group")
  } else if (called_from == "cleaned") {
    select_cols = c("paper","paper2","dataset","group","prior_doses_group","prior_status","variant","variant_group")
  }
  
  #print(paste(r,":",paired))
  if (paired) {
    combined_immunogen_rows_temp = inner_join(immunogen_data_A,immunogen_data_B,
                                              by = select_cols, suffix = suffixes, relationship="many-to-many") 
  } else {
    combined_immunogen_rows_temp = full_join(immunogen_data_A,immunogen_data_B,
                                             by = select_cols, suffix = suffixes, relationship="many-to-many") 
  }
  combined_immunogen_rows = combined_immunogen_rows_temp %>%
    select(any_of(c(select_cols))) %>%
    mutate(
           pairing_index = as.numeric(paste(row_number(),r,sep=".")),
           old_immunogen = comparison_table[r,"old_immunogen"],
           updated_immunogen = comparison_table[r,"updated_immunogen"],
           future_variant = comparison_table[r,"future_variant"],
           comparison_no = comparison_table[r,"comparison_no"]
    ) 
  
  combined_immunogen_data_temp = rbind(
    inner_join(combined_immunogen_rows, immunogen_data_A, by = c(select_cols),relationship = "many-to-many"),
    inner_join(combined_immunogen_rows, immunogen_data_B, by = select_cols,relationship = "many-to-many")
  ) %>%
    mutate(
      immunogen_chronology = factor(case_when(immunogen_group == comparison_table[r,"old_immunogen"] ~"old_immunogen",
                                                      immunogen_group == comparison_table[r,"updated_immunogen"] ~"updated_immunogen",
                                                      immunogen_group == comparison_table[r,"future_variant"] ~"matched_immunogen",
                                                      TRUE~"not sure"), levels = immunogen_chronology_levels),
      crhd_and_comparison_index = as.numeric(paste(crhd_index, comparison_no, sep="."))
    ) 
  if (called_from == "orig") {
    combined_immunogen_data_temp = combined_immunogen_data_temp %>%
    filter(!is.na(data_group))
  }
  
  # Now we need to remove any doubles.
  # Doubles are those where the only thing that is different if the pairing index
  combined_immunogen_data = rbind(combined_immunogen_data,combined_immunogen_data_temp) %>%
    arrange(pairing_index) %>%
    distinct(pick(!matches("pairing_index")),.keep_all=T)
  
}

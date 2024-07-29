file_name = NA
file_names = NULL
file_starts = c("all_old_updated_matched_data_v","all_old_updated_matched_data_long_v")
# If we have not specified a file to open find the appropriate one.
if(is.na(file_name)){
  for(file_start in file_starts){
  # First try for today's file
  file_name = paste0("output/",file_start,version_no,"_",Sys.Date(),".csv")

  test_this_name = try(read.csv(file_name, fileEncoding="UTF-8-BOM"), silent=T)
  
    # otherwise get the most recent file for this version
    if (inherits(test_this_name, "try-error")) {
      this_version_file_list = file.info(list.files("output/", full.names = T))[grep(paste0(".*",file_start,version_no,".*csv"),list.files("output/", full.names = T)),]
      file_name = rownames(this_version_file_list[which.max(this_version_file_list$mtime),])
      test_this_name = try(read.csv(file_name, fileEncoding="UTF-8-BOM"), silent=T)
    } 
    if (inherits(test_this_name, "try-error")) {
      all_versions_file_list = file.info(list.files("output/", full.names = T))[grep(paste0(".*all_old_updated_matched_data.*csv"),list.files("output/", full.names = T)),]
      file_name = rownames(all_versions_file_list[which.max(all_versions_file_list$mtime),])
      test_this_name = try(read.csv(file_name, fileEncoding="UTF-8-BOM"), silent=T)
    }
    if (inherits(test_this_name, "try-error")) {
      stop("Cannot find a file containing data to open - check there is a file called \"all_old_updated_matched_data.*csv\" in your output directory")
    }
  file_names = c(file_names,file_name)
  }
}
  
regression_data <- read.csv(file_names[1], fileEncoding="UTF-8-BOM") %>%
  mutate(is_conflicted=factor(is_conflicted),
         comparison_no = factor(comparison_no), 
         prior_doses = as.numeric(prior_exposures) - min_prior_doses, #prior doses is prior exposures (infection plus vax doses) - this has been pre-calculated 
         has_more_than_2_exposures = ifelse(prior_doses>0, TRUE, FALSE),
         # subtract off the minimum, which is 2
         fold_rise = log10(fold_rise), 
         post_neut = log10(post_neut), 
         has_old_updated_data = paper %in% old_updated_paper_list
  )
# This just does the old_updated fitting for now
regression_data_old_updated = regression_data[!is.na(regression_data$old_updated_pairing_index),]

regression_data_long <- read.csv(file_names[2], fileEncoding="UTF-8-BOM") %>%
  mutate(is_conflicted=factor(is_conflicted),
         comparison_no = factor(comparison_no), 
         prior_doses = as.numeric(prior_exposures) - min_prior_doses, #prior doses is prior exposures (infection plus vax doses) - this has been pre-calculated 
         has_more_than_2_exposures = ifelse(prior_doses>0, TRUE, FALSE),
         # subtract off the minimum, which is 2
         fold_rise = log10(fold_rise), 
         post_neut = log10(post_neut), 
         has_old_updated_data = paper %in% old_updated_paper_list
  )

# This just does the old_updated fitting for now
regression_data_old_updated_long = regression_data_long[!is.na(regression_data_long$old_updated_pairing_index),]


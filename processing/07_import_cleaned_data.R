updated_immunogen_neut_data = read_cleaned_neut_data(cleaned_data_file,cleaned_data_sheet, cleaned_data_range, "cleaned", cleaned_data_summary)

# These are the ones we previously decided to remove
# !(paper=="Moderna" & immunogen_group =="BA.5"), # this is in chalkias
# !(paper=="Moderna" & prior_status != "mixed"), # This is because separated into inf / uninf exists
# !(paper=="Moderna" & str_detect(group,"subgroup")), # Decided not to include subgroup analysis - need to check why
# !(paper=="Pfizer" & immunogen_group=="BA.5" & variant_group=="BA.5"),
# !(paper=="Carr" & prior_status != "mixed"), # check - assume there was separated into infected / uninfected
# !(paper=="Zou" & prior_status != "mixed")) # check - assume there was separated into infected / uninfected

#!(str_starts(paper,"4Chalkias") | str_starts(paper,"14Davis"))) %>% # Remove 4Chalkias and 14Davis_Gardiner as there are an updated versions from NEJM now included
#!(str_detect(paper,"PressRelease"))) %>%# Remove press releases
#mutate(paper = ifelse(str_starts(paper,"3Chalkias"), "3ChalkiasNatMed", paper)) %>% # rename as it has now been published in Nat med
# filter(prior_doses>1) remove those with 1 or less vaccine doese



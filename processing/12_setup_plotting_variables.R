distinct_ref_labels = distinct(combined_dataset_to_use[,c("paper","reference")])
manuscript_reference_labels = distinct_ref_labels$reference
names(manuscript_reference_labels) = distinct_ref_labels$paper
pattern = "[0-9]+"

for (m in c(1:length(manuscript_reference_labels))) {
  ref_no = as.numeric(regmatches(manuscript_reference_labels[m], regexpr(pattern,manuscript_reference_labels[m])))
  extra_match_start = regexpr("Miller[a-zA-Z .()]+[0-9]+)",manuscript_reference_labels[m])[1]
  extra_match = regmatches(manuscript_reference_labels[m], regexpr("Miller[a-zA-Z .()]+[0-9]+)",manuscript_reference_labels[m]))
  ref_no2 = 0
  new_ref_no = ref_no + ref_start-1
  manuscript_reference_labels[m] = gsub(pattern,new_ref_no,manuscript_reference_labels[m])
  if (length(extra_match)>0){
    new_ref_no2 = as.numeric(regmatches(extra_match, regexpr(pattern,extra_match))) + ref_start-1
    extra_match = gsub(pattern,new_ref_no2,extra_match)
    manuscript_reference_labels[m] = paste0(substr(manuscript_reference_labels[m], 1, extra_match_start-1),extra_match)
  }
}
manuscript_reference_labels = sort(manuscript_reference_labels)

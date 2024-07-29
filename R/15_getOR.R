getOR = function(k="severe", benefit, returnCIs=F,k_lo, k_hi) {
  params = get_natmed2021_params()
  if (k=="severe"){
    k = exp(params$logk_Severe)
    k_lo = exp(params$logk_Severe-1.96*params$logk_severe_sd)
    k_hi = exp(params$logk_Severe+1.96*params$logk_severe_sd)
  } else if (k %in% c("symptomatic", "sympt")){
    k = exp(get_natmed2021_params()$logk)
    k_lo = exp(params$logk-1.96*params$logk_sd)
    k_hi = exp(params$logk+1.96*params$logk_sd)
  }  
  OR = exp(k*log10(benefit))
  if(returnCIs){
    CIs = c(exp(k_lo*log10(benefit)),exp(k_hi*log10(benefit)))
    OR = c(OR,CIs)
  }
  OR
  
}
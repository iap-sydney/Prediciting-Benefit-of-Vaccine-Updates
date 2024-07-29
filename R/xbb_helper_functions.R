single_ttest_pval = function(x){
  if (length(x) <2){
    ans=NA
    
  } else {
      a = t.test(x)
      ans = a$p.value
  }
  ans
}

single_ttest_pvalL = function(x){
  ans = log10(single_ttest_pval(x))
}

weighted.geomean = function(x, w=NULL){

  if (is.null(w) | (length(w)!=length(x))){
    w = rep(1,length(x))
  } 
  
  m = weighted.mean(log10(x[!is.na(x)&!is.na(w)]),w[!is.na(x)&!is.na(w)], na.rm=T)
  ans = 10^m
}


print_significance = function(x){
  cuts = c(1,.1,.05,.01,.001,.0001,0)
  labels = rev(c("ns",".","*","**","***","****"))
  val = cut(x,cuts, labels)
  #print(as.character(val))
  #ans = as.character(val)
}

print_single_ttest_pval = function(x){
  ans = print_significance(single_ttest_pval(x))
}
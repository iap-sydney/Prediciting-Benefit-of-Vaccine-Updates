
mean_plotter <- function(x, y_pos=30, rounding = 0, y_pos_quantile= NULL){
  if (!is.null(y_pos_quantile)){
    y_pos = 10^quantile(x,y_pos_quantile)
  }
  return(data.frame(y = log10(y_pos), label = round(10^mean(x),rounding)))
}

n_plotter <- function(x, y_pos=1){
  return(data.frame(y = log10(y_pos), label = paste0("n=",length(x))))
}

manuscript_plot_by_exposures = function(y = "fold_rise", colouring = "paper", shaping = "comparison_no", plot_data = all_old_updated_matched_data, 
                                        show_numbers = TRUE, ignore_matched_data = TRUE, fix_limits=TRUE, n_label=1) {
  if (shaping == "comparison_no"){
    shaping_str = "factor(comparison_no)"
  } else {
    shaping_str = shaping
  }
  
  if (colouring == "comparison_no"){
    colouring_str = "factor(comparison_no)"
  } else {
    colouring_str = colouring
  }
  if(ignore_matched_data){
    plot_data = filter(plot_data, immunogen_chronology != "matched_immunogen")
  }
  this_plot_data = plot_data %>% 
    filter(round(prior_exposures,0)<=5)
  
  eval(parse(text = paste0(
    "this_plot_data = this_plot_data %>%",
    "filter(!is.na(",y,")) %>%",
    "mutate(yval = ",y,",",
    "colouring = ",colouring_str,",",
    "shaping = ",shaping_str,")"
  )))
  if (y=="fold_rise"){
    y_scale_name = "Fold Rise after Boost"
    stat_comparisons = list(c(2,3),c(1,2))
    text_rounding =1
  } else if (y=="post_neut") {
    y_scale_name = "Neutralisation Titres after Boost"
    stat_comparisons = list(c(1,2),c(2,3))
    text_rounding = 0
  }
  plot_by_exposures_manuscript = ggplot(this_plot_data,
                                        aes(x = factor(round(prior_exposures,0)),y=yval))+
    geom_boxplot(aes(group = factor(round(prior_exposures,0))), outlier.colour = NA, fatten=NULL)+
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
                 width = 0.75, size = 1, linetype = "solid", colour="grey50")+
    stat_summary(fun.data = n_plotter, geom = "text",colour="grey50",vjust = -.5, hjust = .5, fun.args=list(y_pos=n_label)) +
    
    geom_point(aes(fill= colouring, colour=colouring, alpha = log10(log10(as.numeric(n))), shape = shaping), 
               size = point_size, position = position_jitter(width=.1))+
    stat_compare_means(comparisons = stat_comparisons,method="t.test", bracket.size = .3, size = 5)+

    scale_x_discrete(breaks = c(2,3,4), name = "Number of Prior Exposures", expand = expansion(mult = c(0.4,.4)), labels = exposure_labels)+
    scale_alpha_continuous(name = "Subjects in Cohort", breaks = brks, labels = brk_labels, range = alpha_range ) + 
    labs(x="Prior Exposures")+
    manuscript_theme
  if(show_numbers){
    plot_by_exposures_manuscript = plot_by_exposures_manuscript + 
                                    stat_summary(fun.data = mean_plotter, geom = "text",colour="grey40", 
                                               fun.args=list(y_pos=1, rounding = text_rounding, y_pos_quantile = .5), size=5, vjust = -.2, hjust = -.2)
  }
  if(fix_limits){
    plot_by_exposures_manuscript = plot_by_exposures_manuscript +
      scale_y_log10(name = y_scale_name, expand = expansion(mult = c(0.08,.08)), breaks = 10^c(0:3),limits = 10^c(0,2.5))
  } else {
    plot_by_exposures_manuscript = plot_by_exposures_manuscript +
      scale_y_log10(name = y_scale_name, expand = expansion(mult = c(0.08,.1)), breaks = 10^c(0:5))
  }
  # Colouring
  if(this_plot_data$colouring[1]==this_plot_data$paper[1]){  
    plot_by_exposures_manuscript = plot_by_exposures_manuscript +
      scale_colour_manual(name = "Reference", values = paper_colours_cleaned, labels = manuscript_reference_labels)+
      scale_fill_manual(name = "Reference", values = paper_colours_cleaned, labels = manuscript_reference_labels)
  } else if (this_plot_data$colouring[1]==this_plot_data$comparison_no[1]){
    plot_by_exposures_manuscript = plot_by_exposures_manuscript +
      scale_colour_manual(name = "Comparison", values=comparison_colours, labels = names(comparison_list))+
      scale_fill_manual(name = "Comparison", values=comparison_colours, labels = names(comparison_list))
  }
  # Shaping
  if (this_plot_data$shaping[1]==this_plot_data$paper[1]){
    plot_by_exposures_manuscript = plot_by_exposures_manuscript +
      scale_shape_manual(name = "Study", values = paper_shapes_cleaned, labels = manuscript_reference_labels)
  } else if (this_plot_data$shaping[1]==this_plot_data$comparison_no[1]){
    plot_by_exposures_manuscript = plot_by_exposures_manuscript +
      scale_shape_manual(name = "Comparison", values = comparison_shapes, labels = names(comparison_list))
  }
  if (y == "post_neut") {
    plot_by_exposures_manuscript = plot_by_exposures_manuscript + 
      facet_wrap(~future_variant, ncol=1)
  }
  
  #paper_legend = get_legend(plot_by_exposures_manuscript)$grobs[[3]]
  #comparison_legend = get_legend(plot_by_exposures_manuscript)$grobs[[1]]
  #subjects_legend = get_legend(plot_by_exposures_manuscript)$grobs[[2]]
  
  #this_legend = plot_grid(comparison_legend, NULL, paper_legend,NULL, subjects_legend,ncol=1)
  
  #plot_by_exposures = plot_grid(plot_by_exposures_manuscript+theme(legend.position = "none"),
  #                              NULL,
  #                              this_legend, nrow=1, rel_widths = c(1.7,.07,1), rel_heights = c(1,1.2))
  #plot_by_exposures
  plot_by_exposures_manuscript
}



plot_regression_data = function(dataset = all_old_updated_matched_data, yvar="fold_rise", grouping="immunogen_chronology", facet1 = "comparison_no", facet2="prior_exposures", add_shading=T, show_pvals = T, n_label = 10){
  y_variable = case_match(yvar,
                          "norm_neut"~"post_neut_normalised",
                          .default = yvar)
  max_doses = max(dataset$prior_exposures)
  
  # Add some NA data for each category
  dummy_data = as.data.frame(expand.grid(
    immunogen_chronology = unique(dataset$immunogen_chronology),
    comparison_no = unique(dataset$comparison_no),
    prior_exposures = unique(round(dataset$prior_exposures,0)),
    prior_doses = unique(round(dataset$prior_doses,0))
    )) 
  dummy_data = left_join(dummy_data, comparison_table,by="comparison_no") 
  
  compare_rounding_str_start = ""
  compare_rounding_str_end = ""
  if(grouping %in% c("prior_exposures","prior_doses")){
    compare_rounding_str_start = "round("
    compare_rounding_str_end = ",0)"
  }
  facet2_rounding_str_start = ""
  facet2_rounding_str_end = ""
  if(facet2 %in% c("prior_exposures","prior_doses")){
    facet2_rounding_str_start = "round("
    facet2_rounding_str_end = ",0)"
  } else if (facet2 %in% c("comparison_no")){
    
  }
  facet1_rounding_str_start = ""
  facet1_rounding_str_end = ""
  if(facet1 %in% c("prior_exposures","prior_doses")){
    facet1_rounding_str_start = "round("
    facet1_rounding_str_end = ",0)"
  } else if (facet1 %in% c("comparison_no")){
    
  }
  expr_str = paste0("this_data = dataset %>% ",
                    "mutate(compare_group  = factor(",compare_rounding_str_start,grouping,compare_rounding_str_end,"),",
                    "facet1_group  = factor(",facet1_rounding_str_start,facet1,facet1_rounding_str_end,"),",
                    "facet2_group  = factor(",facet2_rounding_str_start,facet2,facet2_rounding_str_end,"),",
                    "y_val = ",y_variable,")")
  eval(expr = parse(text=expr_str))
  this_data = this_data %>%
    mutate(x_group = compare_group,  
           pair_group = paste0(compare_group,old_updated_pairing_index),
           alpha_val = log10(log10(as.numeric(n))),
           colouring = factor(comparison_no),
           shaping = paper) %>%
    filter(!is.na(y_val) | (is.na(y_val)&is.na(paper)),immunogen_chronology != "matched_immunogen", !is.na(x_group))
  
  max_doses = max(this_data$prior_exposures)
  dummy_data = dummy_data %>%
    filter(prior_exposures<=max_doses, prior_doses<=max_doses)
  this_data = this_data %>%
    rbind(dummy_data, fill=TRUE)
    
  
  y_label = case_match(yvar,
                       "fold_rise"~"Fold Rise after Boost",
                       "post_neut"~"Neutralisaion Titre after Boost",
                       "norm_neut"~"Neut cf Ancestral after Boost")
  x_label = case_match(grouping,
                       "prior_exposures"~"Prior Exposures",
                       "prior_doses"~"Prior Doses",
                       "immunogen_chronology"~"Immunogen")
  
  # Add the facetting strings
  comparison_facet_string = with(this_data,paste0("Comparison ",comparison_no,"\n",old_immunogen," vs ",updated_immunogen,"\n for ",future_variant))
  if(facet1 %in% c("comparison_no")){
    this_data$facet1_string = comparison_facet_string
  } else if (facet2 %in% c("comparison_no")) {
    this_data$facet2_string = comparison_facet_string
  }
  exposure_facet_string = with(this_data,paste0("",round(prior_exposures,0)," Exposures"))
  if(facet1 %in% c("prior_exposures")){
    this_data$facet1_string = exposure_facet_string
  } else if (facet2%in% c("prior_exposures")) {
    this_data$facet2_string = exposure_facet_string
  } 

  #this_data = filter(this_data, !is.na(x_group))
  this_plot = ggplot(filter(this_data, !is.na(x_group)), aes(x = x_group,y=y_val))+
    geom_boxplot(aes(group = x_group), outlier.colour = NA, fatten=NULL, width=.3)+
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
                 width = 0.3, size = 1, linetype = "solid", colour="grey50")+
    #stat_summary(fun.data = mean_plotter, geom = "text",colour="grey50", fun.args=list(y_pos=1, rounding = 1, y_pos_quantile = .75), size=5, vjust = -.2, hjust = -.3)+
    geom_point(aes(fill= colouring, colour=colouring, alpha = alpha_val, shape = shaping), size = point_size, position = position_jitter(width=.1))+
    geom_line(aes(group = old_updated_pairing_index, fill= colouring, colour=colouring, alpha = alpha_val, shape = shaping) )+
    stat_summary(fun.data = n_plotter, geom = "text",colour="grey40",vjust = -.5, hjust = .5, fun.args=list(y_pos=n_label)) +
    
    #stat_summary(geom="text",fun=mean,aes(label=round(after_stat(y),1)), vjust = -.1, hjust = -.2,  size=4, colour="grey40")+
    
    scale_y_log10(name = y_label, expand = expansion(mult = c(0.08,.2)))+
    scale_alpha_continuous(name = "# Subjects", breaks = brks, labels = brk_labels, range = c(.1,1) )+
    
    scale_colour_manual(name = "Comparison", values=comparison_colours, labels = names(comparison_list))+
    scale_fill_manual(name = "Comparison", values=comparison_colours, labels = names(comparison_list))+
    scale_shape_manual(name = "Study", breaks = names(paper_shapes_cleaned)[names(paper_shapes_cleaned) %in% names(manuscript_reference_labels)], values = paper_shapes_cleaned,labels = manuscript_reference_labels)+
    
    
    labs(x=x_label)+
    manuscript_theme +
    theme(strip.text.x = element_text(size=14))
  
  if(show_pvals){
    this_plot =this_plot + 
      stat_compare_means(comparisons = list(c(4,5)),method="t.test", bracket.size = .3, size = 5)+
      stat_compare_means(comparisons = list(c(3,4)),method="t.test", bracket.size = .3, size = 5)+
      stat_compare_means(comparisons = list(c(2,3)),method="t.test", bracket.size = .3, size = 5)+
      stat_compare_means(comparisons = list(c(1,2)),method="t.test", bracket.size = .3, size = 5)
  }
  
  if(grouping=="immunogen_chronology"){
    this_plot = this_plot+
      scale_x_discrete(labels = c("old","updated"))
  }
  
  if(!is.na(facet1) | !is.na(facet2)){
    if(is.na(facet2)){
      this_plot = this_plot+
        facet_wrap(~facet1_string,ncol=1)
    } else if(is.na(facet1)){
      this_plot = this_plot + 
        facet_wrap(facet2_string, ncol = length(unique(this_data$facet2_group)))
    }
    else {
      #this_plot = this_plot+
      #  facet_wrap(~paste(facet1_string,facet2_string,sep="\n"), ncol = length(unique(this_data$facet2_group)))
      this_plot = this_plot+
        facet_grid(facet1_string~facet2_string, switch = "y") +
        guides(shape=guide_legend(ncol=1))+
        theme(strip.placement = "outside",
              strip.background = element_rect(fill="white", colour="white"),
              strip.text.x = element_text(face="bold"))
    }
  }
  if(add_shading){
    # Add in the colouring
    ylim = layer_scales(this_plot)$y$range$range
    ycuts = c(seq(floor(ylim[1]),ylim[2],.5),ylim[2])
    xlim = as.numeric(layer_scales(this_plot)$x$range$range)
    x_min = .5
    x_max =length(xlim)+.5
    for (i in seq(2,length(ycuts),2)){
      this_plot = this_plot+
        annotate("rect", xmin = x_min-.1, xmax = x_max+.1, ymin = 10^ycuts[i], ymax = 10^ycuts[i+1],alpha = .2, fill="grey") 
    }
  }
  this_plot
}

manuscript_plot_by_immunogen = function(y = "fold_rise", colouring = "comparison_no", shaping = "paper", plot_data = all_old_updated_matched_data, hjustval = -0.6, vjustval = -.1) {
  if (shaping == "comparison_no"){
    shaping_str = "factor(comparison_no)"
  } else {
    shaping_str = shaping
  }
  
  if (colouring == "comparison_no"){
    colouring_str = "factor(comparison_no)"
  } else {
    colouring_str = colouring
  }
  this_plot_data = plot_data %>% 
    filter(immunogen_chronology != "matched_immunogen", round(prior_exposures,0)<=4)
  eval(parse(text = paste0(
    "this_plot_data = this_plot_data %>%",
    "filter(!is.na(",y,")) %>%",
    "mutate(yval = ",y,",",
    "colouring = ",colouring_str,",",
    "shaping = ",shaping_str,")"
  )))
  if (y=="fold_rise"){
    y_scale_name = "Fold Rise after Boost"
    stat_comparisons = list(c(2,3),c(1,2))
    text_rounding =1
  } else if (y=="post_neut") {
    y_scale_name = "Neutralisation Titres after Boost"
    stat_comparisons = list(c(1,2),c(2,3))
    text_rounding = 0
  }
  plot_by_immungen_manuscript = ggplot(this_plot_data,
                                        aes(x = variant_group,y=yval))+
    geom_boxplot(aes(group = variant_group), outlier.colour = NA, fatten=NULL)+
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
                 width = 0.75, size = 1, linetype = "solid", colour="grey50")+
    geom_point(aes(fill= colouring, colour=colouring, alpha = log10(log10(as.numeric(n))), shape = shaping), 
               size = point_size, position = position_jitter(width=.1))+
    stat_compare_means(comparisons = stat_comparisons,method="t.test", bracket.size = .3)+
    stat_summary(geom="text",fun=mean,aes(label=round(after_stat(y),text_rounding)), vjust = vjustval, hjust = hjustval,  size=4, colour="grey45")+
    
    scale_y_log10(name = y_scale_name, expand = expansion(mult = c(0.08,.08)))+
    scale_x_discrete(name = "Variant Tested", expand = expansion(mult = c(0.4,.4)))+
    scale_alpha_continuous(name = "Subjects in Cohort", breaks = brks, labels = brk_labels, range = c(.1,1) ) + 
    
    manuscript_theme
  
  # Colouring
  if(this_plot_data$colouring[1]==this_plot_data$paper[1]){  
    plot_by_immungen_manuscript = plot_by_immungen_manuscript +
      scale_colour_manual(name = "Reference", values = paper_colours_cleaned, labels = manuscript_reference_labels)+
      scale_fill_manual(name = "Reference", values = paper_colours_cleaned, labels = manuscript_reference_labels)
  } else if (this_plot_data$colouring[1]==this_plot_data$comparison_no[1]){
    plot_by_immungen_manuscript = plot_by_immungen_manuscript +
      scale_colour_manual(name = "Comparison", values=comparison_colours, labels = names(comparison_list))+
      scale_fill_manual(name = "Comparison", values=comparison_colours, labels = names(comparison_list))
  }
  # Shaping
  if (this_plot_data$shaping[1]==this_plot_data$paper[1]){
    plot_by_immungen_manuscript = plot_by_immungen_manuscript +
      scale_shape_manual(name = "Paper", values = paper_shapes_cleaned, labels = manuscript_reference_labels)
  } else if (this_plot_data$shaping[1]==this_plot_data$comparison_no[1]){
    plot_by_immungen_manuscript = plot_by_immungen_manuscript +
      scale_shape_manual(name = "Comparison", values = comparison_shapes, labels = names(comparison_list))
  }
  if (y == "post_neut") {
    plot_by_immungen_manuscript = plot_by_immungen_manuscript + 
      facet_wrap(~future_variant)
  }
  plot_by_immungen_manuscript
}

get_legend_names = function(legends){
  names = NULL
  for (i in c(1:(length(legends)-1))){
    # This extracts te etitle
    names = c(names,legends[[i]]$grobs[[2]]$children[1][[1]]$children[[1]]$label)
  }
  names
}

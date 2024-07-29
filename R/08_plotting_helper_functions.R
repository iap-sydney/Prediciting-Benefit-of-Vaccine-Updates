# Plotting helper functions -----------------------------------------------

make_dose_exposure_plot_0a = function(dataset, yvar="fold_rise", grouping="prior_exposures"){
  y_variable = case_match(yvar,
                          "norm_neut"~"post_neut_normalised",
                          .default = yvar)
  
  if(grouping %in% c("prior_exposures","prior_doses")){
    compare_rounding_str_start = "round("
    compare_rounding_str_end = ",0)"
  }
  expr_str = paste0("this_data = dataset %>% ",
                    "mutate(compare_group  = factor(",compare_rounding_str_start,grouping,compare_rounding_str_end,"),",
                    "y_val = ",y_variable,")")
  
  eval(expr = parse(text=expr_str))
  this_data = this_data %>%
    mutate(x_group = compare_group,  
           pair_group = paste0(compare_group,old_updated_pairing_index),
           alpha_val = log10(log10(as.numeric(n)))) %>%
    filter(!is.na(y_val))
  
  y_label = case_match(yvar,
                       "fold_rise"~"Fold Rise after Boost",
                       "post_neut"~"Neuts after Boost",
                       "norm_neut"~"Neut cf Ancestral after Boost")
  x_label = case_match(grouping,
                       "prior_exposures"~"Prior Exposures",
                       "prior_doses"~"Prior Doses",
                       "immunogen_chronology"~"Immunogen")
  
  this_plot = ggplot(this_data, aes(x = x_group,y=y_val))+
    geom_boxplot(aes(group = x_group), outlier.colour = NA)+
    geom_point(aes(colour=paper, alpha = alpha_val), shape = point_shape, size = point_size, position = position_jitter(width=.1))+
    stat_compare_means(comparisons = list(c(4,5)),method="t.test", bracket.size = .3)+
    stat_compare_means(comparisons = list(c(3,4)),method="t.test", bracket.size = .3)+
    stat_compare_means(comparisons = list(c(2,3)),method="t.test", bracket.size = .3)+
    stat_compare_means(comparisons = list(c(1,2)),method="t.test", bracket.size = .3)+
    scale_y_log10(name = y_label)+
    scale_alpha_continuous(name = "Subjects", breaks = brks, labels = brk_labels, range = c(.1,1) )+
    scale_colour_manual(values = paper_colours_cleaned)+
  
    facet_grid(~immunogen_chronology, scales = "free_y")+
    theme_bw()+
    theme(text = element_text(size=16))+
    labs(x=x_label)
  if(grouping=="immunogen_chronology"){
    this_plot = this_plot+
      scale_x_discrete(labels = c("old","updated"))
  }
  
  this_plot
}

make_dose_exposure_plot_bc = function(dataset, yvar="fold_rise", grouping="prior_exposures", comparison_facet = FALSE, facet2=NA){
  y_variable = case_match(yvar,
                          "norm_neut"~"post_neut_normalised",
                          .default = yvar)
  
  compare_rounding_str_start = ""
  compare_rounding_str_end = ""
  if(grouping %in% c("prior_exposures","prior_doses")){
    compare_rounding_str_start = "round("
    compare_rounding_str_end = ",0)"
  }
  facet_rounding_str_start = ""
  facet_rounding_str_end = ""
  if(facet2 %in% c("prior_exposures","prior_doses")){
    facet_rounding_str_start = "round("
    facet_rounding_str_end = ",0)"
  }
  expr_str = paste0("this_data = dataset %>% ",
                    "mutate(compare_group  = factor(",compare_rounding_str_start,grouping,compare_rounding_str_end,"),",
                    "facet_group  = factor(",facet_rounding_str_start,facet2,facet_rounding_str_end,"),",
                    "y_val = ",y_variable,")")
  eval(expr = parse(text=expr_str))
  this_data = this_data %>%
    mutate(x_group = compare_group,  
           pair_group = paste0(compare_group,old_updated_pairing_index),
           alpha_val = log10(log10(as.numeric(n)))) %>%
    filter(!is.na(y_val),immunogen_chronology != "matched_immunogen")
  
  y_label = case_match(yvar,
                       "fold_rise"~"Fold Rise after Boost",
                       "post_neut"~"Neuts after Boost",
                       "norm_neut"~"Neut cf Ancestral after Boost")
  x_label = case_match(grouping,
                       "prior_exposures"~"Prior Exposures",
                       "prior_doses"~"Prior Doses",
                       "immunogen_chronology"~"Immunogen")
  if (!is.na(facet2)){
    # Add some NA data for each category
    dummy_data = as.data.frame(expand.grid(
        x_group = unique(this_data$x_group),
        comparison_no = unique(this_data$comparison_no),
        facet_group = unique((this_data$facet_group))
      )) %>%
      left_join(comparison_table,by="comparison_no")
    
    this_data = this_data %>%
      rbind(dummy_data, fill=TRUE)
  }
  this_plot = ggplot(this_data, aes(x = x_group,y=y_val))+
    geom_boxplot(aes(group = x_group), outlier.colour = NA)+
    geom_point(aes(fill= factor(comparison_no), colour=factor(comparison_no), alpha = alpha_val, shape = paper), size = point_size, position = position_jitter(width=.1))+
    stat_compare_means(comparisons = list(c(4,5)),method="t.test", bracket.size = .3)+
    stat_compare_means(comparisons = list(c(3,4)),method="t.test", bracket.size = .3)+
    stat_compare_means(comparisons = list(c(2,3)),method="t.test", bracket.size = .3)+
    stat_compare_means(comparisons = list(c(1,2)),method="t.test", bracket.size = .3)+
    stat_summary(geom="text",fun=mean,aes(label=round(after_stat(y),1)), vjust = -.1, hjust = -.2,  size=4, colour="brown")+
    
    scale_y_log10(name = y_label, expand = expansion(mult = c(0,.1)))+
    scale_alpha_continuous(name = "# Subjects", breaks = brks, labels = brk_labels, range = c(.1,1) )+
    scale_colour_manual(name = "Comparison", values=comparison_colours, labels = names(comparison_list))+
    scale_fill_manual(name = "Comparison", values=comparison_colours, labels = names(comparison_list))+
    scale_shape_manual(name = "Paper", values = paper_shapes_cleaned)+
    
    labs(x=x_label)+
    xbb_manuscript_theme
  if(grouping=="immunogen_chronology"){
    this_plot = this_plot+
      scale_x_discrete(labels = c("old","updated"))
  }
  
  if(comparison_facet | !is.na(facet2)){
    sep1=" "
    sep2="\n"
    sep3 = " "
    if(is.na(facet2)){
      this_plot = this_plot+
        facet_wrap(~paste0("(",comparison_no,")",sep1,old_immunogen," vs ",updated_immunogen,sep2,"for ",future_variant),ncol=1)
    } else if(!comparison_facet){
      this_plot = this_plot + 
        facet_wrap(~paste0("(",facet_group," exposures)"), ncol = length(unique(this_data$facet_group)))
    }
    else if (comparison_facet & !is.na(facet2)){
      this_plot = this_plot+
        facet_wrap(~paste0("(",comparison_no,")",sep1,old_immunogen," vs ",updated_immunogen,sep2,"for ",future_variant,sep3, "(",facet_group," exposures)"), ncol = length(unique(this_data$facet_group)))
    }
    ylim = layer_scales(this_plot)$y$range$range
    ycuts = c(seq(floor(ylim[1]),ylim[2],.5),ylim[2])
    xlim = as.numeric(layer_scales(this_plot)$x$range$range)
    x_min = .5
    x_max =length(xlim)+.5
    for (i in seq(1,length(ycuts),2)){
      this_plot = this_plot+
      annotate("rect", xmin = x_min, xmax = x_max, ymin = 10^ycuts[i], ymax = 10^ycuts[i+1],alpha = .2, fill="grey")
    }
  }
  ylim = layer_scales(this_plot)$y$range$range
  this_plot
}

make_dose_exposure_plot_ef = function(dataset, yvar="fold_rise", grouping="prior_exposures"){
  y_variable = case_match(yvar,
                          "norm_neut"~"post_neut_normalised",
                          .default = yvar)
  
  if(grouping %in% c("prior_exposures","prior_doses")){
    compare_rounding_str_start = "round("
    compare_rounding_str_end = ",0)"
  }
  expr_str = paste0("this_data = dataset %>% ",
                    "mutate(compare_group  = factor(",compare_rounding_str_start,grouping,compare_rounding_str_end,"),",
                    "y_val = ",y_variable,")")
  eval(expr = parse(text=expr_str))
  this_data = this_data %>%
    mutate(x_group = factor(paste(compare_group,immunogen_chronology,sep="_")),  
           pair_group = paste0(compare_group,old_updated_pairing_index),
           alpha_val = log10(log10(as.numeric(n)))) %>%
    filter(!is.na(y_val))#, immunogen_chronology!="matched_immunogen")
  title_str = "Using All Data"
  if (grouping == "prior_exposures"){
    this_data = filter(this_data,prior_status!="mixed")
    title_str = "No Mixed Infection Data"
  }
  
  
  this_data_text = this_data %>%
    group_by(x_group, compare_group) %>%
    summarise(y_val_meanL = mean(log10(y_val)), .groups = "keep") %>%
    group_by(compare_group)%>%
    summarise(y_val_ratio = 10^diff(y_val_meanL), .groups = "keep")
  
  y_label = case_match(yvar,
                       "fold_rise"~"Fold Rise after Boost",
                       "post_neut"~"Neuts after Boost",
                       "norm_neut"~"Neut cf Ancestral after Boost")
  x_label = case_match(grouping,
                       "prior_exposures"~"Prior Exposures",
                       "prior_doses"~"Prior Doses")
  this_plot = ggplot(this_data, aes(x = x_group,y=y_val))+
    geom_boxplot(aes(group = x_group),outlier.colour = NA)+
    geom_point(aes(colour=paper, alpha = alpha_val), shape = point_shape, size = point_size, position = position_jitter(width=.1))+
    geom_line(aes(colour=paper, alpha = alpha_val, group = pair_group))+
    #stat_compare_means(comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10)),method="t.test")+
    stat_summary(aes(label = round(after_stat(y),1)), geom="text",fun=mean)+
    geom_text(data = this_data_text, aes(x=paste0(compare_group,"_old_immunogen"), y=max(this_data$y_val),label=paste("Ratio =",round(y_val_ratio,2))),nudge_x=.5)+
    
    
    scale_y_log10(name = y_label)+
    scale_alpha_continuous(name = "Subjects", breaks = brks, labels = brk_labels, range = c(.1,1) )+
    scale_colour_manual(values = paper_colours_cleaned)+
    
    theme_bw()+
    theme(text = element_text(size=16),
          axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
    labs(x=x_label, title = title_str)  
  this_plot
}

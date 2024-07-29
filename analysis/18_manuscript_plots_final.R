manuscript_plot_XBB_variant = ggplot(filter(combined_dataset_to_use, !is.na(immunogen_group), variant_group == "XBB", !is.na(fold_rise), immunogen_group != "XBB"),
                                     aes(x=immunogen_group, y = fold_rise))+
  geom_boxplot(aes(group=immunogen_group), outlier.color = NA)+
  #stat_summary(geom="text",fun=mean,aes(label=round(after_stat(y),2)), vjust = -.1, hjust = -.2,  size=4)+
  stat_summary(fun.data = mean_plotter, geom = "text",colour="grey40", fun.args=list(y_pos=1, rounding = 1, y_pos_quantile = .5), size=4,vjust = -.5, hjust = -.25)+
  stat_summary(fun.data = n_plotter, geom = "text",colour="grey40", size=6,vjust = -.5, hjust = .5) +
  
  stat_compare_means(comparisons = list(c(1,2),c(2,3)),method="t.test", size = 5)+
  geom_point(aes(colour=immunogen_group, shape = paper), size = 2, position = position_jitter(width=.3))+
  scale_colour_manual(values = c("grey50", "grey50", "grey50"), guide="none")+
  scale_y_log10(breaks = c(1,3,10,30))+
  scale_shape_manual(name = "Paper", values = paper_shapes_cleaned,labels = manuscript_reference_labels)+
  manuscript_theme+
  labs(y = "Fold Rise to XBB",x = "Immunogen in Booster")+ 
  coord_cartesian(ylim = c(1,50))
manuscript_plot_XBB_variant
ggsave("output/plots/18_XBB_variant.pdf",manuscript_plot_XBB_variant, width=mw*1.2, height = mh*1.5) 


fold_rise_regression_plot = plot_regression_data(yvar = "fold_rise", facet2="comparison_no",facet1 = "prior_exposures", show_pvals=F, n_label = .3)

legends = get_legend(fold_rise_regression_plot + guides(shape=guide_legend(ncol=3))+theme(legend.text = element_text(size=11)))$grobs
legend_names = get_legend_names(legends)

paper_legend_no = which(legend_names == "Study")
comparison_legend_no = which(legend_names == "Comparison")
subjects_legend_no = which(endsWith(legend_names,"Subjects"))
comparison_legend = legends[[comparison_legend_no]]
paper_legend_3col = legends[[paper_legend_no]]
paper_legend_1col = get_legend(fold_rise_regression_plot + guides(shape=guide_legend(ncol=1))+theme(legend.text = element_text(size=11)))$grobs[[paper_legend_no]]
subjects_legend =  legends[[subjects_legend_no]]

horizontal_legend = plot_grid(NULL,comparison_legend,NULL,subjects_legend,paper_legend_3col,NULL,nrow=1, rel_widths = c(.1,1.4,.1,.8,3.5,.1), align = "hv")
vertical_legend = plot_grid(comparison_legend, NULL, paper_legend_1col,NULL, subjects_legend,ncol=1)

dummy_data = generate_dummy_data()
data_points = generate_plot_data_points()

model_plot_no_interaction_fold_rise_only = ggplot(filter(dummy_data, interaction=="No", variable=="fold_rise"), aes(x=prior_exposures+jitter,y=10^value, colour=factor(comparison_no), lty = immunogen_chronology))+
  geom_point(data = filter(data_points, variable=="fold_rise"), 
             aes(y=value, shape = paper), size=2.5) +
  geom_line(size=1.3) +
  facet_grid(variable_label~imm_chron_label, scales = "free_y", switch="y")+
  scale_y_log10(name = "Fold Rise after Boost") +
  scale_shape_manual(name = "Paper", values = paper_shapes_cleaned, labels = manuscript_reference_labels)+
  scale_linetype_manual(name = "Model Fit",values = c("old_immunogen"="solid","updated_immunogen" = "dashed"), labels = c("Old Immunogen Fit","Updated Immunogen Fit"))+
  scale_colour_manual(name = "Comparison" , values = c(comparison_colours,"any"="black"), labels = c(names(comparison_list),"Model Fit")) +
  scale_x_continuous(name = "Prior Exposures", breaks = c(2:4), labels = head(exposure_labels,-1))+
  labs(y= "")+
  manuscript_theme +
  theme(strip.text.y = element_blank())

ylim = c(0,4)
xlim = as.numeric(layer_scales(model_plot_no_interaction_fold_rise_only)$x$range$range)
ycuts = c(seq(floor(ylim[1]),ylim[2],.5),ylim[2])
for (i in seq(1,length(ycuts),2)){
  model_plot_no_interaction_fold_rise_only = model_plot_no_interaction_fold_rise_only+
    annotate("rect", xmin = xlim[1], xmax = xlim[2], ymin = 10^ycuts[i], ymax = 10^ycuts[i+1],alpha = .2, fill="grey")
}

model_plot_no_interaction_fold_rise_only = model_plot_no_interaction_fold_rise_only +
  coord_panel_ranges(panel_ranges = list(
    list(x=xlim,y=c(-0.15,2.3)), # Panel 1
    list(x=xlim,y=c(-0.15,2.3)) # Panel 2
  ))

model_plot_no_interaction_fold_rise_only_tagged = tag_facet(model_plot_no_interaction_fold_rise_only, size=5)

regression_fold_rise_and_model_plot = plot_grid(fold_rise_regression_plot+theme(legend.position = "none"),
                                                NULL,
                                                model_plot_no_interaction_fold_rise_only+theme(legend.position = "none"),
                                                horizontal_legend, ncol=1, rel_heights = c(5.5,.2,3,1), align = "h", labels = c("A","","B"), label_x=0,label_y = c(1,1,1.2), label_size = 30)
ggsave("output/plots/18_regression_fold_rise_data_and_model.pdf", regression_fold_rise_and_model_plot,width=mw*2.1, height=mh*4.5)


fold_rise_plot_by_exposures_manuscript = manuscript_plot_by_exposures (y="fold_rise", colouring = "comparison_no", shaping = "paper", n_label=0)

fold_rise_plot_by_exposures_manuscript_legend = plot_grid(fold_rise_plot_by_exposures_manuscript+theme(legend.position = "none"),
                                                          NULL,
                                                          vertical_legend, nrow=1, rel_widths = c(1.7,.07,1.1), rel_heights = c(1,1.2))
ggsave("output/plots/18_fold_rise_by_exposures.pdf",fold_rise_plot_by_exposures_manuscript_legend, width=mw*1.45, height=mh*1.75)


post_neut_regression_plot = plot_regression_data(yvar = "post_neut", facet2="comparison_no",facet1 = "prior_exposures", show_pvals=F, n_label = 10)

absolute_legends = get_legend(post_neut_regression_plot + guides(shape=guide_legend(ncol=3))+theme(legend.text = element_text(size=11)))$grobs
absolute_legend_names = get_legend_names(absolute_legends)

paper_absolute_legend_no = which(absolute_legend_names=="Study")
paper_legend_3col_absolute = absolute_legends[[paper_absolute_legend_no]]

horizontal_aboslute_legend = plot_grid(NULL,comparison_legend,NULL,subjects_legend,paper_legend_3col_absolute,NULL,nrow=1, rel_widths = c(.1,1.4,.1,.8,3.5,.1), align = "hv")
regression_absolute_data = plot_grid(post_neut_regression_plot+theme(legend.position = "none"),
                                     NULL,
                                     horizontal_aboslute_legend, ncol=1, rel_heights = c(5.5,.1,1), 
                                     align = "h")
ggsave("output/plots/18_regression_absolute_data.pdf",regression_absolute_data,width=mw*2.2, height=mh*4)


eff_old_updated_plot = ggplot(mutate(effs_old_updated_df, letterLabel = head(intToUtf8(2*as.numeric(data_type)+(as.numeric(outcome)-1)+63))),
                              aes(x=old_boost_eff*100, y = 100*value, colour = factor(immunogen_chronology)))+
  geom_line(size = 1.5)+
  labs(x = "Old Booster Symptomatic Protection\n(vs Naive)", y = "")+# y = "Updated Booster")+
  scale_colour_manual(name="Immunogen",values = c("orange","chocolate3"), labels = c("Old Immunogen", "Updated Immunogen"))+
  annotate("rect",xmin=30, xmax = 80, ymin = 0, ymax = 100, colour = "grey90", fill = "grey90")+
  annotate("segment",x=60,y = 0, xend=60,yend = 100, colour = "grey40", lty="dashed")+
  #annotate("text", label = round(filter(effs_old_updated_df,data_type==eff_plot_factors[3], old_boost_eff==.3)$value,2)*100, x = 10, y = 10)+
  
  geom_text(aes(label=letter_label), x = 2, y = 98, size = 8, colour="black")+
  geom_line(size = 1.5)+
  facet_grid(data_type~outcome, switch="y", scales = "free_y")+
  
  xbb_manuscript_theme+
  theme(legend.key.width = unit(40,"pt"))
ggsave("output/plots/18_efficacies_old_updated.pdf",eff_old_updated_plot,width = mw*2.1, height = mh*3.6)

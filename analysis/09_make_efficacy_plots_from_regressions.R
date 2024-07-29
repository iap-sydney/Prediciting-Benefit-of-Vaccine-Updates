eff_plot = ggplot(filter(eff_df,data_type!="Difference from Baseline"),
                  aes(x=baseline_eff*100, y = 100*value, colour = factor(doses), linetype = immunogen_chronology))+
  geom_line(size = 1.5)+
  labs(x = "Baseline Population Protection", y = "Effectiveness")+
  scale_linetype_manual(name = "Scenario",
                        values = c("dashed","solid"),
                        labels = c("Older Immunogen","Updated Immunogen"))+
  scale_colour_manual(name="Prior Doses",values = c("grey","orange","darkorange1","chocolate3"), labels = c("Baseline Protection","Primary Course","3 Exposures","4 Exposures"))+

  facet_grid(data_type~outcome, switch="y", scales = "free_y")+
  
  xbb_manuscript_theme+
  theme(legend.key.width = unit(40,"pt"))
ggsave("output/plots/91_efficacies.pdf",eff_plot,width = w*1.3, height = h*1.3)


eff_old_updated_plot = ggplot(effs_old_updated_df,
                    aes(x=old_boost_eff*100, y = 100*value, colour = factor(immunogen_chronology)))+
    geom_line(size = 1.5)+
    labs(x = "Old Booster Symptomatic Protection\n(vs Naive)", y = "")+# y = "Updated Booster")+
    scale_colour_manual(name="Immunogen",values = c("orange","chocolate3"), labels = c("Old Immunogen", "Updated Immunogen"))+
    annotate("rect",xmin=30, xmax = 80, ymin = 0, ymax = 100, colour = "grey90", fill = "grey90")+
    annotate("segment",x=60,y = 0, xend=60,yend = 100, colour = "grey40", lty="dashed")+
    
    geom_line(size = 1.5)+
    facet_grid(data_type~outcome, switch="y", scales = "free_y")+
    
    xbb_manuscript_theme+
    theme(legend.key.width = unit(40,"pt"))
ggsave("output/plots/92_efficacies_old_updated.pdf",eff_old_updated_plot,width = w*1.3, height = h*1.8)
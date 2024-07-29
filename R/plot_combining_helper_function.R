combine_plots = function(a,b,rel_widths=1,nrow=1,ncol=NULL){
  combined_plot = cowplot::plot_grid(a+theme(legend.position = "none"),
                     b+theme(legend.position = "none"),
                     get_legend(b),
                     nrow = nrow, 
                     rel_widths = rel_widths,
                     labels=c("A","B",""))
}
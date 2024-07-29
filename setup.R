library(iaputils)
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(stringr)
library(ggpmisc)
library(glmmTMB)

iaputils::generate_default_project_structure()
# This version incldues the new beta_BA.1 data and delta_BA.1 data as BA.1 immunogens...
version_no = "2.1"

cleaned_data_file = paste0("raw-data/UpdatedImmunogen_NeutData_v",version_no,".xlsx")
cleaned_data_summary = paste0("raw-data/UpdatedImmunogen_MetaData_v",version_no,".xlsx")
cleaned_data_sheet = "2023"
cleaned_data_range = "C1:Y10000"

min_prior_doses = 2

manufacturer_colours = c("Pfizer"= "blue",
                  "Moderna"= "red",
                  "Novavax"= "green4")

w = 8
h=6
mh = 3
mw = 5
jitter_width = .2
small_jitter_width = .08
neut_norm_factor = 1
wsize=15
wround=1
include_abs_neuts_in_plots = FALSE
alpha_range = c(.25,1)
max_dec_pl = 12


variant_levels = c("Ancestral","BA.5","XBB.1.5","XBB.1.16","XBB.2.3.2","BQ.1.1","CH.1.1")
variant_group_levels = c("Ancestral","Beta","Gamma","Delta","BA.1","BA.2","BA.4.6","BA.5","XBB","other Omicron","bivalent","BQ","BF")
variant_group_combined_levels = c("Ancestral","early VOCs","BA.1","BA.5","XBB","other Omicron")
immunogen_levels = c("Ancestral","BA.1","BA.5","bivalent","XBB","XBB.1.5","XBB.1.16","BA.5_XBB.1.5","BA.5_XBB.1.16", "Beta", "Delta", "early_VOC")
immunogen_group_levels = c("Ancestral","early VOC","BA.1","BA.5","XBB","bivalent")

variant_group_levels_reduced = c("Ancestral","BA.1","BA.5","XBB","other BA Omicron","other non BA Omicron","early VOC")
immunogen_group_levels_reduced = c("Ancestral","early VOC","BA.1","BA.5","XBB")

immunogen_chronology_levels = c("old_immunogen","updated_immunogen","matched_immunogen","not sure")

species_levels = c("Human","Mouse")
dataset_levels = c("FDA"="fda","VIEW-hub"="viewhub","VSV"="vsv", "Extra"="extra")

variant_modified_variant_groups = c("Ancestral"="Ancestral","Beta"="Beta","Delta"="Delta","BA.1"="Omicron BA.1",
                                          "BA.5"="Omicron BA.4/5","Gamma"="Gamma","BA.2"="Omicron BA.2",     
                                          "BA.4.6"= "Omicron BA.4.6",
                                          "other Omicron"="BA.2.75","other Omicron"="Omicron BA.2.75","other Omicron"="Omicron BA.2.75.2","other Omicron"="BA.2.75.2",
                                          "BF"="OmicronBF.7","BQ"="OmicronBQ.1.1",
                                          "XBB"="OmicronXBB.1","BQ"="BQ.1.1"
                                          )

variant_shapes = c(16,98,103,100,7,8,3,15,17,13,96,10,12)
names(variant_shapes) = variant_group_levels
variant_shapes = c(variant_shapes,"XBB.1.5"=24,"XBB.1.16"=25, "XBB.2.3.2"=23)
paper_list = c("Branche (1) CID", "Branche (2) RS", "Carr", "Lasrado", "Zou","Chalkias","Moderna","Pfizer")
paper_colours = c("pink","skyblue","palegreen","goldenrod","plum","indianred","red","blue")
names(paper_colours) = paper_list

ref_start = 1
paper_list_cleaned = c("ModernaFDA2023", "ChalkiasNEJM", "ChalkiasMedrxiv","ChalkiasNatMed2022","ChalkiasNatMed2023","ChalkiasMedrxivXBB.1.5", 
                       "BrancheCID", "BrancheNatMed", 
                       "PfizerFDA2022","PfizerFDA2023","ZouNEJM", 
                       "Collier+MillerNEJM",
                       "TanLancetID", "KawasujiMicroSpect","DavisGardinerNEJM", "CarrLancetID","RosslerNatComm","ChoiNatMed", 
                       "AddetiaNature", "LaunayNEJM",  "PajonNEJM", "WangNEJM","KurhadeNatMed",
                       "HoffmannLancetID","GravensteinMedrxiv","HeJClinMed","UrschelMedrxiv","JiangBiorxiv","WangCell","WangJID")
paper_colours_cleaned = c("red3", "tomato3","coral","indianred","darkred","red",
                          "purple","chocolate1",
                          "blue","royalblue","steelblue",
                          "seagreen",
                          "lightblue", "slategray3","plum","cyan","deeppink","lightgreen",
                          "yellow2","greenyellow","skyblue","bisque3","darkgoldenrod1",
                          "magenta","green","pink","aquamarine","grey","yellow","gold")
potential_shapes = c(c(1:25),c(35:38),c(33,34,47,63,64))
paper_shapes_cleaned = potential_shapes[c(1:length(paper_list_cleaned))]

comparison_list = c(1:4)
long_comparison_list = c(1:6)
names(comparison_list) = c("(1) Anc VS Early VOC for BA.1","(2) Anc VS BA.1 for BA.5",
                          "(3) Early VOC VS BA.1 for BA.5","(4) BA.1 VS BA.5 for XBB")
names(long_comparison_list) = c(names(comparison_list),"(5) Ancestral VS BA.1 for XBB","(6) Ancestral VS BA.5 for XBB")

comparison_colours = c("red","blue","darkgreen","yellow3")
comparison_colours = c("#4B9BD5","#B8AED6","#8876B6","#EE2C24")
comparison_colours = c("#EE2C24","#4B9BD5","#14B794","#F6A71C")
long_comparison_colours = c("#EE2C24","#4B9BD5","#14B794","#F6A71C","#B8AED6","#8876B6")

names(comparison_colours) = comparison_list
names(long_comparison_colours) = long_comparison_list

comparison_shapes = c(21:24)
long_comparison_shapes = c(21:26)
#comparison_shapes = c(19,15,18,17)
names(comparison_shapes) = comparison_list
names(long_comparison_shapes) = long_comparison_list


comparison_table = data.frame(old_immunogen = c("Ancestral","Ancestral","early VOC","BA.1"),
                              updated_immunogen = c("early VOC","BA.1","BA.1","BA.5"),
                              future_variant = c("BA.1","BA.5","BA.5","XBB")) %>%
  mutate(comparison_no = row_number())

long_comparison_table = data.frame(old_immunogen = c("Ancestral","Ancestral","early VOC","BA.1","Ancestral","Ancestral"),
                              updated_immunogen = c("early VOC","BA.1","BA.1","BA.5","BA.1","BA.5"),
                              future_variant = c("BA.1","BA.5","BA.5","XBB","XBB","XBB")) %>%
  mutate(comparison_no = row_number())
  
names(paper_colours_cleaned) = paper_list_cleaned
names(paper_shapes_cleaned) = paper_list_cleaned
paper_order = c(19,7,8,16,18,2,4,5,3,6,12,15,25,26,24,14,28,23,20,1,21,9,10,17,13,27,30,22,29,11)
paper_colours_cleaned = paper_colours_cleaned[paper_order]
paper_shapes_cleaned = paper_shapes_cleaned[paper_order]

paper_jitter_denom = length(paper_list_cleaned)/(2*jitter_width)

xbb_manuscript_theme = theme_bw()+
                             theme(text = element_text(size=24), 
                                   legend.title = element_text(size=18),
                             legend.text = element_text(size=12),
                             legend.key.size = unit(15, units = "pt"),
                             legend.spacing = unit(4,units="pt"))
xbb_manuscript_theme_angled = xbb_manuscript_theme+
  theme(axis.text.x = element_text(angle=45, vjust=1,hjust=1))


exposure_labels = c(expression(2~dose~(1^o)),expression(1^o+1),expression(1^o+2),expression(1^o+3))

manuscript_theme = theme_bw()+
  theme(legend.text = element_text(size=10), 
        legend.key.height = unit(12,"pt"), 
        legend.key.width = unit(4,"pt"),
        text = element_text(size=14),
        #legend.box.margin = margin(t=1,b=1,unit="pt"),
        legend.spacing.y = unit(2,"pt"),
        axis.text.y = element_text(size=20, colour = "black"),
        axis.text.x = element_text(size=20, colour = "black", angle = 25, vjust = 1, hjust = 1),
        axis.title = element_text(size=20, colour = "black"),
        axis.title.x = element_text(vjust = -.5, face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.line = element_line(linewidth = 1),
        panel.border = element_rect(colour=NA),
        panel.grid = element_line(colour=NA),
        strip.background = element_rect(fill="white", colour="white"),
        strip.placement = "outside",
        strip.text = element_text(size=20, colour = "black")
        )


brk_ns = c(10,30,100,300)
brks = log10(log10(brk_ns))
brk_labels = 10^10^brks
point_shape = 16
point_size = 3
# Helper functions --------------------------------------------------------

get_log10_p_symbol = function(x){
  if (length(x)<2){
    ans = NA
  } else {
    ans = t.test(log10(x), mu = 0)$p.value
    cutpoints = c(0.001, 0.01, 0.05, 0.1, 1)
    symbols = c("***", "**", "*", ".", "-")
    ans = symbols[ans<=cutpoints][1]
  }
  ans
}

# List of papers that include neut data for a future variant for both old and updated immunogens
old_updated_paper_list = c("AddetiaNature",
                           "BrancheNatMed",
                           "ChalkiasNatMed2022",
                           "HeJClinMed",
                           "LaunayNEJM",
                           "PajonNEJM",
                           "ChalkiasNEJM",
                           "KawasujiMicroSpect",
                           "TanLancetID",
                           "BrancheCID",
                           "RosslerNatComm")

mixed_eff_params_colnames = c("model_name","dependant_var","with_exposure_cts","with_exposure_2_cat","includes_interaction","with_conflicted","includes_valency","includes_old","includes_updated","includes_matched",
                              "AIC","intercept","intercept_pval","updated_benefit_param","updated_benefit_pval","exp_cts_param","exp_cts_pval","int_cts_coeff","int_cts_pval",
                              "exp_cat_param","exp_cat_pval","int_cat_coeff","int_cat_pval","conflicted_coeff","conflicted_pval","is_bivalent_coeff","is_bivalent_pval","intercept_stderr","updated_benefit_stderr","exp_cts_stderr","int_cts_stderr","int_cat_stderr","is_bivalent_stderr")

eff_plot_factors = c("Protection\n(vs Naive)","Increase in\nProtection","Relative Efficacy\n(updated vs old)")


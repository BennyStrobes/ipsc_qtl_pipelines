args = commandArgs(trailingOnly=TRUE)
library(edgeR)
library(ggplot2)
library(ggthemes)
library(reshape)
library(cowplot)





eqtl_comparison_to_reference_bar_plot <- function(file_stem, output_file) {
    # First extract data. And get into nice data format
    pvalues <- c()
    version <- c()
    time_step <- c()
    for (temp_time_step in 0:15) {
        ipsc_file <- paste0(file_stem, temp_time_step,"_real_v_matched_controls.txt")
        data <- read.table(ipsc_file,header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("real", length(data$real_pvalue))))
        pvalues <- c(pvalues,data$matched_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$matched_pvalue)))
        version <- c(version, as.character(rep("matched", length(data$matched_pvalue))))
    }
    df <- data.frame(pvalues = as.numeric(pvalues), version = factor(version,c("real","matched")), time_step = factor(time_step))

    box_plot <- ggplot(df, aes(x=time_step, y=pvalues, fill=version)) + geom_boxplot(width=.54)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(fill= "Test version",x = "time step", y = "pvalue")
    ggsave(box_plot, file=output_file,width = 20,height=10.5,units="cm")

}










visualization_directory = args[1]
distance = args[2]
maf_cutoff = args[3]


###############################################
# Plot distribution of pvalues found in our eqtl data
# But only those variant-gene pairs that are found in:
## a. Nick Banovich's ipsc data
ipsc_file_stem <- paste0(visualization_directory, "ipsc_banovich_comparison_distance_",distance,"_maf_cutoff_",maf_cutoff,"_time_step_")
ipsc_plot_file <- paste0(visualization_directory, "ipsc_banovich_comparison_distance_",distance,"_maf_cutoff_",maf_cutoff, "_bar_plot.png")
eqtl_comparison_to_reference_bar_plot(ipsc_file_stem, ipsc_plot_file)

## b. GTEx Heart left ventricle data
gtex_file_stem <- paste0(visualization_directory, "gtex_v7_heart_left_ventricle_comparison_distance_",distance,"_maf_cutoff_",maf_cutoff,"_time_step_")
gtex_plot_file <- paste0(visualization_directory, "gtex_v7_heart_left_ventricle_comparison_distance_",distance,"_maf_cutoff_",maf_cutoff, "_bar_plot.png")
eqtl_comparison_to_reference_bar_plot(gtex_file_stem, gtex_plot_file)


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



make_pvalue_histogram_one_time_step <- function(input_file, time_step) {
    all_eqtl_nominal_pvalues <- read.table(input_file, header=TRUE)
    pvalues <- all_eqtl_nominal_pvalues$pvalue
    # pvalues <- runif(30000, 0.0, 1.0)
    # Colors for histogram
    barfill <- "#4271AE"
    barlines <- "#1F3552"
    # Title of plot
    title <- paste0("t = ", time_step)
    # put data into data frame for plotting
    df <- data.frame(pvalues = pvalues)
    # make plots
    p <- ggplot(df, aes(x = pvalues)) +
        geom_histogram(aes(y = ..count..), colour = barlines, fill = barfill,binwidth = 0.01, boundary = 0) +
        scale_x_continuous(name = "pvalue", breaks = round(seq(0,1,by=.5),1)) + theme(axis.text=element_text(size=13)) +
        scale_y_continuous(name = "counts") + 
        ggtitle(title)
    return(p)
}

pvalue_histogram_across_time_steps <- function(input_stem, output_file) {
    time_step <- 0
    p0 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 1
    p1 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 2
    p2 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 3
    p3 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 4
    p4 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 5
    p5 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 6
    p6 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 7
    p7 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 8
    p8 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 9
    p9 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 10
    p10 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 11
    p11 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 12
    p12 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 13
    p13 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 14
    p14 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)
    time_step <- 15
    p15 <- make_pvalue_histogram_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"), time_step)
    print(time_step)

    pdf(output_file)
    gg <- plot_grid(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,nrow=4,ncol=4,label_size=13)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1)
    print(combined_gg)
    dev.off()

}




visualization_directory = args[1]
distance = args[2]
maf_cutoff = args[3]
normalization_method = args[4]
independent_time_step_eqtl_dir = args[5] 


#############################################
# Make histogram of nominal pvalue distribution for each time step
# Then put each time step on one plot
file_stem <- paste0(independent_time_step_eqtl_dir, "eqtl_prepare_eqtl_distance_",distance,"_maf_cutoff_",maf_cutoff,"_normalization_meth_",normalization_method,"_time_step_")
output_file <- paste0(visualization_directory, "nominal_eqtl_pvalue_histogram_",distance,"_maf_cutoff_",maf_cutoff,"_normalization_meth_",normalization_method, "_histogram_across_time_steps.pdf")
pvalue_histogram_across_time_steps(file_stem, output_file)




###############################################
# Plot distribution of pvalues found in our eqtl data
# But only those variant-gene pairs that are found in:
## a. Nick Banovich's ipsc data
ipsc_file_stem <- paste0(visualization_directory, "ipsc_banovich_comparison_distance_",distance,"_maf_cutoff_",maf_cutoff,"_normalization_meth_",normalization_method,"_time_step_")
ipsc_plot_file <- paste0(visualization_directory, "ipsc_banovich_comparison_distance_",distance,"_maf_cutoff_",maf_cutoff,"_normalization_meth_",normalization_method, "_bar_plot.png")
#eqtl_comparison_to_reference_bar_plot(ipsc_file_stem, ipsc_plot_file)

## b. GTEx Heart left ventricle data
gtex_file_stem <- paste0(visualization_directory, "gtex_v7_heart_left_ventricle_comparison_distance_",distance,"_maf_cutoff_",maf_cutoff,"_normalization_meth_",normalization_method,"_time_step_")
gtex_plot_file <- paste0(visualization_directory, "gtex_v7_heart_left_ventricle_comparison_distance_",distance,"_maf_cutoff_",maf_cutoff,"_normalization_meth_",normalization_method, "_bar_plot.png")
#eqtl_comparison_to_reference_bar_plot(gtex_file_stem, gtex_plot_file)


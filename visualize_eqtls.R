args = commandArgs(trailingOnly=TRUE)
library(edgeR)
library(ggplot2)
library(ggthemes)
library(reshape)
library(cowplot)




# Visualize distribution of best nominal pvalues per gene
best_nominal_pvalue_per_gene_histogram <- function(pvalues, output_file,time_step) {
    # Colors for histogram
    barfill <- "#4271AE"
    barlines <- "#1F3552"
    # Title of plot
    title <- paste0("best nominal eqtl pvalues per gene. T = ", time_step)
    # put data into data frame for plotting
    df <- data.frame(pvalues = pvalues)
    # make plots
    p <- ggplot(df, aes(x = pvalues)) +
        geom_histogram(aes(y = ..count..), colour = barlines, fill = barfill) +
        scale_x_continuous(name = "pvalues") +
        scale_y_continuous(name = "counts") + 
        ggtitle(title)
    ggsave(p, file=output_file,width = 15,height=10.5,units="cm")
}

# Visualize distribution of best corrected pvalues per gene
best_corrected_pvalue_per_gene_histogram <- function(pvalues, output_file,time_step) {
    # Colors for histogram
    barfill <- "#4271AE"
    barlines <- "#1F3552"
    # Title of plot
    title <- paste0("best corrected eqtl pvalues per gene. T = ", time_step)
    # put data into data frame for plotting
    df <- data.frame(pvalues = pvalues)
    # make plots
    p <- ggplot(df, aes(x = pvalues)) +
        geom_histogram(aes(y = ..count..), colour = barlines, fill = barfill) +
        scale_x_continuous(name = "pvalues") +
        scale_y_continuous(name = "counts") + 
        ggtitle(title)
    ggsave(p, file=output_file,width = 15,height=10.5,units="cm")
}


# Visualize distribution nominal pvalues
nominal_pvalues_histogram <- function(pvalues, output_file,time_step) {
    # Colors for histogram
    barfill <- "#4271AE"
    barlines <- "#1F3552"
    # Title of plot
    title <- paste0("nominal eqtl pvalues. T = ", time_step)
    # put data into data frame for plotting
    df <- data.frame(pvalues = pvalues)
    # make plots
    p <- ggplot(df, aes(x = pvalues)) +
        geom_histogram(aes(y = ..count..), colour = barlines, fill = barfill,binwidth = 0.01, boundary = 0) +
        scale_x_continuous(name = "pvalues") +
        scale_y_continuous(name = "counts") + 
        ggtitle(title)
    ggsave(p, file=output_file,width = 15,height=10.5,units="cm")
}


eqtl_results_file = args[1]
all_eqtl_results = args[2]
output_dir = args[3]
time_step = args[4]
# Load in data
eqtl_data <- read.table(eqtl_results_file, header=TRUE)
all_eqtl_nominal_pvalues <- read.table(all_eqtl_results, header=TRUE)



###########################################################
# Visualize distribution of nominal pvalues
#############################################################
output_file <- paste0(output_dir,time_step,"_nominal_pvalues_hist.png")
nominal_pvalues_histogram(all_eqtl_nominal_pvalues$pvalue,output_file,time_step)



###########################################################
# Visualize distribution of best nominal pvalues per gene
#############################################################
output_file <- paste0(output_dir,time_step,"_best_nominal_pvalue_per_gene_hist.png")
#best_nominal_pvalue_per_gene_histogram(eqtl_data$pvalue, output_file,time_step)

###########################################################
# Visualize distribution of best corrected pvalue per gene
############################################################
output_file <- paste0(output_dir,time_step,"_best_corrected_pvalue_per_gene_hist.png")
#best_corrected_pvalue_per_gene_histogram(eqtl_data$BF_pvalue, output_file,time_step)



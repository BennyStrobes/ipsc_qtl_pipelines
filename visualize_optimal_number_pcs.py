import os
import sys
import pdb
import gzip
import random
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import seaborn as sns
sns.set_style("white")




def line_plot_number_of_genes_below_thresh(file_prefix, threshold, output_file):
    num_pcs = range(7)
    num_hits = []
    for num_pc in num_pcs:
        counts = 0
        used_genes = {}
        input_file = file_prefix + '_num_pcs_' + str(num_pc) + 'eqtl_results.txt'
        f = open(input_file)
        head_count = 0
        for line in f:
            line = line.rstrip()
            data = line.split()
            if head_count == 0:
                head_count = head_count + 1
                continue
            gene_id = data[1]
            pvalue = float(data[-1])
            if gene_id not in used_genes and pvalue < threshold:
                counts = counts + 1
                used_genes[gene_id] = 1
        f.close()
        num_hits.append(counts)
    plt.clf()
    fig = plt.figure()
    plt.scatter(num_pcs,num_hits)
    plt.plot(num_pcs,num_hits)
    axes = plt.gca()
    #axes.set_xlim([0,25])
    #axes.set_ylim([-5,700])
    plt.xlabel('Number of PCs')
    plt.ylabel('Number of genes that have association < ' + str(threshold))
    plt.title('time step 0 / threshold = ' + str(threshold))
    plt.savefig(output_file)







file_prefix = sys.argv[1]

threshold = .001
output_file = file_prefix + '_number_genes_below_thresh_' + str(threshold) + '_line_plot.png'
line_plot_number_of_genes_below_thresh(file_prefix, threshold, output_file)

threshold = .0005
output_file = file_prefix + '_number_genes_below_thresh_' + str(threshold) + '_line_plot.png'
line_plot_number_of_genes_below_thresh(file_prefix, threshold, output_file)

threshold = .0001
output_file = file_prefix + '_number_genes_below_thresh_' + str(threshold) + '_line_plot.png'
line_plot_number_of_genes_below_thresh(file_prefix, threshold, output_file)
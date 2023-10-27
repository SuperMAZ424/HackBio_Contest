# Figures from g to j.

# define the plotting function
plot_this <- function(x, y, group, xlab = NULL, ylab = NULL, main = NULL,
                      scatter = FALSE){
  
  # extract the grouping factors
  groups = unique(group)
  cl = rainbow(length(groups))
  # generate the plot
  plot(1, xlab = xlab, ylab = ylab, 
       type = "n", main = main, xlim=c(min(na.omit(x))*0.9,max(na.omit(x))*1.1),
       ylim=c(min(na.omit(y)*0.9),max(na.omit(y))*1.1), las = 1)
  if (scatter == FALSE){
    for (i in 1:length(groups)) {
      lines(x[group == groups[i]], 
            y[group == groups[i]],
            type = "l", col = cl[i])
    }
  }else{
    for (i in 1:length(groups)){
      points(x[group == groups[i]], 
             y[group == groups[i]],pch=16,
             col = cl[i])
    }
    
  }
}



pdf("Figs_g_to_j.pdf", width = 16, height = 9)
par(mfrow = c(2,2))


fig1gdat = read.table("fig_one_g.txt", sep = "\t", header = T)
plot_this(x=fig1gdat$Median_unique_nr_frag, y = fig1gdat$fmx_delta_donor_llk,
          group = fig1gdat$technology, main = "Fig 1g", xlab="Median Unique Fragments",
          ylab="Median Freemuxlet", scatter = T)


fig1hdat = read.table("fig_one_h.txt", sep="\t", header = T)
plot_this(y=fig1hdat$seurat_score, x = fig1hdat$log_median_unique_nr_frag_in_regions,
          group = fig1hdat$technology, main = "Fig 1h", xlab="Median Unique Fragments",
          ylab="Median Seurat Score", scatter = T)


fig1idat = read.table("fig_one_i.txt", sep="\t", header = T)
plot_this(x=fig1idat$Median_tss_enrichment, y = fig1idat$Median_frip,
          group = fig1gdat$technology, main = "Fig 1i", xlab="Median TSS Enrichment",
          ylab="Median FRIP", scatter = T)



##### The first trial to plot Fig. J
# fig1jdat = read.table("fig_one_j.txt", sep="\t", header = T)
# plot_this(x=fig1jdat$top2kdars_median_dar_tss_dist/1000, y = fig1jdat$top2kdars_median_dar_logfc,
#           group = fig1gdat$technology, main = "Fig 1j", ylab="Median Strength",
#           xlab="Median Distance to TSS (kbp)", scatter = T)



#From the previous trial it is obvious that there is an outlier in the data. 
#So, the next step is to plot it after removing the outliers.



#get the data where outliers exist
out = fig1jdat$top2kdars_median_dar_tss_dist

#create new data where outliers are removed by using quantile and IQR
fig1jdat1=fig1jdat[fig1jdat$top2kdars_median_dar_tss_dist <= (1.5*IQR(out)+quantile(out,0.75)),]

# do the same for the y axis data
out2 = fig1jdat1$top2kdars_median_dar_logfc
fig1jdat1 = fig1jdat1[fig1jdat1$top2kdars_median_dar_logfc <=quantile(out2, 0.75)+1.5*IQR(out2),]


#plot the new data
plot_this(x=fig1jdat1$top2kdars_median_dar_tss_dist/1000, y = fig1jdat1$top2kdars_median_dar_logfc,
          group = fig1gdat$technology, main = "Fig 1j", ylab="Median Strength",
          xlab="Median Distance to TSS (kbp)", scatter = T)

dev.off()

# HackBio Contest module 1

data = read.table("fig_One_a_e.dat", sep="\t", header = T)


# extract x axis and y axis data

x = data$depth
y = data$Unique_nr_frag_in_regions
techs = unique(data$tech)

par(mfrow = c(2,3))


plot_this <- function(x,y,group, xlab = NULL, ylab = NULL, main = NULL,
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
  # add the legend to the curve
  #legend("topleft", legend = groups, col = cl, pch = 16, box.lty = 0,
  #       bg = "transparent")

  
  }

pdf("Figs_a_to_f.pdf", width = 16, height =9)
par(mfrow = c(2,3))

plot_this(x,y, data$tech, main = "Fig 1a", 
          ylab = "Unique Fragments in Peaks", xlab = "Sequencing Depth")
plot_this(x, data$TSS_enrichment,data$tech, main="Fig 1b", 
          ylab="Sequencing Depth", xlab = "Unique Fragments in Peaks")
plot_this(x, data$X._unique_nr_frag_in_regions_in_cells, data$tech,
          main="Fig 1c",ylab="Sequencing Efficiency",xlab = "Sequencing Depth")

plot_this(x,data$median_cell_type_pred_score,data$tech,
          main="Fig 1d", xlab="Sequencing Depth",
          ylab = "Median Seurat Score")

plot_this(x, data$fc__B_cell, data$tech, main="Fig 1e", 
          xlab = "Sequencing Depth", ylab="BCell Strength (FC)")

fig1fdat = read.table("fig_one_f.txt", sep = "\t", header = T)

plot_this(x=fig1fdat$Median_Unique_nr_frag_in_regions, y=fig1fdat$Mean_scrublet_doublet_scores_fragments,
          group = fig1fdat$technology, main = "Fig 1f", ylab = "Median Scrublet Score",
          xlab="Median Unique Fragments", scatter = TRUE)

dev.off()


  

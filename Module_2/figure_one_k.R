fig1kdat = read.table("figure_one_k.txt", sep="\t", header = T)
pdf("Fig_K.pdf", width = 16, height = 9)
par(mfrow = c(2,5))

 lay1 = unique(fig1kdat$tech)
 cl = rainbow(length(lay1))
 for (i in 1:length(lay1)){

   dtfm = fig1kdat[fig1kdat$tech == lay1[i], ]
   
   
   bp = boxplot(data = dtfm, seurat_cell_type_pred_score~sample_id, col= cl[i], 
           main = lay1[i],xaxt = "n", ylab= "Seurat Score")

 }
dev.off()
 
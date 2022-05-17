myheatmap <- function(res, contrast, intgroup, ddsKEEP) {
  ### sub: top genes ###
  resD <- res
  resD
  summary(resD)
  resDSort <- resD[order(resD$padj, na.last = NA, decreasing = F),]
  resDSort
  topDESeq2 <- resDSort[1:100,]
  topDESeq2
  topDESeq2 <- topDESeq2[order(topDESeq2$log2FoldChange, na.last = NA, decreasing = T),]
  topDESeq2
  write.csv(topDESeq2, file=paste(resultsdir, contrast, "_cmydas_topDESeq2_100_.csv", sep = ""))
  rld <- vst(ddsKEEP)
  head(assay(rld))
  topgenes <- head(rownames(resDSort),40)
  mat <- assay(rld)[topgenes,]
  (mat <- mat -rowMeans(mat))
  col.pan <- colorpanel(100, "blue","white","red")
  #Non scaled heatmap for topgenes
  heatmap.2(mat, col=col.pan, Rowv=TRUE, scale="none", trace="none")
  #Scaled heatmap for topgenes
  #intgroup <- Sex
  (scaled.mat<-t(scale(t(mat))))
  #intgroup=ddsKEEP$Cadmium
  pheatmap(scaled.mat, fontsize_row = 6, labels_col = intgroup,
           #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
           color = col.pan, 
           scale = "row",
           main = "Heatmap for the top 40 differentially expressed genes",treeheight_col = 25, cluster_cols = T, 
           cluster_rows = F, clustering_distance_cols = "euclidean" )
  heatmap.2(scaled.mat, col=col.pan, Rowv=T, scale="row", #labRow = ddsKEEP$Sex,
            trace="none" ,margins = c(10,8),cexRow=0.5, cexCol=1, keysize=1,labCol = intgroup, 
            main =  "Heatmap for the top 40 DEGs", dendrogram = "both", cex.main = 10)
  dev.off()
  pdf(paste0(resultsdir, contrast, "_VSD_scaled_topgenes_heatmap.pdf"), width = 7, height = 7)
  heatmap.2(scaled.mat, col=col.pan, Rowv=T, scale="row", #labRow = ddsKEEP$Sex,
            trace="none" ,margins = c(10,8),cexRow=0.5, cexCol=1, keysize=1,labCol = intgroup, 
            main =  "Heatmap for the top 40 DEGs", dendrogram = "both", cex.main = 10)
  dev.off()
  graphics.off()
  ###MAPlots
  plotMA(res,ylim=c(-5,5))
  
}

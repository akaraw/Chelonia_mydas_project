mytopgene <- function(res, intgroup, contrast, ddsKEEP) {
  topGene <- rownames(res)[which.min(res$log2FoldChange)]
  mygene <- topGene
  d <- plotCounts(ddsKEEP,gene = mygene, intgroup = intgroup, #,ntgroup, 
                  main = paste( "feature -", contrast, '_', mygene), returnData=TRUE)
  p <- ggplot(d, aes(x= d[,2], y=count)) + 
    geom_boxplot(colour = "red", fill = "orange", alpha = 0.2, 
                 outlier.colour="black", outlier.shape=8, outlier.size=2, notch=F, notchwidth = 0.5) + 
    #geom_dotplot(binwidth = 50, binaxis='y', stackdir='center', dotsize=1)
    geom_point(position=position_jitter(w=0.1,h=0), colour = 'purple', size = 1) + 
    scale_y_log10(breaks=c(25,100,400)) + 
    theme(
      #panel background elements
      panel.background = element_rect(
        fill = "grey90",
        colour = "black",
        size = 1,
      ),
      legend.position= "bottom",
      plot.title = element_text(color="black", size=12, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=11, face="bold"),
      axis.title.y = element_text(color="#993333", size=11, face="bold")
    ) + 
    ggtitle(paste0("Condition: ", contrast, " gene - ", mygene)) + xlab(contrast) + 
    ylab("Noramlized gene count") +
    labs(fill = d$intgroup) +
    stat_summary(fun=mean, geom="point", shape=23, size=4) + scale_color_grey() +
    scale_fill_manual(values=c("#999999", "#E69F00"))
  print(p)
  ggsave(p, file=paste0(resultsdir, contrast,"_", mygene,".png", sep = ""), width = 7, height = 7, units = "cm")
  
}

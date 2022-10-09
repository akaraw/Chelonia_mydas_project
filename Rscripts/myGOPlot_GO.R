library(GOplot)

myGOPlot <- function(res, species, contrast, goBPRes, resultsdir) {

#CM for chelina mydas
goID <- goBPRes$category
len=length(goID)
CM <- data.frame(Category=rep("BP", len), ID=goID, Term=goBPRes$term, Genes=NA)
CM
head(CM <- merge(CM, GOdf, by.x= "ID", by.y= "GO"))
CM
CM$Genes <- CM$Genes.y
CM <- dplyr::select(CM, -c("Genes.x", "Genes.y"))  
head(CM)
goBPRes$over_represented_pvalue
CM$adj_pval <- goBPRes$over_represented_pvalue

res
resGOplot <- as.data.frame(res)
ID <- data.frame(ID=rownames(resGOplot))
resGOplot <- cbind(ID, resGOplot)
head(resGOplot)
resSig
resGOplot <- dplyr::select(resGOplot, ID, log2FoldChange, baseMean, lfcSE, pvalue, padj, stat)
rownames(resGOplot) <- NULL
colnames(resGOplot) <- c("ID","logFC","AveExpr","t","P.Value", "adj.P.Val", "B")
head(resGOplot)
head(CM)

#resGOplot[resGOplot$ID == "PRPF19",]
#start GOPlot
circ <- circle_dat(CM, resGOplot)
circ
keep_circ <- !duplicated(circ$ID)
zscore <- circ[keep_circ,]
zscore
#zscore <- dplyr::select(zscore, -c("genes", "logFC", "count", "adj_pval", "category"))
GOBar(subset(circ, category == 'BP'))
unique(circ$term)
process <- unique(circ$term)

circ

pdf(file = paste0(resultsdir, "/", contrast, species, "_GO_circ.pdf"), width = 16, height = 16)
GOCircle(circ, nsub = process)
dev.off()


genes <- dplyr::select(resGOplot, ID, logFC)

# Create the plot
process <- unique(circ$term)
chord <- chord_dat(data = circ, genes = genes, process = process)
(chord)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)

pdf(file = paste0(resultsdir, "/", contrast, species, "_GO_plot.pdf"), width = 16, height = 16)
fig <- GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
print(fig)
dev.off()
}

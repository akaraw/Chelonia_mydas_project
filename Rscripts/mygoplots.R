mygoplot <- function(contrast, goResults, myspecies, resultsdir, geneID2GO, method = "Wallenius", 
                     repcnt = 2000, pwf){
  graphics.off()
  pdf(paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".pdf"), width = 14, height = 14)
  fig <- goResults %>% 
    top_n(30, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
    ggplot(aes(x=hitsPerc, y=term)) + 
    geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
    theme(text = element_text(size = 12, family = "serif")) +
    labs(x="Hits (%)", y=NULL, colour="p value", size="Count")
  print(fig)
  dev.off()
  graphics.off()
  goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)
  goResults = goResults[goResults$under_represented_pvalue<0.05,]
  if (nrow(goResults) > 0) {
    
    pdf(paste0(resultsdir, "/", contrast, "_enriched_go_", "under_represented_", 
               myspecies, ".pdf"))
    
    fig <- goResults %>% 
      top_n(30, wt=-under_represented_pvalue) %>% 
      mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
      ggplot(aes(x=hitsPerc, y=term)) + 
      geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
      theme(text = element_text(size = 12, family = "serif")) +
      labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 
    print(fig)
    dev.off()
    
  }
}

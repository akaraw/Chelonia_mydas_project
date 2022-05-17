library(rtracklayer)
library(GenomicFeatures)

gtf <- "C:/Users/kar131/OneDrive - CSIRO/15.Chritable_green_turtle/NCBI/GCF_015237465.2_rCheMyd1.pri.v2_genomic.gtf"
gtf <- rtracklayer::readGFF(gtf)
featuret <- read.table("C:/Users/kar131/OneDrive - CSIRO/15.Chritable_green_turtle/NCBI/GCF_015237465.2_rCheMyd1.pri.v2_feature_table.txt",
                      sep = "\t", fill = T, header = T)

featuret <- featuret[featuret$feature == "mRNA",]
df <- featuret %>%
  group_by(symbol) %>%
  filter(product_length == max(product_length)) %>%
  arrange(related_accession, symbol, name)

list <- df$related_accession
write(list, "C:/Users/kar131/OneDrive - CSIRO/15.Chritable_green_turtle/NCBI/longest_isoforms.txt")
df1 <- read.table("C:/Users/kar131/OneDrive - CSIRO/15.Chritable_green_turtle/NCBI/cmydas_go",
                  fill = T, sep = "\t")
df1$V14 <- gsub("[|]",",",df1$V14)
df1 <- dplyr::select(df1, V1, V14)

merged <- merge(df, df1, by.x = "related_accession", by.y = "V1")
merged <- dplyr::select(merged, symbol, V14)
merged <- merged[which(merged$V14 != "-"),]
(gene2go <- merged[which(merged$V14 != ""),])
gene2go <- separate_rows(gene2go, V14, sep = ",") %>%
  mutate(symbol) %>%
  as.data.frame()

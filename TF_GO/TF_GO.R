library(clusterProfiler)
load("TF_GO.Rdata")
gene.df <- bitr(genelist, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(
  gene          = gene.df$ENTREZID,
  keyType = "ENTREZID",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory = 10)
dotplot(ego, showCategory = 10)
emapplot(ego, showCategory = 30)




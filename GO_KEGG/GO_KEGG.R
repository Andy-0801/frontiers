library(org.Hs.eg.db)
library(clusterProfiler)
keytypes(org.Hs.eg.db) 
library(enrichplot)
load("GO.Rdata")
genelist = degenes$V1
head(genelist)
gene.df <- bitr(genelist, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
gene.df
go_BP <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05,keyType = 'SYMBOL')


dotplot(go_BP, showCategory=20)
write.csv(go_BP,"GO-enrich.csv",row.names =F)


KEGG <- enrichKEGG(gene = gene.df$ENTREZID, 
           organism = "hsa",
           keyType = "kegg",
           pvalueCutoff = 0.05, 
           pAdjustMethod = "BH",
           minGSSize = 10, 
           maxGSSize = 500,
           qvalueCutoff = 0.1, 
           use_internal_data = FALSE)
dotplot(KEGG, showCategory=15)
barplot(KEGG, showCategory=15)
write.csv(summary(kk),"KEGG-enrich.csv",row.names =F)


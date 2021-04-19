
load("tur_module.Rdata")
keggSig = KEGG[KEGG$PValue < 0.05,]
library(tidyr)
keggSig = separate(keggSig, Term, sep = ":",
                   into = c("ID", "Term"))
library(ggplot2)
ggplot(keggSig,aes(x=Fold.Enrichment,y=Term)) + 
  geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_colour_gradient(low="green",high="red")+
  labs(
    color=expression(-log[10](P.value)),
    size="Gene number",
    x="Fold enrichment",
    title="Pathway enrichment"
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )

ggsave('tur_david.pdf',width = 7,height = 4)

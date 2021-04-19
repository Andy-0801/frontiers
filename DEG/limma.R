library(limma)
library(dplyr)
load("DEG.Rdata")
group <- phenodata[, 1]
design <- model.matrix(~0+factor(group))
rownames(design) <- colnames(datExpr)
colnames(design) <- levels(factor(group))
head(design)
contrast.matrix<-makeContrasts("Tumor-Normal",
                               levels = design)

deg = function(datExpr,design,contrast.matrix){
  ##step1
  fit <- lmFit(datExpr,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  head(nrDEG)
  return(nrDEG)
}

deg = deg(datExpr,design,contrast.matrix)
head(deg)
dim(deg)
write.csv(deg,file = "gene2.csv")
save(deg,file = 'deg.Rdata')
library(ggplot2)
diff <- deg
logFC <-diff$logFC
adj <- diff$ adj.P.Val
data <- data.frame(logFC=diff$logFC,padj=diff$ adj.P.Val)
data$sig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < 1)& data$logFC > -1] <- "no"
 data$sig[data$padj <= 0.05 & data$logFC >= 1] <- "up"
data$sig[data$padj <= 0.05 & data$logFC <= -1] <- "down"
x_lim <- max(logFC,-logFC)
library(ggplot2)
library(RColorBrewer)
theme_set(theme_bw())
p <- ggplot(data,aes(logFC,-1*log10(padj),
                     color = sig))+geom_point()+
  xlim(-3,6) +  labs(x="log2(FoldChange)",y="-log10(P.Value)")+ylim(0,125)
p <- p + scale_color_manual(values =c('blue','black','red'))+ 
  geom_hline(yintercept=-log10(0.01),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))
p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p
print(p)

getwd()

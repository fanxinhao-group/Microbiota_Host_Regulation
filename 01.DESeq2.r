library(tidyverse)
library(DESeq2)
library(dplyr)

counts<-read.csv("gene_count_matrix.csv",header = T)
colDatas <- read.csv("conditions.csv",header = T)
rownames(counts)<-counts[,1]
counts<-counts[,-1]

mycounts <- counts[,1:12]
colData <- colDatas[1:12,]
condition <- factor(c(rep("control",6),rep("treat",6)), levels = c("control","treat"))
condition
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)

write.csv(res,file="DESeq2_results.csv")
table(res$padj<0.05)

diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
diff_gene_deseq2 <- data.frame(diff_gene_deseq2)
diff_gene_deseq2 <- diff_gene_deseq2 %>% mutate(type= if_else(.$log2FoldChange > 0, 'Up', 'Down'))
diff_gene_deseq2 <- data.frame(diff_gene_deseq2, geneID = rownames(diff_gene_deseq2))
diff_gene_deseq2 <- diff_gene_deseq2 %>% separate(col = geneID, into = c("GeneID", "Symbol"), sep = "\\|", remove = TRUE)
diff_gene_deseq2 <- diff_gene_deseq2 %>% mutate(Symbol = replace(Symbol, is.na(Symbol), "-"))

dim(diff_gene_deseq2)

write.csv(diff_gene_deseq2,file= "DEG.csv",quote = F,row.names = T)
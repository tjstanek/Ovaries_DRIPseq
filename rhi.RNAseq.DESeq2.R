library("DESeq2")

# Processing RNAseq
sampleInfo <- read.csv("rhi.RNAseq.te.sampleInfo.csv")
countdata <- read.table("rhi.RNAseq.te.count.txt", row.names=1)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                       colData = sampleInfo,
                       design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- relevel(dds$condition, ref = "Ctrl")
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)

write.table(res, "rhi_RNAseq_te_DESeq.txt", sep = "\t", row.names = TRUE, quote = FALSE, col.names = TRUE)
library("DESeq2")

##rlog-transform htseq-rev counts & TE counts combined
RNA <- read.table("rhi.RNAseq.htseq.te.combined.mat", header = FALSE, stringsAsFactors = FALSE)
keep <- rowSums(counts(RNA)) > 0
keep <- data.matrix(keep)
rld <- rlog(keep, blind=TRUE)
head(rld)
write.table(rld, "rhi.RNAseq.htseq.te.rlog_normalized.txt", sep = "\t", row.names = TRUE, quote = FALSE, col.names = FALSE)

##rlog-transform k-seek counts
DRIP <- read.table("rhi.DRIPseq.kseek.compiled.mat", header = FALSE, stringsAsFactors = FALSE)
head(DRIP, 3)
DRIP <- data.matrix(DRIP)
head(DRIP, 3)
rld <- rlog(DRIP, blind=TRUE)
head(rld)
write.table(rld, "rhi.DRIPseq.kseek.rlog_normalized.txt", sep = "\t", row.names = TRUE, quote = FALSE, col.names = FALSE)

RNA <- read.table("rhi.RNAseq.kseek.compiled.mat", header = FALSE, stringsAsFactors = FALSE)
head(RNA, 3)
RNA <- data.matrix(RNA)
head(RNA, 3)
rld <- rlog(RNA, blind=TRUE)
head(rld)
write.table(rld, "rhi.RNAseq.kseek.rlog_normalized.txt", sep = "\t", row.names = TRUE, quote = FALSE, col.names = FALSE)



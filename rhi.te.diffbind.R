library(DiffBind)

#Diffbind for bedtools-intersected sets(IP/Inp&IP/RH+) IDR peaks (TE Only)
idr.T1 <- read.csv("rhi-dripseq-teall-samples-idr-ipinprh.csv")
idr.T1 <- dba(sampleSheet=idr.T1)
idr.T1 <- dba.count(idr.T1, bUseSummarizeOverlaps=TRUE,summits=FALSE)
dba.plotPCA(idr.T1,  attributes=DBA_TREATMENT, label=DBA_ID)
idr.PCA <- plot(idr.T1)
png("rhi.DRIP.teonly.idr.ipinprh.PCA.png", width=400, height=400)
plot(idr.PCA)
dev.off()

th=0.05
idr.T1 <- dba.normalize(idr.T1, normalize=DBA_NORM_LIB,library=DBA_LIBSIZE_FULL)
idr.T1.contrast <- dba.contrast(idr.T1, categories=DBA_TREATMENT, minMembers=2)
idr.T1.contrast
idr.T1.contrast <- dba.analyze(idr.T1.contrast, method=DBA_ALL_METHODS)
dba.show(idr.T1.contrast, bContrasts=TRUE)
idr.T1.out <- dba.report(idr.T1.contrast, method=DBA_DESEQ2)
idr.T1.df <- as.data.frame(idr.T1.out)
write.table(idr.T1.df, file="results/rhi_DRIPseq_teonly_ipinprh_idr_bytreatment_deseq_th.05.txt", sep="\t", quote=F, row.names=F)
idr.T1.out <- dba.report(idr.T1.contrast, method=DBA_DESEQ2, th=1)
idr.T1.df <- as.data.frame(idr.T1.out)
write.table(idr.T1.df, file="results/rhi_DRIPseq_teall_ipinprh_idr_bytreatment_deseq.th.all.txt", sep="\t", quote=F, row.names=F)
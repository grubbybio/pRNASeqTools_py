options(warn=-1)
args = commandArgs(trailingOnly = T)
norm = args[1]
pvalueoo <- as.numeric(args[2])
fdroo <- as.numeric(args[3])
foldchangeoo <- as.numeric(args[4])
args[-c(1:6)] -> args
args[c(TRUE, FALSE)] -> genotype
args[c(FALSE, TRUE)] -> replicates
message("Parameters OK! loading...")
suppressMessages(library(DESeq2))
suppressMessages(library(NMF))
library(RColorBrewer)
library(pheatmap)
nmf.options(grid.patch=TRUE)
p <- sum(as.numeric(replicates))
data.frame(rep("x",p), stringsAsFactors = FALSE) -> b
colnames(b) <- "genotype"
p <- 0
for (n in 1:length(genotype)){
  for (k in 1:replicates[n]){
    p = p + 1
    paste(genotype[n],"_",k,sep="") -> rownames(b)[p]
    genotype[n] -> b[p,]
    assign(rownames(b)[p], read.table(paste(genotype[n],"_",k,".txt",sep = ""), as.is = T, row.names = 1,header = T))
  }
}
message("Loading completed.")

# Extract gene lengths (bp) — same across all replicates, take from first file
first_rep_data <- eval(parse(text=rownames(b)[1]))
gene_lengths <- first_rep_data$Length
names(gene_lengths) <- rownames(first_rep_data)

eval(parse(text=rownames(b)[1]))[,1] -> a
for(q in 2:p){
  cbind(a,eval(parse(text=rownames(b)[q]))[,1]) -> a
}
rownames(a) <- rownames(eval(parse(text=rownames(b)[1])))
colnames(a) <- rownames(b)
apply(a,2,sum)/1000000 -> rtotal
a[apply(a,1,mean)>=1,] -> a
DESeqDataSetFromMatrix(countData = a, colData = b, design = ~ genotype) -> dds
if(norm == "RPM"){
  sizeFactors(dds) <- rtotal
}
DESeq(dds) -> dds
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("genotype")])
rownames(b)[1:p] -> rownames(df)
colnames(df) <- "Genotype"
rld <- rlog(dds)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pdf(file="DEG_overview.pdf",6,6)
aheatmap(log2.norm.counts, Rowv=NA, Colv=NA, annCol=df)
message("Top 1000 completed!")
plotPCA(rld,ntop=1000,intgroup=c("genotype"))
message("PCA completed!")
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,show_colnames = F,col=colors)
message("Dist completed!")
dev.off()
for(j in 2:length(genotype)){
  results(dds,contrast = c("genotype",genotype[j],genotype[1])) -> res
  out <- as.data.frame(res)

  # Add per-sample TPM columns: TPM = (counts / len_kb) / sum(counts / len_kb) * 1e6
  gene_len_kb <- gene_lengths[rownames(out)] / 1000               # bp → kb
  rpkm_mat <- as.matrix(counts(dds)) / gene_len_kb                 # RPK (counts / kb)
  tpm_mat <- sweep(rpkm_mat, 2, colSums(rpkm_mat) / 1e6, "/")      # TPM normalize
  colnames(tpm_mat) <- paste0(colnames(counts(dds)), "_TPM")
  out <- cbind(out, as.data.frame(tpm_mat))
  if(fdroo < 1){
    subset(out, padj < fdroo & log2FoldChange >= log2(foldchangeoo)) -> resup
    subset(out, padj < fdroo & log2FoldChange <= -log2(foldchangeoo)) -> resdown
  }else{
    subset(out, pvalue < pvalueoo & log2FoldChange >= log2(foldchangeoo)) -> resup
    subset(out, pvalue < pvalueoo & log2FoldChange <= -log2(foldchangeoo)) -> resdown
  }
  write.csv(out,paste(genotype[j],"vs",genotype[1],"total","csv",sep="."),quote=F)
  write.csv(resup,paste(genotype[j],"vs",genotype[1],"total","upregulated","csv",sep="."),quote=F)
  write.csv(resdown,paste(genotype[j],"vs",genotype[1],"total","downregulated","csv",sep="."),quote=F)

  # ── Volcano plot ──────────────────────────────────────────────────────────
  volcano_pdf <- paste(genotype[j],"vs",genotype[1],"volcano","pdf",sep=".")
  message(paste("Generating volcano plot:", volcano_pdf))
  pval_col <- if(fdroo < 1) "padj" else "pvalue"
  # Replace NA padj with 1 for plotting
  plot_df <- out[!is.na(out$log2FoldChange), ]
  plot_pval <- plot_df[[pval_col]]
  plot_pval[is.na(plot_pval)] <- 1
  plot_fc   <- plot_df$log2FoldChange
  sig_up   <- plot_fc >=  log2(foldchangeoo) & plot_pval < (if(fdroo < 1) fdroo else pvalueoo)
  sig_down <- plot_fc <= -log2(foldchangeoo) & plot_pval < (if(fdroo < 1) fdroo else pvalueoo)
  pdf(file=volcano_pdf, width=7, height=6)
  plot(plot_fc, -log10(plot_pval + 1e-300),
       pch=16, cex=0.5, col="grey60",
       xlab=expression(log[2]~Fold~Change),
       ylab=expression(-log[10]~Pvalue),
       main=paste(genotype[j], "vs", genotype[1]))
  points(plot_fc[sig_up],   -log10(plot_pval[sig_up]   + 1e-300),
         pch=16, cex=0.6, col="#E41A1C")
  points(plot_fc[sig_down], -log10(plot_pval[sig_down] + 1e-300),
         pch=16, cex=0.6, col="#377EB8")
  abline(v=c(-log2(foldchangeoo), log2(foldchangeoo)),
         h=-log10(if(fdroo < 1) fdroo else pvalueoo),
         lty=2, col="grey30")
  legend("topright",
         legend=c(paste0("Up (", sum(sig_up), ")"),
                  paste0("Down (", sum(sig_down), ")"),
                  paste0("NS (", sum(!sig_up & !sig_down), ")")),
         col=c("#E41A1C", "#377EB8", "grey60"), pch=16, cex=0.8)
  dev.off()
  message(paste("Volcano plot saved:", volcano_pdf))
}

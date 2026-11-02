#!/usr/bin/env Rscript
# Independent coefficient-based sensitivity; does not replace primary edgeR discovery.
suppressPackageStartupMessages({library(readxl); library(edgeR); library(limma)})
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)) args[1] else "/Volumes/wang/paper14/ais_paraspinal_ecology/raw/bulk/GSE254300_readcount.xlsx"
out <- "results/tissue/model_dependence"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- as.data.frame(read_excel(input))
x <- as.matrix(d[,-1]); storage.mode(x) <- "numeric"
rownames(x) <- make.unique(as.character(d[[1]]))
id <- colnames(x)
stopifnot(all(is.finite(x)),all(x>=0),all(x==floor(x)),!anyDuplicated(id))
donor <- factor(sub("[AT]$","",id))
side <- as.integer(grepl("T$",id))
ais <- as.integer(grepl("^IS",id))
stopifnot(all(table(donor,side)==1))
design <- cbind(model.matrix(~0+donor),side=side,interaction=ais*side)
y <- DGEList(x)
universe <- read.delim("results/tissue/gse254300_bulk/retained_gene_universe.tsv")[[1]]
stopifnot(all(universe %in% rownames(y)))
y <- calcNormFactors(y[universe,,keep.lib.sizes=FALSE])
v <- voom(y,design,plot=FALSE)
C <- cbind(mean=c(rep(0,nlevels(donor)),-1,-.5),
           interaction=c(rep(0,nlevels(donor)),0,-1))
f <- eBayes(contrasts.fit(lmFit(v,design),C),robust=TRUE)
se <- f$stdev.unscaled*sqrt(f$s2.post)
z <- data.frame(gene_symbol=rownames(f),mean_logFC=f$coefficients[,1],
                mean_FDR=p.adjust(f$p.value[,1],"BH"),
                interaction_logFC=f$coefficients[,2],interaction_SE=se[,2],
                interaction_FDR=p.adjust(f$p.value[,2],"BH"))
legacy <- read.delim(gzfile("results/tissue/gse254300_bulk/cross_etiology_equivalence/gene_cross_etiology_equivalence.tsv.gz"))
legacy <- legacy[match(z$gene_symbol,legacy$gene_symbol),]
stopifnot(identical(z$gene_symbol,legacy$gene_symbol))
summary <- lapply(c(.25,.5,1),function(m){
 p <- pmax(pt((z$interaction_logFC+m)/z$interaction_SE,f$df.total,lower.tail=FALSE),
           pt((z$interaction_logFC-m)/z$interaction_SE,f$df.total))
 q <- p.adjust(p,"BH")
 z[[paste0("equivalence_FDR_",m)]] <<- q
 data.frame(margin=m,voom_equivalent_all=sum(q<.05),
            voom_equivalent_primary_mean_hits=sum(q<.05 & legacy$shared_fdr05),
            legacy_strict_confirmed=sum(q<.05 & legacy$strictly_conserved_margin_0p50))
})
write.table(z,file.path(out,"voom_gene_sensitivity.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,summary),file.path(out,"equivalence_model_comparison.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
metrics <- data.frame(metric=c("genes","voom_mean_FDR05","primary_hits_voom_FDR05","primary_hits_same_sign","mean_effect_spearman","voom_interaction_FDR10"),
 value=c(nrow(z),sum(z$mean_FDR<.05),sum(z$mean_FDR<.05 & legacy$shared_fdr05),
 sum(sign(z$mean_logFC)==sign(legacy$shared_logFC) & legacy$shared_fdr05),
 cor(z$mean_logFC,legacy$shared_logFC,method="spearman"),sum(z$interaction_FDR<.1)))
write.table(metrics,file.path(out,"model_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(metrics); print(do.call(rbind,summary))

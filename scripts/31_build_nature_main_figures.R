#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})
set.seed(14016)

root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else ".")
outdir <- file.path(root, "results", "manuscript_figures_nature")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ink <- "#24292F"; muted <- "#6B7280"; grid <- "#E5E7EB"
blue <- "#356C9B"; teal <- "#2A8C82"; ochre <- "#C78A2C"
red <- "#B64E4A"; pale_blue <- "#BFD6E6"; pale <- "#F4F6F8"
class_cols <- c("#2D7F72", "#6E9F8D", "#C4A46B", "#A7A9AC", "#7B5B59")

theme_nc <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = ink), axis.text = element_text(colour = ink),
      axis.title = element_text(colour = ink), plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(colour = muted, size = base_size),
      plot.caption = element_text(colour = muted, hjust = 0, size = base_size - 1),
      legend.title = element_text(face = "bold"), legend.key.height = unit(3.2, "mm"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

short_path <- function(x) {
  x <- sub("^(HALLMARK|REACTOME)__", "", x)
  x <- sub("^(HALLMARK|REACTOME)_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  x <- paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))
  x <- gsub("complex iii", "Complex III", x, ignore.case=TRUE)
  x <- gsub("complex iv", "Complex IV", x, ignore.case=TRUE)
  x <- gsub("complex i ", "Complex I ", x, ignore.case=TRUE)
  x <- gsub("tca", "TCA", x, fixed=TRUE)
  x <- gsub("cdh1", "CDH1", x, fixed=TRUE)
  x
}

save_figure <- function(p, stem, width = 7.20, height = 7.00) {
  ggsave(file.path(outdir, paste0(stem, ".pdf")), p, width = width, height = height,
         device = grDevices::pdf, bg = "white", useDingbats = FALSE)
  ggsave(file.path(outdir, paste0(stem, ".svg")), p, width = width, height = height,
         device = svglite::svglite, bg = "white")
  ggsave(file.path(outdir, paste0(stem, ".png")), p, width = width, height = height,
         device = ragg::agg_png, dpi = 600, bg = "white")
  ggsave(file.path(outdir, paste0(stem, ".tiff")), p, width = width, height = height,
         device = ragg::agg_tiff, dpi = 600, compression = "lzw", bg = "white")
}

read_tsv <- function(...) read.delim(file.path(root, ...), check.names = FALSE, stringsAsFactors = FALSE)

# ---- Figure 1: design, estimands and the failure of list subtraction ----
gene <- read_tsv("results", "integration", "orthogonal_axes", "gene_evidence_axis_map.tsv")
gene <- gene[complete.cases(gene[, c("ais_logFC", "cs_logFC", "shared_FDR", "ais_FDR", "cs_FDR", "interaction_FDR")]), ]
# Annotation joins may repeat a gene; plot the primary gene-level unit once.
gene_fields <- c("gene_symbol","shared_logFC","shared_FDR","ais_logFC","ais_FDR","cs_logFC","cs_FDR","interaction_FDR")
stopifnot(nrow(unique(gene[,gene_fields]))==length(unique(gene$gene_symbol)))
gene <- gene[!duplicated(gene$gene_symbol),]
gene$shared <- !is.na(gene$shared_FDR) & gene$shared_FDR < 0.05
gene$separate <- ifelse(gene$ais_FDR < 0.05 & gene$cs_FDR < 0.05, "Both separate contrasts",
                 ifelse(gene$ais_FDR < 0.05, "AIS only",
                 ifelse(gene$cs_FDR < 0.05, "Congenital only", "Neither")))

# Figure 1 retains only the cross-etiology effect scatter.
lab <- gene[gene$shared,]
lab <- lab[order(lab$shared_FDR),][1:8,]
p1b <- ggplot(gene, aes(ais_logFC, cs_logFC)) +
  geom_hline(yintercept=0, colour=grid, linewidth=.35) + geom_vline(xintercept=0, colour=grid, linewidth=.35) +
  geom_abline(slope=1, intercept=0, linetype="22", colour=muted, linewidth=.4) +
  geom_point(aes(colour=shared), size=.7, alpha=.58) +
  geom_text_repel(data=lab, aes(label=gene_symbol), size=2.2, family="sans", colour=ink,
                  min.segment.length=0, segment.color="#B9BEC3", box.padding=.2, max.overlaps=Inf, seed=14016) +
  scale_colour_manual(values=c(`FALSE`="#C9CDD1",`TRUE`=blue), labels=c("Other genes","Mean-side FDR < 0.05"), name=NULL) +
  coord_equal() + theme_nc() + theme(legend.position="inside",legend.position.inside=c(.23,.88),
    legend.background=element_rect(fill=alpha("white",.9),colour=NA),legend.key.width=unit(3,"mm")) +
  labs(x="AIS side effect (log2 fold change)",y="Congenital-scoliosis side effect (log2 fold change)")

voom <- read_tsv("results","tissue","model_dependence","voom_gene_sensitivity.tsv")
model_comparison <- merge(gene,voom,by="gene_symbol",sort=FALSE)
stopifnot(nrow(model_comparison)==nrow(gene),!anyDuplicated(model_comparison$gene_symbol))
model_rho <- cor(model_comparison$shared_logFC,model_comparison$mean_logFC,method="spearman")
p1_model <- ggplot(model_comparison,aes(shared_logFC,mean_logFC))+
  geom_abline(slope=1,intercept=0,colour=muted,linetype="22",linewidth=.35)+
  geom_point(colour="#BCC4CB",size=.45,alpha=.35)+
  geom_point(data=model_comparison[model_comparison$shared,],colour=blue,size=.65,alpha=.7)+
  annotate("text",x=-Inf,y=Inf,label=sprintf("Spearman rho = %.3f",model_rho),
           hjust=-.06,vjust=1.4,size=2.55,colour=ink)+
  coord_equal()+theme_nc()+
  labs(x="edgeR mean-side effect (log2)",y="voom mean-side effect (log2)")
primary_voom <- model_comparison[model_comparison$shared,]
retained <- sum(primary_voom$mean_FDR<.05)
p1_retention <- ggplot(primary_voom,aes(mean_FDR))+
  stat_ecdf(geom="step",colour=blue,linewidth=.65)+
  geom_vline(xintercept=.05,linetype="22",colour=ochre,linewidth=.45)+
  annotate("point",x=.05,y=retained/nrow(primary_voom),size=2.1,colour=ochre)+
  annotate("text",x=.07,y=.77,label=sprintf("%d / %d retain\nvoom FDR < 0.05",retained,nrow(primary_voom)),
           hjust=0,size=2.55,lineheight=1.05,colour=ink)+
  scale_x_continuous(limits=c(0,1),breaks=c(.05,.25,.5,1))+
  scale_y_continuous(limits=c(0,1),labels=percent_format(accuracy=1))+
  theme_nc()+labs(x="voom BH-adjusted P value",y="Cumulative fraction of\nedgeR discoveries")
fig1 <- (p1b + theme(legend.position="bottom") | (p1_model / p1_retention))+
  plot_layout(widths=c(1.35,1))+
  plot_annotation(tag_levels="a",theme=theme(plot.tag=element_text(family="sans",face="bold",size=10,colour=ink)))
save_figure(fig1,"Figure1_orthogonal_design",height=4.85)
write.table(model_comparison,file.path(outdir,"Figure1_model_comparison.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

# ---- Figure 2: equivalence and donor stability ----
eq <- read.delim(gzfile(file.path(root,"results/tissue/gse254300_bulk/cross_etiology_equivalence/gene_cross_etiology_equivalence.tsv.gz")),check.names=FALSE)
eqs <- read_tsv("results","tissue","gse254300_bulk","cross_etiology_equivalence","equivalence_margin_sensitivity.tsv")
lodo <- read_tsv("results","tissue","gse254300_bulk","leave_one_donor_out_correlations.tsv")
eq_hit <- eq[eq$shared_fdr05,]
eq_hit$class <- ifelse(eq_hit$strictly_conserved_margin_0p50,"Equivalent at ±0.50","Equivalence inconclusive")
p2a <- ggplot(eq_hit,aes(ais_logFC,cs_logFC)) +
  annotate("rect",xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,fill="white") +
  geom_abline(slope=1,intercept=c(-.5,.5),colour="#B9BEC3",linetype="22",linewidth=.35) +
  geom_abline(slope=1,intercept=0,colour=muted,linewidth=.4) +
  geom_point(aes(colour=class),size=.9,alpha=.72) + coord_equal() +
  scale_colour_manual(values=c("Equivalence inconclusive"=pale_blue,"Equivalent at ±0.50"=red),name=NULL) +
  theme_nc()+theme(legend.position="bottom") +
  labs(x="AIS side effect (log2 fold change)",y="Congenital-scoliosis side effect (log2 fold change)")
eqs$margin=factor(eqs$equivalence_margin_log2,levels=eqs$equivalence_margin_log2,labels=paste0("±",eqs$equivalence_margin_log2))
n_shared <- sum(eq$shared_fdr05, na.rm=TRUE)
stopifnot(length(unique(eqs$shared_hits)) == 1, unique(eqs$shared_hits) == n_shared)
voom_eq <- read_tsv("results","tissue","model_dependence","equivalence_model_comparison.tsv")
eq_plot <- rbind(
  data.frame(margin=eqs$equivalence_margin_log2,n=eqs$shared_hits_equivalent_fdr05,model="Approximate QL"),
  data.frame(margin=voom_eq$margin,n=voom_eq$voom_equivalent_primary_mean_hits,model="voom"))
eq_plot$fraction <- eq_plot$n/n_shared
p2b <- ggplot(eq_plot,aes(margin,fraction,colour=model,shape=model,group=model))+
  geom_line(linewidth=.6)+geom_point(size=2.1,fill="white",stroke=.7)+
  geom_text(data=eq_plot[eq_plot$model=="Approximate QL",],aes(label=n),
            vjust=-1,size=2.5,show.legend=FALSE)+
  geom_text(data=eq_plot[eq_plot$model=="voom",],aes(label=n),
            vjust=1.7,size=2.5,show.legend=FALSE)+
  scale_colour_manual(values=c("Approximate QL"=blue,voom=ochre),name=NULL)+
  scale_shape_manual(values=c("Approximate QL"=21,voom=24),name=NULL)+
  scale_x_continuous(breaks=c(.25,.5,1),labels=c("±0.25","±0.50","±1.00"))+
  scale_y_continuous(labels=percent_format(accuracy=1),limits=c(-.12,1),breaks=c(0,.25,.5,.75,1))+
  theme_nc()+theme(legend.position="bottom",legend.text=element_text(size=6.3))+
  labs(x="Interaction margin (log2)",y="Fraction of mean-side discoveries")
lodo$etiology <- ifelse(grepl("^CS_", lodo$omitted_donor), "Congenital", "AIS")
lodo$donor_label <- gsub("_", " ", lodo$omitted_donor, fixed=TRUE)
p2c <- ggplot(lodo,aes(reorder(donor_label,spearman_with_full),spearman_with_full,colour=etiology))+
  geom_segment(aes(xend=donor_label,y=.95,yend=spearman_with_full),linewidth=.55)+geom_point(size=2.2)+
  coord_flip()+scale_colour_manual(values=c(AIS=teal,Congenital=ochre),labels=c(AIS="AIS donor",Congenital="Congenital donor"),name=NULL)+
  scale_y_continuous(limits=c(.95,.975),breaks=c(.95,.96,.97),expand=expansion(mult=c(0,.03)))+
  theme_nc()+theme(legend.position="bottom",plot.margin=margin(4,5,4,10))+
  labs(x="Omitted donor",y="Spearman correlation with full-data effects")
fig2 <- (p2a | (p2b / p2c)) + plot_layout(widths=c(1.12,.88)) +
  plot_annotation(tag_levels="a",theme=theme(plot.tag=element_text(family="sans",face="bold",size=10,colour=ink)))
save_figure(fig2,"Figure2_equivalence_robustness",height=5.75)

# ---- Figure 3: pathway architecture and counterfactual separation ----
cam_shared <- read_tsv("results","tissue","gse254300_bulk","camera_shared_concave_minus_convex.tsv")
cam_ais <- read_tsv("results","tissue","gse254300_bulk","camera_AIS_concave_minus_convex.tsv")
cam_cs <- read_tsv("results","tissue","gse254300_bulk","camera_CS_concave_minus_convex.tsv")
cam_msc <- read_tsv("results","tissue","gse110359_msc","camera_AIS_minus_healthy.tsv")
top <- cam_shared$pathway[order(cam_shared$FDR)][1:14]
mk <- function(x,axis){z<-x[x$pathway %in% top,c("pathway","Direction","FDR")];z$axis<-axis;z}
pmat <- rbind(mk(cam_shared,"Mean side"),mk(cam_ais,"AIS"),mk(cam_cs,"Congenital"))
pmat$score <- pmin(-log10(pmax(pmat$FDR,1e-15)),10)*ifelse(pmat$Direction=="Down",-1,1)
pmat$label <- factor(short_path(pmat$pathway),levels=rev(short_path(top)))
pmat$axis <- factor(pmat$axis,levels=c("Mean side","AIS","Congenital"))
p3a <- ggplot(pmat,aes(axis,label))+
  geom_tile(fill="#F5F6F7",colour="white",linewidth=.4)+
  geom_point(aes(size=-log10(FDR),fill=Direction),shape=21,colour="white",stroke=.3)+
  scale_fill_manual(values=c(Down=blue,Up=red))+scale_size_area(max_size=6,name="-log10 FDR")+
  guides(size=guide_legend(title.position="top",order=1,override.aes=list(fill=blue)),fill=guide_legend(title.position="top",order=2))+
  theme_nc()+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position="bottom",legend.box="horizontal",
    legend.text=element_text(size=6),legend.title=element_text(size=6.5),legend.key.width=unit(2.5,"mm"),legend.spacing.x=unit(1,"mm"))+
  labs(x=NULL,y=NULL,fill="Direction")
axis_map <- read_tsv("results","integration","orthogonal_axes","pathway_evidence_axis_map_camera.tsv")
axis_map$x <- -log10(pmax(axis_map$FDR_shared_side,1e-15))*ifelse(axis_map$Direction_shared_side=="Down",-1,1)
axis_map$y <- -log10(pmax(axis_map$FDR_msc,1e-15))*ifelse(axis_map$Direction_msc=="Down",-1,1)
axis_map$joint <- ifelse(axis_map$both_fdr05 & axis_map$direction_concordant,"Joint concordant",
                         ifelse(axis_map$both_fdr05,"Joint opposed","Not jointly significant"))
bcaa_rows <- which(grepl("BRANCHED_CHAIN_AMINO_ACID",axis_map$pathway) & axis_map$both_fdr05)
bcaa <- seq_len(nrow(axis_map)) == bcaa_rows[1]
p3b <- ggplot(axis_map,aes(x,y))+
  geom_hline(yintercept=c(log10(.05),-log10(.05)),colour=grid,linewidth=.35)+
  geom_vline(xintercept=c(log10(.05),-log10(.05)),colour=grid,linewidth=.35)+
  geom_point(aes(colour=joint),size=.9,alpha=.7)+
  geom_text_repel(data=axis_map[bcaa,],aes(label="Branched-chain amino-acid\ncatabolism"),size=2.3,family="sans",min.segment.length=0,seed=14016)+
  scale_colour_manual(values=c("Not jointly significant"="#C7CBD0","Joint opposed"=ochre,"Joint concordant"=teal),name=NULL)+
  theme_nc()+theme(legend.position="bottom")+
  labs(x="Paraspinal mean-side evidence\nsigned -log10 FDR",y="BM-MSC AIS-control evidence\nsigned -log10 FDR")
fig3 <- (p3a | p3b) + plot_layout(widths=c(1.25,1))+
  plot_annotation(tag_levels="a",theme=theme(plot.tag=element_text(family="sans",face="bold",size=10,colour=ink)))
save_figure(fig3,"Figure3_pathway_architecture",height=5.6)

# ---- Figure 4: composition sensitivity and single-cell localization ----
sens <- read_tsv("results","integration","composition_precision_sensitivity","camera_celltype_specificity_sensitivity.tsv")
top12 <- cam_shared$pathway[order(cam_shared$FDR)][1:12]
sd <- sens[sens$pathway %in% top12,]
sd$filter <- factor(sd$gene_filter,levels=c("all_retained_genes","exclude_specificity_fraction_ge_0p70","exclude_specificity_fraction_ge_0p50"),
                    labels=c("All genes","Exclude\n>=0.70","Exclude\n>=0.50"))
sd$label <- factor(short_path(sd$pathway),levels=rev(short_path(top12)))
p4a <- ggplot(sd,aes(filter,label))+
  geom_tile(fill="#F5F6F7",colour="white",linewidth=.4)+
  geom_point(aes(size=-log10(FDR),fill=Direction,alpha=FDR<.05),shape=21,colour="white",stroke=.3)+
  scale_fill_manual(values=c(Down=blue,Up=red))+scale_alpha_manual(values=c(`TRUE`=1,`FALSE`=.28),guide="none")+
  scale_size_area(max_size=6,name="-log10 FDR")+
  guides(size=guide_legend(title.position="top",order=1,override.aes=list(fill=blue)),fill=guide_legend(title.position="top",order=2))+
  theme_nc()+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position="bottom",legend.box="horizontal",
    legend.text=element_text(size=6),legend.title=element_text(size=6.5),legend.key.width=unit(2.5,"mm"),legend.spacing.x=unit(1,"mm"))+
  labs(x=NULL,y=NULL,fill="Direction")
sc <- read_tsv("results","integration","orthogonal_axes","celltype_pathways","donor_paired_celltype_pathway_effects.tsv")
focus <- top12[1:8]; sc <- sc[sc$pathway %in% focus,]
sc$label <- factor(short_path(sc$pathway),levels=rev(short_path(focus)))
cell_type_short <- c(
  Endothelial="Endo", FAP="FAP", Immune="Immune",
  Myogenic="Myo", Pericyte_SMC="Peri/SMC"
)
sc$cell_type_short <- unname(cell_type_short[sc$cell_type])
stopifnot(!anyNA(sc$cell_type_short))
sc$column <- paste(sc$donor,sc$cell_type_short,sep=" · ")
sc$column <- factor(sc$column,levels=unique(sc$column[order(sc$donor,sc$cell_type)]))
lim <- max(abs(sc$concave_minus_convex),na.rm=TRUE)
p4b <- ggplot(sc,aes(column,label,fill=concave_minus_convex))+
  geom_tile(colour="white",linewidth=.35)+
  scale_fill_gradient2(low=blue,mid="white",high=red,midpoint=0,limits=c(-lim,lim),name="Concave - convex\nmodule score")+
  guides(fill=guide_colourbar(barwidth=unit(27,"mm"),barheight=unit(2.4,"mm"),title.position="top",title.hjust=.5))+
  theme_nc()+theme(axis.line=element_blank(),axis.ticks=element_blank(),
                   axis.text.x=element_text(angle=48,hjust=1,size=6.2),
                   legend.position="bottom",legend.title=element_text(hjust=.5),
                   plot.margin=margin(4,12,4,5))+
  labs(x=NULL,y=NULL)
surv <- aggregate(FDR~gene_filter,data=sens[sens$pathway %in% cam_shared$pathway[cam_shared$FDR<.05],],FUN=function(z)sum(z<.05))
n_primary_pathways <- sum(cam_shared$FDR < 0.05, na.rm=TRUE)
surv$gene_filter <- factor(surv$gene_filter,levels=c("all_retained_genes","exclude_specificity_fraction_ge_0p70","exclude_specificity_fraction_ge_0p50"),labels=c("Original","Exclude >=0.70","Exclude >=0.50"))
p4c <- ggplot(surv,aes(gene_filter,FDR))+
  geom_col(width=.58,fill=teal)+geom_text(aes(label=paste0(FDR," / ",n_primary_pathways)),vjust=-.5,size=2.8,family="sans",fontface="bold")+
  scale_y_continuous(limits=c(0,28),breaks=c(0,10,20),expand=expansion(mult=c(0,.02)))+theme_nc()+
  labs(x=NULL,y=NULL)
fig4 <- (p4a | p4b) / p4c + plot_layout(widths=c(1.02,.98),heights=c(2.6,.72))+
  plot_annotation(tag_levels="a",theme=theme(plot.tag=element_text(family="sans",face="bold",size=10,colour=ink)))
save_figure(fig4,"Figure4_cellular_localization",height=7.25)

# ---- Figure 5: conservative tissue-to-plasma evidence funnel ----
counts <- read_tsv("results","integration","origin_annotation","origin_evidence_class_counts.tsv")
origin_labels <- c(
  "1_paraspinal_secretory_compatible_not_proven"="Secretory-compatible\nnot proven",
  "2_paraspinal_expressed_intracellular_or_damage_release"="Intracellular /\ndamage-compatible",
  "3_plasma_detected_paraspinal_low_or_absent"="Paraspinal low\nor absent",
  "4_broad_or_unresolved"="Broad /\nunresolved",
  "5_liver_or_hematologic_source_priority"="Liver or blood\nprioritized"
)
counts$short <- factor(
  unname(origin_labels[counts$cell_origin_evidence_class]),
  levels=rev(unname(origin_labels))
)
stopifnot(!anyNA(counts$short))
counts$colour <- class_cols
bridge <- read_tsv("results","integration","orthogonal_axes","shared_adaptation_plasma_bridge_candidates.tsv")
high_confidence_source <- read_tsv(
  "results","integration","origin_annotation",
  "high_confidence_paraspinal_secretory_compatible.tsv"
)
n_catalogue <- sum(counts$n_proteins)
n_source_compatible <- counts$n_proteins[
  counts$cell_origin_evidence_class == "1_paraspinal_secretory_compatible_not_proven"
]
p5a_dat <- data.frame(stage=factor(c("Plasma\ncatalogue","Source-\ncompatible","High-confidence\nsource-compatible","Tissue\nbridge"),
                                  levels=c("Plasma\ncatalogue","Source-\ncompatible","High-confidence\nsource-compatible","Tissue\nbridge")),
                      n=c(n_catalogue,n_source_compatible,nrow(high_confidence_source),nrow(bridge)),
                      x=c(1,2.15,3.25,4.2))
p5a_dat$y <- c(3,2,3,1)
branch_edges <- data.frame(x=c(1,2.15,2.15),xend=c(2.15,3.25,4.2),
                           y=c(3,2,2),yend=c(2,3,1))
p5a <- ggplot(p5a_dat,aes(x,n))+
  geom_segment(data=branch_edges,aes(x=x,xend=xend,y=y,yend=yend),inherit.aes=FALSE,
               linewidth=1.8,colour="#DCE1E5",lineend="round")+
  geom_point(aes(y=y),size=10,shape=21,fill=teal,colour="white",stroke=1)+
  geom_text(aes(y=y,label=n),colour="white",fontface="bold",family="sans",size=2.8)+
  geom_text(aes(y=y-.42,label=stage),family="sans",size=2.05,lineheight=.9)+
  coord_cartesian(xlim=c(.55,4.65),ylim=c(.25,3.5),clip="off")+theme_void()
p5b <- ggplot(counts,aes(n_proteins,short,fill=short))+
  geom_col(width=.62)+geom_text(aes(label=paste0(n_proteins,"  (",sprintf("%.1f",percent),"%)")),hjust=-.08,size=2.6,family="sans")+
  scale_fill_manual(values=setNames(class_cols,levels(counts$short)),guide="none")+
  scale_x_continuous(limits=c(0,365),breaks=c(0,100,200,300),expand=expansion(mult=c(0,.01)))+theme_nc()+
  labs(x="Proteins",y=NULL)
eq_bridge <- eq[match(bridge$gene_symbol,eq$gene_symbol),]
stopifnot(!anyNA(eq_bridge$gene_symbol),eq_bridge$gene_symbol == bridge$gene_symbol)
ev <- rbind(
  data.frame(gene=bridge$gene_symbol,evidence="Mean-side FDR < 0.05",value=bridge$shared_FDR < 0.05,detail=sprintf("%.2f",bridge$shared_logFC)),
  data.frame(gene=bridge$gene_symbol,evidence="Same direction in AIS and congenital",value=sign(bridge$ais_logFC)==sign(bridge$cs_logFC),detail="yes"),
  data.frame(gene=bridge$gene_symbol,evidence="Secretory support",value=bridge$plasma_classical_or_explicit_secretory_support,detail="yes"),
  data.frame(gene=bridge$gene_symbol,evidence="No prioritized systemic source",value=!bridge$plasma_competing_systemic_source,detail="yes"),
  data.frame(gene=bridge$gene_symbol,evidence="Equivalent at ±0.50",
             value=eq_bridge$strictly_conserved_margin_0p50,
             detail=ifelse(eq_bridge$strictly_conserved_margin_0p50,"yes","no"))
)
ev$evidence <- factor(ev$evidence,levels=unique(ev$evidence)); ev$gene <- factor(ev$gene,levels=rev(bridge$gene_symbol))
p5c <- ggplot(ev,aes(evidence,gene))+
  geom_tile(aes(fill=value),colour="white",linewidth=.7)+
  geom_text(aes(label=detail),size=2.55,family="sans",fontface="bold",colour=ink)+
  scale_fill_manual(values=c(`TRUE`="#CFE7DF",`FALSE`="#ECEDEF"),guide="none")+
  scale_x_discrete(labels=c("Mean-side effect\nFDR < 0.05","Concordant\nAIS and congenital","Secretory-pathway\nsupport","No prioritized\nsystemic source","Approximate QL\nequivalence ±0.50"))+
  theme_nc()+theme(axis.line=element_blank(),axis.ticks=element_blank(),
                   axis.text.x=element_text(angle=0,hjust=.5,size=6.8,lineheight=1.1),
                   plot.margin=margin(4,8,6,8))+
  labs(x=NULL,y=NULL)
fig5 <- (p5a | p5b) / p5c + plot_layout(heights=c(1.2,.9),widths=c(1.02,.98))+
  plot_annotation(tag_levels="a",theme=theme(plot.tag=element_text(family="sans",face="bold",size=10,colour=ink)))
save_figure(fig5,"Figure5_plasma_bridge",height=6.75)

p5a_export <- p5a_dat
p5a_export$stage <- gsub("\n"," ",as.character(p5a_export$stage),fixed=TRUE)
write.table(p5a_export,file.path(outdir,"Figure5_funnel.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c(
  paste("generated",format(Sys.time(),tz="UTC",usetz=TRUE),sep="\t"),
  paste("R",R.version.string,sep="\t"),
  paste("ggplot2",as.character(packageVersion("ggplot2")),sep="\t"),
  paste("patchwork",as.character(packageVersion("patchwork")),sep="\t"),
  "visual_standard\t183-mm double-column composites; sans 7.5 pt base; colourblind-conscious direction palette",
  "claim_rule\tgraphics summarize analysis-native tables and do not introduce new inferential tests"
),file.path(outdir,"figure_manifest.tsv"))

cat("Built five Nature Communications-style main figures in",outdir,"\n")

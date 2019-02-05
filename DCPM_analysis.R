# load TCGA datasets and map sample IDs
CNVtable.tcga <- readRDS("TCGA_OV_CBS.rds")
meth.tcga <- readRDS("TCGA_OV_meth.rds")
mapped <- sort(intersect(colnames(CNVtable.tcga),colnames(meth.tcga)))
CNVtable.tcga.mapped <- CNVtable.tcga[,mapped]
meth.tcga.mapped <- meth.tcga[,mapped]
cpg.annot.tcga <- readRDS("Illumina27k_annotation.rds")
cpg.genes.tcga <- as.character(cpg.annot.tcga[rownames(meth.tcga.mapped),1])
mapped.cpgs.tcga <- which(cpg.genes.tcga %in% rownames(CNVtable.tcga.mapped))
CNVtable.tcga.mapped <- CNVtable.tcga.mapped[cpg.genes.tcga[mapped.cpgs.tcga],]
meth.tcga.mapped <- meth.tcga.mapped[mapped.cpgs.tcga,]
gx.tcga <- readRDS("TCGA_OV_GX.rds")
mapped2.tcga <- intersect(colnames(meth.tcga),intersect(colnames(CNVtable.tcga),colnames(gx.tcga)))
CNV.order.tcga2 <- sapply(mapped2.tcga,function(x)which(colnames(CNVtable.tcga)==x))
meth.order.tcga2 <- sapply(mapped2.tcga,function(x)which(colnames(meth.tcga)==x))
gx.order.tcga2 <- sapply(mapped2.tcga,function(x)which(colnames(gx.tcga)==x))
mapped.cpgs.tcga2 <- which(cpg.genes.tcga %in% intersect(rownames(CNVtable.tcga.mapped),rownames(gx.tcga)))

# load ICGC datasets (already mapped to sample IDs)
icgc.IDs <- read.table("sample.OV-AU.tsv",sep="\t",head=T,row.names=1)
CNVtable.icgc <- readRDS("ICGC_OV-AU_CNVbyCpG.rds")
meth.icgc <- readRDS("ICGC_OV-AU_Meth.rds")
cpg.annot2 <- readRDS("full450k.rds")
gx.icgc <- readRDS("ICGC_OV-AU_GX.rds")

# find promoter probes with min 10% copy-gain
allpromoter.450k.probes = rownames(cpg.annot2)[which(cpg.annot2$UCSC_RefGene_Group=="TSS200")]
icgc.amprate <- apply(CNVtable.icgc,1,function(x)sum(x>log(2.4/2,base=2),na.rm=T)/sum(!is.na(x)))
icgc.promoter.minamp <- intersect(rownames(CNVtable.icgc)[which(icgc.amprate>0.1)],allpromoter.450k.probes)

cpg.annot.tcga2 <- readRDS("Illumina27k_annotation2.rds")
allpromoter.27k.probes = rownames(cpg.annot.tcga2)[which(cpg.annot.tcga2$Distance_to_TSS<200)]
tcga.amprate <- apply(CNVtable.tcga.mapped,1,function(x)sum(x>2.5,na.rm=TRUE)/sum(!is.na(x)))
tcga.minamp.genes <- which(tcga.amprate>0.1)
tcga.promoter.minamp <- allpromoter.27k.probes[which(cpg.annot.tcga2[allpromoter.27k.probes,"Symbol"] %in% rownames(CNVtable.tcga.mapped)[tcga.minamp.genes] & allpromoter.27k.probes %in% rownames(meth.tcga.mapped))]

# compute meth-CNV correlations

icgc.promoter.compMeth <- rep(NA,length(icgc.promoter.minamp))
icgc.promoter.compMeth.p <- rep(NA,length(icgc.promoter.minamp))
for(i in 1:length(icgc.promoter.minamp)){
  icgc.promoter.compMeth.p[i] <- cor.test(x=as.numeric(meth.icgc[icgc.promoter.minamp[i],]),y=as.numeric(CNVtable.icgc[icgc.promoter.minamp[i],]),method="spearman")$p.value
  icgc.promoter.compMeth[i] <- cor(x=as.numeric(meth.icgc[icgc.promoter.minamp[i],]),y=as.numeric(CNVtable.icgc[icgc.promoter.minamp[i],]),method="spearman",use="complete")
}

tcga.promoter.compMeth <- rep(NA,length(tcga.promoter.minamp))
tcga.promoter.compMeth.p <- rep(NA,length(tcga.promoter.minamp))
for(i in 1:length(tcga.promoter.minamp)){
  tcga.promoter.compMeth.p[i] <- cor.test(x=as.numeric(meth.tcga.mapped[tcga.promoter.minamp[i],]),y=as.numeric(CNVtable.tcga.mapped[as.character(cpg.annot.tcga2[tcga.promoter.minamp[i],"Symbol"]),]),method="spearman")$p.value
  tcga.promoter.compMeth[i] <- cor(x=as.numeric(meth.tcga.mapped[tcga.promoter.minamp[i],]),y=as.numeric(CNVtable.tcga.mapped[as.character(cpg.annot.tcga2[tcga.promoter.minamp[i],"Symbol"]),]),method="spearman",use="complete")
}

# Fig 1A
plot(x=exp(seq(from=log(0.1),to=log(0.001),by=(-0.2))),y=log(sapply(seq(from=log(0.1),to=log(0.001),by=(-0.2)),function(x)sum(icgc.promoter.compMeth.p<exp(x) & icgc.promoter.compMeth<0,na.rm=T)/((sum(!is.na(icgc.promoter.compMeth.p))*exp(x))/2)),base=10),xlab="Spearman correlation p-value",ylab="log10 fold-enrichment",type="l",col="blue",lwd=2,ylim=c(0,3),main="Enrichment of CNV-meth associated loci: ICGC")
points(x=exp(seq(from=log(0.1),to=log(0.001),by=(-0.2))),y=log(sapply(seq(from=log(0.1),to=log(0.001),by=(-0.2)),function(x)sum(icgc.promoter.compMeth.p<exp(x) & icgc.promoter.compMeth>0,na.rm=T)/((sum(!is.na(icgc.promoter.compMeth.p))*exp(x))/2)),base=10),type="l",col="green",lwd=2)
abline(h=log(2,base=10),lty=2,lwd=2)> legend("topright",legend=c("Meth~1/CNV","Meth~CNV","2-fold enrichment"),col=c("blue","green","black"),lwd=2,lty=c(1,1,2))

# Fig 1B
plot(x=exp(seq(from=log(0.1),to=log(0.001),by=(-0.2))),y=log(sapply(seq(from=log(0.1),to=log(0.001),by=(-0.2)),function(x)sum(tcga.promoter.compMeth.p<exp(x) & tcga.promoter.compMeth<0,na.rm=T)/((sum(!is.na(tcga.promoter.compMeth.p))*exp(x))/2)),base=10),xlab="Spearman correlation p-value",ylab="log10 fold-enrichment",type="l",col="blue",lwd=2,ylim=c(0,3),main="Enrichment of CNV-meth associated loci: TCGA")
points(x=exp(seq(from=log(0.1),to=log(0.001),by=(-0.2))),y=log(sapply(seq(from=log(0.1),to=log(0.001),by=(-0.2)),function(x)sum(tcga.promoter.compMeth.p<exp(x) & tcga.promoter.compMeth>0,na.rm=T)/((sum(!is.na(tcga.promoter.compMeth.p))*exp(x))/2)),base=10),type="l",col="green",lwd=2)
abline(h=log(2,base=10),lty=2,lwd=2)> legend("topright",legend=c("Meth~1/CNV","Meth~CNV","2-fold enrichment"),col=c("blue","green","black"),lwd=2,lty=c(1,1,2))

# Fig 1C
plot(x=sapply(intersect(tcga.promoter.minamp,icgc.promoter.minamp),function(x)tcga.promoter.compMeth[which(tcga.promoter.minamp==x)]),y=sapply(intersect(tcga.promoter.minamp,icgc.promoter.minamp),function(x)icgc.promoter.compMeth[which(icgc.promoter.minamp==x)]),xlab="TCGA CNV-Meth correlation (rho)",ylab="ICGC CNV-Meth correlation (rho)")
abline(v=0,lty=2)
abline(h=0,lty=2)

# compute meth-cellularity correlations

icgc.cellularity <- icgc.IDs[colnames(meth.icgc),"percentage_cellularity"]
icgc.promoter.cellularity <- rep(NA,length(icgc.promoter.minamp))
icgc.promoter.cellularity.p <- rep(NA,length(icgc.promoter.minamp))
for(i in 1:length(icgc.promoter.minamp)){
  icgc.promoter.cellularity.p[i] <- cor.test(x=as.numeric(meth.icgc[icgc.promoter.minamp[i],]),y=as.numeric(icgc.cellularity),method="spearman")$p.value
  icgc.promoter.cellularity[i] <- cor(x=as.numeric(meth.icgc[icgc.promoter.minamp[i],]),y=as.numeric(icgc.cellularity),method="spearman",use="complete")
}

tcga.slide.data <- read.table("clinical_slide_all_ov.txt",sep="\t",head=T)
tcga.slide.data <- tcga.slide.data[which(substr(as.character(tcga.slide.data[,1]),start=14,stop=15)=="01"),]
tcga.slide.barcodes <- substr(as.character(tcga.slide.data[,1]),start=1,stop=12)
tcga.cellularity <- sapply(colnames(meth.tcga.mapped),function(x)median(as.numeric(as.character(tcga.slide.data[which(tcga.slide.barcodes==gsub(x,pattern=".",replace="-",fixed=T)),"percent_tumor_cells"]))))

tcga.promoter.cellularity <- rep(NA,length(tcga.promoter.minamp))
tcga.promoter.cellularity.p <- rep(NA,length(tcga.promoter.minamp))
for(i in 1:length(tcga.promoter.minamp)){
  tcga.promoter.cellularity.p[i] <- cor.test(x=as.numeric(meth.tcga.mapped[tcga.promoter.minamp[i],]),y=tcga.cellularity,method="spearman")$p.value
  tcga.promoter.cellularity[i] <- cor(x=as.numeric(meth.tcga.mapped[tcga.promoter.minamp[i],]),y=tcga.cellularity,method="spearman",use="complete")
}

# list genes with DCPM in both cohorts

icgc.dcpm.genes <- as.character(cpg.annot2[icgc.promoter.minamp[icgc.promoter.compMeth>0 & icgc.promoter.compMeth.p<0.1 & icgc.promoter.cellularity.p>0.1],"UCSC_RefGene_Name"])
tcga.dcpm.genes <- as.character(cpg.annot.tcga[tcga.promoter.minamp[tcga.promoter.compMeth>0 & tcga.promoter.compMeth.p<0.1 & tcga.promoter.cellularity.p>0.1],"GENESYMBOL"])
joint.dcpm.genes <- intersect(icgc.dcpm.genes,tcga.dcpm.genes)

# estimate statistical significance of overlap
length(joint.dcpm.genes)
#[1] 24
length(intersect(tcga.dcpm.genes,intersect(tcga.mappedgenes,mappedgenes.icgc)))
#[1] 187
length(setdiff(intersect(tcga.mappedgenes,mappedgenes.icgc),tcga.dcpm.genes))
#[1] 3035
length(intersect(icgc.dcpm.genes,intersect(tcga.mappedgenes,mappedgenes.icgc)))
#[1] 58
1-phyper(23,187,3035,58)
#[1] 1.110223e-15

# compute CNV-gx correlations

gx.mappedloci.tss1500 <- rownames(meth.icgc)[which(as.character(cpg.annot2[rownames(meth.icgc),"UCSC_RefGene_Name"]) %in% rownames(gx.icgc) & as.character(cpg.annot2[rownames(meth.icgc),"UCSC_RefGene_Group"]) %in% c("TSS200","TSS1500") & !duplicated(as.character(cpg.annot2[rownames(meth.icgc),"UCSC_RefGene_Name"])))]
mappedgenes.icgc = as.character(cpg.annot2[gx.mappedloci.tss1500,"UCSC_RefGene_Name"])
icgc.cnvtogx.rho = rep(NA,length(mappedgenes.icgc))
icgc.cnvtogx.p = rep(NA,length(mappedgenes.icgc))
for(i in 1:length(mappedgenes.icgc)){
  if(sum(!is.na(as.numeric(CNVtable.icgc[gx.mappedloci.tss1500[i],])))>10){
    icgc.cnvtogx.rho[i] <- cor(x=as.numeric(CNVtable.icgc[gx.mappedloci.tss1500[i],]),y=gx.icgc[mappedgenes.icgc[i],],method="spearman")
    icgc.cnvtogx.p[i] <- cor.test(x=as.numeric(CNVtable.icgc[gx.mappedloci.tss1500[i],]),y=gx.icgc[mappedgenes.icgc[i],],method="spearman")$p.value
  }
}

tcga.mappedgenes = intersect(intersect(rownames(gx.tcga),rownames(CNVtable.tcga)),as.character(cpg.annot.tcga2[rownames(meth.tcga.mapped),"Symbol"]))
tcga.cnvtogx.rho = rep(NA,length(tcga.mappedgenes))
tcga.cnvtogx.p = rep(NA,length(tcga.mappedgenes))
for(i in 1:length(tcga.mappedgenes)){
  tcga.cnvtogx.rho[i] <- cor(x=CNVtable.tcga[tcga.mappedgenes[i],CNV.order.tcga2],y=gx.tcga[tcga.mappedgenes[i],gx.order.tcga2],method="spearman")
  tcga.cnvtogx.p[i] <- cor.test(x=CNVtable.tcga[tcga.mappedgenes[i],CNV.order.tcga2],y=gx.tcga[tcga.mappedgenes[i],gx.order.tcga2],method="spearman")$p.value
}



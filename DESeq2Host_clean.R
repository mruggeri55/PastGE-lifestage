setwd('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/')

library(DESeq2)
library(gplots)
library(pheatmap)
library(vsn)
library(RColorBrewer)
library(ggplot2)

cts=read.table("computational/CG_AllCountsHost_new.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample
head(cts) 
length(cts[,1])  #34805 isogroups

names(cts)<-gsub("X","",names(cts))
names(cts)

######## read in metadata

colData=read.csv('downstream_analyses/Mar2020/DESeq/input_files/CG_colData.csv',row.names = 1)
colData$genotype=as.factor(colData$genotype)
head(colData)

######### Quality control

#cds=estimateDispersions(cts,method="blind")
#vsdBlind=varianceStabilizingTransformation(cds)

#v="/Users/maria/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/downstream_analyses/Mar2020/DESeq/arrayQualityMetrics/"

#arrayQualityMetrics(vsdBlind,outdir=v,intgroup=c("genotype"),force=TRUE) #check .html output file in new folder

#47HL has absolutely got to go. It's ticking all 3 outlier boxes. Remove this sample and run again 

counts.nobad1=cts[,-c(36)]  
#conditions.nobad1=colData[-c(36),]

# real=newCountDataSet(counts.nobad1,conditions.nobad1) 
# real=estimateSizeFactors(real)
# 
# cds=estimateDispersions(real,method="blind")
# vsdBlind=varianceStabilizingTransformation(cds)
# 
# arrayQualityMetrics(vsdBlind,outdir=v,intgroup=c("genotype"),force=TRUE) #check .html output file in new folder

#other potential outliers aren't extreme, will retain remaining samples

#go ahead and remove low expression genes - those with count less than 2 in more than 90% of samples
counts.nobad1$low = apply(counts.nobad1[,1:47],1,function(x){sum(x<=2)})  #making new column counting number of samples with counts <=10 within each isogroup (host) <=5 within each dataset
counts.nobad1<-counts.nobad1[-which(counts.nobad1$low>43),1:47] #42.3 is 90% of 47 samples
nrow(counts.nobad1) #20185
head(counts.nobad1)

counts=counts.nobad1
#colData=conditions.nobad1

####################### starting DESeq2 analysis -- explore pca of all genes
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ 1)
dds # 20185 isogroups

# differential expression analysis - DESeq will calculate the log2fold and Wald test p-value for last variable in design formula
vsd_null<-vst(dds,blind=TRUE)
plotPCA(vsd_null,intgroup=c("treatment",'stage'))
vsd_null_adults<-vsd_null[,vsd_null$stage=="Adult"]
plotPCA(vsd_null_adults,intgroup='treatment')
plotPCA(vsd_null_adults,intgroup='origin')
plotPCA(vsd_null_adults,intgroup='genotype')
vsd_null_recruits<-vsd_null[,vsd_null$stage=="Recruit"]
plotPCA(vsd_null_recruits,intgroup='treatment')
plotPCA(vsd_null_recruits,intgroup='origin')


pca_dat=plotPCA(vsd_null,intgroup=c("origin",'stage','treatment'),returnData=T)
percentVar <- round(100 * attr(pca_dat, "percentVar"))
pca_dat$origin_treatment=paste(pca_dat$origin,pca_dat$treatment,sep=' ')
ggplot(pca_dat, aes(x=PC1, y=PC2,colour=origin_treatment,fill=origin_treatment, shape=stage)) +
  geom_point(size = 5,stroke=1) +
  scale_shape_manual(name= 'life stage',values=c(21,24))+
  scale_fill_manual(name='origin treatment',
                    labels=c('inshore control','inshore heat','offshore control','offshore heat'),
                    values=c("orange",alpha("orange",0.3),"aquamarine3",alpha("aquamarine3",0.3)),
                    breaks=c('in Control','in Heat','off Control','off Heat'),
                    guide='legend') +
  scale_colour_manual(name='origin treatment',
                      labels=c('inshore control','inshore heat','offshore control','offshore heat'),
                      values=c("orange",'orange','aquamarine3', "aquamarine3"),
                      breaks=c('in Control','in Heat','off Control','off Heat')) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  guides(fill= guide_legend(override.aes = list(alpha=c(0.3,1,0.3,1))))+
  theme_light() + theme(text=element_text(size=20))

# recruit pca
pca_dat=plotPCA(vsd_null_recruits,intgroup=c("origin",'treatment'),returnData=T)
percentVar <- round(100 * attr(pca_dat, "percentVar"))
pca_dat$origin_treatment=as.factor(paste(pca_dat$origin,pca_dat$treatment,sep=' '))
recruits=ggplot(pca_dat, aes(x=PC1, y=PC2,colour=origin_treatment,fill=origin_treatment)) +
  geom_point(size = 3,stroke=1,shape=24) +
  scale_fill_manual(name='origin treatment',
                    labels=c('inshore control','inshore heat','offshore control','offshore heat'),
                    values=c("orange",alpha("orange",0.3),"aquamarine3",alpha("aquamarine3",0.3)),
                    breaks=c('in Control','in Heat','off Control','off Heat'),
                    guide='legend') +
  scale_colour_manual(name='origin treatment',
                      labels=c('inshore control','inshore heat','offshore control','offshore heat'),
                      values=c("orange",'orange','aquamarine3', "aquamarine3"),
                      breaks=c('in Control','in Heat','off Control','off Heat')) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  guides(fill= guide_legend(override.aes = list(alpha=c(1,0.3,1,0.3))))+
  theme_light() + theme(text=element_text(size=16))+ggtitle('Recruits')

#save legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(recruits)

recruits <- recruits + theme(legend.position = 'none')

# adult pca
pca_dat=plotPCA(vsd_null_adults,intgroup=c("origin",'treatment'),returnData=T)
percentVar <- round(100 * attr(pca_dat, "percentVar"))
pca_dat$origin_treatment=as.factor(paste(pca_dat$origin,pca_dat$treatment,sep=' '))
pca_dat$PC2=-1*pca_dat$PC2
adults=ggplot(pca_dat, aes(x=PC1, y=PC2,colour=origin_treatment,fill=origin_treatment)) +
  geom_point(size = 3,stroke=1,shape=21) +
  scale_fill_manual(name='origin treatment',
                    labels=c('inshore control','inshore heat','offshore control','offshore heat'),
                    values=c("orange",alpha("orange",0.3),"aquamarine3",alpha("aquamarine3",0.3)),
                    breaks=c('in Control','in Heat','off Control','off Heat'),
                    guide='legend') +
  scale_colour_manual(name='origin treatment',
                      labels=c('inshore control','inshore heat','offshore control','offshore heat'),
                      values=c("orange",'orange','aquamarine3', "aquamarine3"),
                      breaks=c('in Control','in Heat','off Control','off Heat')) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  guides(fill= guide_legend(override.aes = list(alpha=c(0.3,1,0.3,1))))+
  theme_light() + theme(text=element_text(size=16),legend.position = 'none')+ggtitle('Adults')

#library(gridExtra)
#grid.arrange(adults,larvae,recruits,legend,nrow=1,widths=c(2.5, 2.5, 2.5, 1.2))

# take a look at PC3
library(ggfortify)
pca=prcomp(t(assay(vsd_null)),scale. = TRUE)
pca=as.data.frame(pca$x)
colData$stage_treatment=paste(colData$stage,colData$treatment)
ggplot(pca,aes(x=PC2, y=PC3, colour=colData$stage_treatment,fill=colData$stage_treatment,shape=colData$stage)) +
  geom_point(size=3)+
  scale_shape_manual(values = c(21, 24)) + scale_fill_manual(values=c('dark blue','light blue','seagreen','darkseagreen1'))+
  scale_color_manual(values=c('dark blue','dark blue','seagreen','seagreen'))+labs(fill='treatment group')
# does not separate by treatment along PC3

############# ~ stage + origin + treatment ###########
#design(dds) <- ~ genotype + stage + origin + treatment # including family makes origin have much less DEGs, about same number for treatment and stage
#design(dds) <- ~ stage + genotype + treatment
design(dds) <- ~ stage + origin + treatment

dds <- DESeq(dds)
plotDispEsts(dds)

# dds <- DESeq(dds, test='LRT', reduced = ~ stage + origin + treatment)
# summary(results(dds)) #this says genotype only accounts for 144 genes
# including genotype (family) greatly reduced genes DE by origin

# dds <- DESeq(dds, test='LRT', reduced = ~ stage + treatment)
# summary(results(dds)) #this says genotype accounts for 1,729 genes
# fam_genes=rownames(subset(results(dds),padj < 0.1))
# origin_genes=rownames(subset(res_origin,padj < 0.1))
# candidates=list('fam'=fam_genes,'origin'=origin_genes)
# venn(candidates) # 984 genes overlap between the two, 318 unique to origin, 745 unique to fam

resultsNames(dds)
res <- results(dds) # this is the treatment effect bc last term in model
res_stage <- results(dds, name='stage_Recruit_vs_Adult')
res_origin <- results(dds, name='origin_off_vs_in' ) #not including fam: 1302, including fam: 135
#res_fam <- results(dds, name='genotype') #146 genes DE by family

metadata(res)$filterTheta #0.35

plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)

p.adj.cutoff=0.1
# add T/F column if gene is sig
res$threshold <- as.logical(res$padj < p.adj.cutoff)
res_stage$threshold <- as.logical(res_stage$padj < p.adj.cutoff)
res_origin$threshold <- as.logical(res_origin$padj < p.adj.cutoff)
# Sort the results tables
res_sorted <- res[order(res$padj), ]
res_stage_sorted <- res_stage[order(res_stage$padj), ]
res_origin_sorted <- res_origin[order(res_origin$padj), ]
# Get significant genes
sigTrmt <- row.names(res_sorted)[which(res_sorted$threshold)]
sigStage <- row.names(res_stage_sorted)[which(res_stage_sorted$threshold)]
sigOrigin <- row.names(res_origin_sorted)[which(res_origin_sorted$threshold)]
sigTrmtOnly <- sigTrmt[!sigTrmt %in% c(sigStage,sigOrigin)]
sigStageOnly <- sigStage[!sigStage %in% c(sigTrmt,sigOrigin)]
sigOriginOnly <- sigOrigin[!sigOrigin %in% c(sigTrmt,sigStage)]
sigTrmtStage <- sigTrmt[sigTrmt %in% sigStage]
sigAll <- sigTrmtStage[sigTrmtStage %in% sigOrigin]

candidates=list('Treatment'=sigTrmt,'Stage'=sigStage,'Origin'=sigOrigin)
venn(candidates)

# make pretty venn diagram
library(VennDiagram)
venn.diagram(
  x=candidates,
  category.names = c('Treatment','Stage','Origin'),
  filename='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/plots/HOST_Venn_stage_origin_trmt.png',
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 2,
  cat.col = c("#440154ff", '#21908dff', 'black'),
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 180),
  cat.dist = c(0.08, 0.08, 0.05)
)

# are these trmt genes the same as DESeq1?
# read in DESeq output and see if sig trmt genes are the same
DESeq1=read.csv('downstream_analyses/Mar2020/DESeq/input_files/VSDandPVALS_SDtheta02.csv',row.names = 1)
DESeq1_trmt_genes=DESeq1[DESeq1$adjp.t<=0.1 & !is.na(DESeq1$adjp.t),]
overlap = DESeq1_trmt_genes[rownames(DESeq1_trmt_genes) %in% sigTrmt,] #559 genes overlap between DESeq 1 and 2 results
# aka almost all detected in DESeq2 are detected in DESeq but DESeq is detecting a lot more
overlap2 = sigTrmt[sigTrmt %in% rownames(DESeq1_trmt_genes)] #559 / 605 DESeq2 trmt DEGs

# transform counts for plotting
vsd=vst(dds,blind=FALSE)
norm_df=assay(vsd)
means=apply(norm_df,1,mean) # means of rows
norm_df_scaled=norm_df-means #rescale expression data so it's up and down
df=as.data.frame(colData(dds)[,c("treatment","stage","origin")])
display.brewer.pal(n = 6, name = 'Set2')
brewer.pal(n = 6, name = "Set2")
ann_colors=list(
  treatment=c(Control='#66C2A5', Heat='#FC8D62'),
  origin=c(off='#A6D854','in'='#FFD92F'),
  stage=c(Adult='#8DA0CB',Recruit='#E78AC3')
)
col=color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.7)(100)
pheatmap(norm_df_scaled[sigTrmt,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment then stage for all sig trmt DEGs
col=color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.4)(100)
pheatmap(norm_df_scaled[sigTrmtOnly,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment, not really stage

col=color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.7)(100)
pheatmap(norm_df_scaled[sigTrmtStage,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors)

col=color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.5)(100)
pheatmap(norm_df_scaled[sigOrigin,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors)
col=color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.2)(100)
pheatmap(norm_df_scaled[sigOriginOnly,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors)

col=color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.4)(100)
pheatmap(norm_df_scaled[sigStage,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors)
col=color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.3)(100)
pheatmap(norm_df_scaled[sigStageOnly,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors)

col=color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.3)(100)
pheatmap(norm_df_scaled[sigAll,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors)

out=as.data.frame(cbind(norm_df,res$padj,res_stage$padj,res_origin$padj,res$log2FoldChange,res_stage$log2FoldChange,res_origin$log2FoldChange))
names(out)[48:53]=c('trmt_padj','stage_padj','origin_padj','trmt_LFC','stage_LFC','origin_LFC')

write.csv(out,file='downstream_analyses/Mar2020/DESeq/DESeq2_vsd_pvals_30June2021.csv')

##### generating GO output
#treatment
res$direction=ifelse(res$log2FoldChange>0,1,-1) #red=higher expression in heat
res$logP<-(-log((res$padj+0.0000000001),10))*(res$direction)

#origin
res_origin$direction=ifelse(res_origin$log2FoldChange>0,1,-1) #red=higher expression in offshore
res_origin$logP<-(-log((res_origin$padj+0.0000000001),10))*(res_origin$direction)

#stage
res_stage$direction=ifelse(res_stage$log2FoldChange>0,1,-1) # red = higher expression in recruit
res_stage$logP<-(-log((res_stage$padj+0.0000000001),10))*(res_stage$direction)

trmt_out<-as.data.frame(cbind("gene"=row.names(res),"logP"=res$logP))
origin_out<-as.data.frame(cbind("gene"=row.names(res_origin),"logP"=res_origin$logP))
stage_out<-as.data.frame(cbind("gene"=row.names(res_stage),"logP"=res_stage$logP))

write.csv(trmt_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankHost_ARtrmt.csv",quote=F,row.names=F)
write.csv(origin_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankHost_ARorigin.csv",quote=F,row.names=F)
write.csv(stage_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankHost_ARstage.csv",quote=F,row.names=F)







################################################################
########### for the interaction of treatment and life stage ###########
###### control for origin, then look at stage-trmt responses
###### first try with interaction in formula -- go to ?results example #2 for how to extract results
design(dds) <- ~ origin + stage + treatment + stage:treatment
dds <- DESeq(dds)
resultsNames(dds)

res_trmt_Adult <- results(dds, contrast=c('treatment','Heat','Control')) #this is treatment effect for Adults, exactly the same as res_Adult
res_trmt_Recruit <- results(dds, contrast=list(c('treatment_Heat_vs_Control','stageRecruit.treatmentHeat'))) # treatment effect for recruits, same as res_Recruit
res_int <- results(dds, name='stageRecruit.treatmentHeat') # this is the interaction term -- 3 up, 6 down genes
int_genes=rownames(subset(res_int,padj<0.1))
iso2gene=read.table('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/past_blast_annot_July2021/past_iso2gene.tab',sep = '\t', colClasses = "character",row.names = 1)
int_gnames=iso2gene[rownames(iso2gene) %in% int_genes,]
int_gnames
# [1] "Mitogen-activated protein kinase 6 OS=Mus musculus OX=10090 GN=Mapk6 PE=1 SV=3 E(blastx)=5e-94"
# [2] "Endothelin-converting enzyme 2 OS=Homo sapiens OX=9606 GN=ECE2 PE=1 SV=1 E(blastx)=1e-143"

#pull out LFCs for int genes in adults and recruits
adult_int_genes=subset(res_trmt_Adult,row.names(res_trmt_Adult) %in% int_genes,select=log2FoldChange)
recruit_int_genes=subset(res_trmt_Recruit,row.names(res_trmt_Recruit) %in% int_genes,select=log2FoldChange)
rownames(adult_int_genes)==rownames(recruit_int_genes)
int_genes_LFCs=merge(as.data.frame(adult_int_genes),as.data.frame(recruit_int_genes),by='row.names')
names(int_genes_LFCs)=c('isogroup','adult_LFC','recruit_LFC')
rownames(int_genes_LFCs)=int_genes_LFCs$isogroup
int_genes_LFCs=merge(int_genes_LFCs,iso2gene,by='row.names',all.x=T)
colnames(int_genes_LFCs)[5]='gene'
int_genes_LFCs=int_genes_LFCs[,2:5]
write.csv(int_genes_LFCs,file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/interesting_genes/TreatbyStage_interaction_HostAR.csv')

write.csv(as.data.frame(res_int),file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_TreatByStage_interaction.csv')




######## combine the factors of interest into a single factor with all combinations of the original factors
dds$group <- factor(paste0(dds$stage, dds$treatment))
# change the design to include just this factor, e.g. ~ group
design(dds) <- ~ origin + group
dds <- DESeq(dds)
res_group=results(dds)
summary(res_group)
resultsNames(dds)

res_Adult <- results(dds, contrast=c("group", "AdultHeat", "AdultControl")) # 102 up, 92 down
res_Recruit <- results(dds, contrast=c("group", "RecruitHeat", "RecruitControl")) # 120 up, 191 down
res_stage_control <- results(dds, contrast=c('group', 'AdultControl', 'RecruitControl')) #1881 up, 2017 down
res_stage_heat <- results(dds, contrast=c('group','AdultHeat','RecruitHeat')) #1901 up, 2233 down
res_origin <- results(dds,contrast=c('origin','in','off')) #775 up, 533 down

Larv_res=read.csv('downstream_analyses/Mar2020/DESeq/DESeq2_June2021/Larv_DESeq2_VsdPvalsLFC_19July2021.csv',row.names = 1)
Larv_trmtGenes=rownames(subset(Larv_res,trmt_padj < 0.1))
Adult_trmtGenes=rownames(subset(res_Adult,padj < 0.1))
Recruit_trmtGenes=rownames(subset(res_Recruit,padj < 0.1))

candidates=list('Adult'=Adult_trmtGenes,'Recruit'=Recruit_trmtGenes,'Larvae'=Larv_trmtGenes)
venn(candidates)

# pull out overlapping genes
AR=Adult_trmtGenes[Adult_trmtGenes %in% Recruit_trmtGenes] # 44 overlapping genes between adults and recruits
AL=Adult_trmtGenes[Adult_trmtGenes %in% Larv_trmtGenes] # 23 genes overlapping between adults and larvae
RL=Recruit_trmtGenes[Recruit_trmtGenes %in% Larv_trmtGenes] # 20 genes overlapping between recruits and larvae
ARL=Adult_trmtGenes[Adult_trmtGenes %in% Recruit_trmtGenes & Adult_trmtGenes %in% Larv_trmtGenes] # 9 genes responding in all life stages

# what are these genes?
iso2gene=read.table('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/past_blast_annot_July2021/past_iso2gene.tab',sep = '\t', colClasses = "character")

AR=iso2gene[iso2gene$V1 %in% AR & !iso2gene$V1 %in% ARL,] # 10 annotated
AL=iso2gene[iso2gene$V1 %in% AL & !iso2gene$V1 %in% ARL,] # 4 annotated
RL=iso2gene[iso2gene$V1 %in% RL & !iso2gene$V1 %in% ARL,] # 5 annotated
ARL=iso2gene[iso2gene$V1 %in% ARL,] # 5 annotated
all=rbind(AR,AL,RL,ARL)

all$stages=c(rep('A&R',times=10),rep('A&L',times=4),rep('R&L',times=5),rep('A&R&L',times=5))

adult_LFC_sub=subset(res_Adult,subset=rownames(res_Adult) %in% all$V1,select=log2FoldChange)
recruit_LFC_sub=subset(res_Recruit,subset=rownames(res_Recruit) %in% all$V1,select=log2FoldChange)
larv_LFC_sub=subset(Larv_res,subset=rownames(Larv_res) %in% all$V1,select=trmt_LFC)

all_order=all[match(rownames(adult_LFC_sub),all$V1),]
merge_df=as.data.frame(cbind(all_order,adult_LFC_sub,recruit_LFC_sub,larv_LFC_sub))
merge_df$gene=as.character(lapply(strsplit(merge_df$V2, split='OS='),'[',1))
merge_df$GN=as.character(lapply(strsplit(merge_df$V2, split='GN='),'[',2))
merge_df$GN=as.character(lapply(strsplit(merge_df$GN, split='PE='),'[',1))
merge_df$V2=NULL
names(merge_df)=c('isogroup','stages','adult (LFC)', 'recruit (LFC)', 'larvae (LFC)','gene','GN')
write.csv(x=merge_df,file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GeneNames_trmt_by_stage.csv')

# make pretty venn diagram
library(VennDiagram)
venn.diagram(
  x=candidates,
  category.names = c('Adult','Recruit','Larvae'),
  filename='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/TrmtResponse_byLifeStageVenn.png',
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 2,
  cat.col = c("#440154ff", '#21908dff', 'black'),
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 180),
  cat.dist = c(0.08, 0.08, 0.05)
)

venn.diagram(
  x=candidates,
  category.names = c('Adult','Recruit','Larvae'),
  filename='test.png',
  lwd=5,
  col=c("#762A83", "#5AAE61", "#FFED6F"),
  fill = c(alpha("#762A83",0.3), alpha("#5AAE61",0.3), alpha("#FFED6F",0.3)),
  cex = 3,
  cat.col = c("black", 'black', 'black'),
  cat.cex = 3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 180),
  cat.dist = c(0.08, 0.08, 0.05)
)

# try plotting heatmap of LFCs between treatment and control
all_genes=c(Adult_trmtGenes,Recruit_trmtGenes,Larv_trmtGenes)
all_genes=unique(all_genes)

adultLFC=as.data.frame(res_Adult$log2FoldChange,row.names = rownames(res_Adult))
recruitLFC=as.data.frame(res_Recruit$log2FoldChange,row.names = rownames(res_Recruit))
larvLFC=as.data.frame(Larv_res$trmt_LFC,row.names = rownames(Larv_res))

# 19169 total genes shared between adult/recruit and larval datasets

master_df=as.data.frame(merge(adultLFC,recruitLFC,by=0,all=TRUE))
rownames(master_df)=master_df$Row.names
master_df$Row.names=NULL
master_df=as.data.frame(merge(master_df,larvLFC,by=0,all=TRUE))
master_sub=subset(master_df,Row.names %in% all_genes)
master_order=master_sub[match(all_genes,master_sub$Row.names),]
rownames(master_order)=master_order$Row.names
master_order$Row.names=NULL

#col=color=colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.1)(100)
col=color = colorRampPalette(c('turquoise','turquoise','black','yellow','yellow'),bias=1.1)(50)
pheatmap(na.omit(master_order),clustering_distance_rows="correlation",cluster_rows = T,
         cluster_cols = F,show_rownames = F,labels_col = c('Adult','Recruit','Larvae'),
         color=col,fontsize = 20)
col=color = colorRampPalette(c('#21908dff','#21908dff','#21908dff','black','yellow','yellow','yellow'),bias=1.1)(50)
pheatmap(na.omit(master_order),cluster_rows = F,
         cluster_cols = F,show_rownames = F,labels_col = c('Adult','Recruit','Larvae'),
         color=col,fontsize = 20)


write.csv(as.data.frame(res_Adult),file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_Adult_trmt_response.csv')
write.csv(as.data.frame(res_Recruit),file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_Recruit_trmt_response.csv')
write.csv(as.data.frame(res_stage_control),file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_stage_control_AR_response.csv')


###### Generating directional GO output
res_Adult=read.csv('downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_Adult_trmt_response.csv',row.names = 1)
res_Recruit=read.csv('downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_Recruit_trmt_response.csv',row.names = 1)
Larv_res=read.csv('downstream_analyses/Mar2020/DESeq/DESeq2_June2021/Larv_DESeq2_VsdPvalsLFC_19July2021.csv',row.names = 1)

#Adult response
res_Adult$direction=ifelse(res_Adult$log2FoldChange>0,1,-1) #red=higher expression in heat
res_Adult$logP<-(-log((res_Adult$padj+0.0000000001),10))*(res_Adult$direction)

#Recruit response
res_Recruit$direction=ifelse(res_Recruit$log2FoldChange>0,1,-1)
res_Recruit$logP<-(-log((res_Recruit$padj+0.0000000001),10))*(res_Recruit$direction)

#Larval response
Larv_res$direction=ifelse(Larv_res$trmt_LFC>0,1,-1)
Larv_res$logP<-(-log((Larv_res$trmt_padj+0.0000000001),10))*(Larv_res$direction)

Adult_out<-as.data.frame(cbind("gene"=row.names(res_Adult),"logP"=res_Adult$logP))
Recruit_out<-as.data.frame(cbind("gene"=row.names(res_Recruit),"logP"=res_Recruit$logP))
Larv_out<-as.data.frame(cbind("gene"=row.names(Larv_res),"logP"=Larv_res$logP))

write.csv(Adult_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankHost_AdultResponseToTreat.csv",quote=F,row.names=F)
write.csv(Recruit_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankHost_RecruitResponseToTreat.csv",quote=F,row.names=F)
write.csv(Larv_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankHost_LarvResponseToTreat.csv",quote=F,row.names=F)

res_Adult=read.csv("downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_Adult_trmt_response.csv",row.names = 1)
res_Recruit=read.csv("downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_Recruit_trmt_response.csv", row.names = 1)



################### for the interaction of treatment and origin #######################
###### first try with interaction in formula -- go to ?results example #2 for how to extract results
design(dds) <- ~ stage + origin + treatment + origin:treatment
dds <- DESeq(dds)
resultsNames(dds)
res_int <- results(dds, name='originoff.treatmentHeat') # this is the interaction term --
int_genes=rownames(subset(res_int,padj<0.1)) # no genes for interaction


####### NOW combine the factors of interest into a single factor with all combinations of the original factors
dds$group <- factor(paste0(dds$origin, dds$treatment))
# change the design to include just this factor, e.g. ~ group
design(dds) <- ~ stage + group
dds <- DESeq(dds)
res_group=results(dds)
summary(res_group) #310 up, 353 down
resultsNames(dds)

res_in <- results(dds, contrast=c("group", "inHeat", "inControl")) # 116 up, 200 down
res_off <- results(dds, contrast=c("group", "offHeat", "offControl")) # 30 up, 13 down
res_origin_control <- results(dds, contrast=c('group','offControl','inControl')) # 95 up,171 down
res_origin_heat <- results(dds, contrast=c('group','offHeat','inHeat')) # 103 up, 111 down

candidates=list('inshore'=rownames(subset(res_in,padj<0.1)),'offshore'=rownames(subset(res_off,padj<0.1)))
venn(candidates)

##### normalize counts and save results
vsd=vst(dds,blind=FALSE)
norm_df=assay(vsd)

out=as.data.frame(cbind(norm_df,res_in$padj,res_off$padj,res_origin_control$padj,res_origin_heat$padj,res_in$log2FoldChange,res_off$log2FoldChange,res_origin_control$log2FoldChange,res_origin_heat$log2FoldChange))
names(out)[48:55]=c('in_padj','off_padj','origin_control_padj','origin_heat_padj','in_LFC','off_LFC','origin_control_LFC','origin_heat_LFC')

write.csv(out,file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2HostAR_origin_response_to_treat_VsdPvalsLFC_Oct2021.csv')

##### generating GO output
#inshore
res_in$direction=ifelse(res_in$log2FoldChange>0,1,-1) #red=higher expression in heat
res_in$logP<-(-log((res_in$padj+0.0000000001),10))*(res_in$direction)

#offshore
res_off$direction=ifelse(res_off$log2FoldChange>0,1,-1) #red=higher expression in offshore
res_off$logP<-(-log((res_off$padj+0.0000000001),10))*(res_off$direction)

in_out<-as.data.frame(cbind("gene"=row.names(res_in),"logP"=res_in$logP))
off_out<-as.data.frame(cbind("gene"=row.names(res_off),"logP"=res_off$logP))

write.csv(in_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankHost_AR_inshore_response_Oct2021.csv",quote=F,row.names=F)
write.csv(off_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankHost_AR_offshore_response_Oct2021.csv",quote=F,row.names=F)

# look at heatmaps
norm_df=assay(vsd)
means=apply(norm_df,1,mean) # means of rows
norm_df_scaled=norm_df-means #rescale expression data so it's up and down
df=as.data.frame(colData(dds)[,c("treatment","stage","origin")])
brewer.pal(n = 6, name = "Set2")
ann_colors=list(
  treatment=c(Control='#66C2A5', Heat='#FC8D62'),
  origin=c(off='#A6D854','in'='#FFD92F'),
  stage=c(Adult='#8DA0CB',Recruit='#E78AC3')
)

col=color = colorRampPalette(c('turquoise','turquoise','turquoise','black','yellow','yellow','yellow'),bias=1.7)(30)
pheatmap(norm_df_scaled[rownames(subset(res_in,padj<0.1)),],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment, not really stage

pheatmap(norm_df_scaled[rownames(subset(res_off,padj<0.1)),],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment, not really stage

#norm_df_order = norm_df_scaled[,c(1,5,9,13,36,40,44,3,7,11,36,38,42,46,2,6,10,14,37,41,45,4,8,12,
                                  #17,21,25,29,33,19,23,27,31,18,22,26,30,34,20,24,28,32)]
pheatmap(norm_df_order[rownames(subset(res_in,padj<0.1)),],color=col,cluster_cols=F,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment, not really stage





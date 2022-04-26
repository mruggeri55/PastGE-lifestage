library(DESeq2)
library(gplots)
library(pheatmap)
library(vsn)
library(RColorBrewer)

setwd('/Users/Maria/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/')

cts=read.table("computational/Dec2019/CG_AllCountsSym.txt",row.names=1) #Reading in the table of counts per isogroup by sample
head(cts) 
length(cts[,1])  #34160 isogroups

host=read.table("computational/CG_AllCountsHost_new.txt",header=TRUE,row.names=1)
colnames=colnames(host)
colnames(cts)=colnames 

names(cts)<-gsub("X","",names(cts))
names(cts)

#######################Creating table of conditions for your experiment -- origin, treatment, and life stage

# origin=treatment=stage=gt=c(1:length(names(cts)))
# stage[grep("A",names(cts))]="Adult"
# stage[grep("L",names(cts))]="Recruit"
# treatment[grep("C",names(cts))]="Control"
# treatment[grep("H",names(cts))]="Heat"
# gt=gsub("X|A|L|C|H","",names(cts))
# gts=unique(gt)#"10" "16" "18" "22" "29" "31" "35" "46" "47" "6"  "7"  "8" 
# ori=c('in','in','in','in','off','off','off','off','off','in','in','in')
# gts_ori=cbind(gts,ori)
# gt=as.data.frame(gt)
# names(gt)='gts'
# origin=merge(gt,gts_ori,by='gts',sort=FALSE)
# names(origin)=c('genotype','origin')
# 
# colData=data.frame(cbind(origin,treatment,stage))
# head(colData)

#write.csv(colData,file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/CG_colData_all.csv')
colData=read.csv('downstream_analyses/Mar2020/DESeq/DESeq2_June2021/CG_colData_all.csv')

################# To sum number of mapped reads
# 
# readsleft=c()
# for (column in names(cts)) {
#   val=sum(cts[,column])
#   readsleft=append(readsleft,val)}
# 
# RLtable=data.frame(cbind(names(cts),readsleft))
# RLtable$readsleft=as.numeric(as.character(RLtable$readsleft))
# write.csv(RLtable,"readsleft_CGsym.csv",quote=F) 

################ original DESeq methods for finding outliers

# real=newCountDataSet(cts,colData) 
# real=estimateSizeFactors(real)
# 
# ##### first, quality control #####
# 
# cds=estimateDispersions(real,method="blind")
# vsdBlind=varianceStabilizingTransformation(cds)
# 
# v="/Users/maria/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/downstream_analyses/Mar2020/DESeq/arrayQualityMetrics_CGsym/"
# 
# arrayQualityMetrics(vsdBlind,outdir=v,intgroup=c("genotype"),force=TRUE) #check .html output file in new folder

# remove 6, 22, 36, 46

counts.nobad1=cts[,-c(6,22,36,46)]  
conditions.nobad1=colData[-c(6,22,36,46),]

# real=newCountDataSet(counts.nobad1,conditions.nobad1) 
# real=estimateSizeFactors(real)
# 
# cds=estimateDispersions(real,method="blind")
# vsdBlind=varianceStabilizingTransformation(cds)
# 
# arrayQualityMetrics(vsdBlind,outdir=v,intgroup=c("genotype"),force=TRUE) #check .html output file in new folder

#other potential outliers aren't extreme, will retain remaining samples

#go ahead and remove low expression genes - those with count less than 2 in more than 90% of samples
counts.nobad1$low = apply(counts.nobad1[,1:44],1,function(x){sum(x<=2)})  #making new column counting number of samples with counts <=10 within each isogroup (host) <=5 within each dataset
counts.nobad1<-counts.nobad1[-which(counts.nobad1$low>40),1:44] #39.6 is 90% of 44 samples
nrow(counts.nobad1) #14124
head(counts.nobad1) 

counts=counts.nobad1
colData=conditions.nobad1

####################### starting DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ 1)
dds # 14124 isogroups

# differential expression analysis - DESeq will calculate the log2fold and Wald test p-value for last variable in design formula, but can pull out other contrasts later
# for now using NULL design to estimate size factors and dispersions and look at patterns in blind data
dds<-DESeq(dds)
vsd_null<-vst(dds,blind=TRUE)
plotDispEsts(dds)
plotPCA(vsd_null,intgroup=c("treatment",'stage')) # clearly separates by life stage, no real pattern by treatment
plotPCA(vsd_null,intgroup=c("origin",'treatment'))
pca_dat=plotPCA(vsd_null,intgroup=c("origin",'stage','treatment'),returnData=T) # woah really separates out by origin along PC2!!!
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

# sub out life stages to see if treatment gets resolved -- not really
vsd_null_adults<-vsd_null[,vsd_null$stage=="Adult"]
plotPCA(vsd_null_adults,intgroup=c('origin','treatment'))
vsd_null_recruits<-vsd_null[,vsd_null$stage=="Recruit"]
plotPCA(vsd_null_recruits,intgroup=c('origin','treatment'))

# take a look at PC3
library(ggfortify)
pca=prcomp(t(assay(vsd_null)))
pca=as.data.frame(pca$x)
colData$origin_treatment=paste(colData$origin,colData$treatment)
ggplot(pca,aes(x=PC1, y=PC3, colour=colData$origin_treatment,fill=colData$origin_treatment,shape=colData$stage)) +
  geom_point(size=3)+
  scale_shape_manual(values = c(21, 24)) + scale_fill_manual(values=c('dark blue','light blue','seagreen','darkseagreen1'))+
  scale_color_manual(values=c('dark blue','dark blue','seagreen','seagreen'))+labs(fill='treatment group')
# does not separate by treatment along PC3

############# ~ stage + origin + treatment ###########
design(dds) <- ~ stage + origin + treatment
dds <- DESeq(dds)
plotDispEsts(dds)

resultsNames(dds)
res <- results(dds,name="treatment_Heat_vs_Control") # this is the treatment effect bc last term in model
res_stage <- results(dds, name='stage_Recruit_vs_Adult')
res_origin <- results(dds, name='origin_off_vs_in' )

metadata(res)$filterTheta #0.116

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
  filename='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/SYM_Venn_stage_origin_trmt.png',
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 2,
  cat.col = c("#440154ff", '#21908dff', 'black'),
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 180),
  cat.dist = c(0.08, 0.08, 0.05)
)

vsd1=vst(dds, blind=FALSE)
plotPCA(vsd1,intgroup=c('stage','origin')) #looks the same as blind = TRUE, now PC1 explains a bit more variance (30% with blind=F vs 28% with blind=T)
plotPCA(vsd1[rownames(vsd1) %in% sigTrmt,],intgroup=c('treatment','origin'))
plotPCA(vsd1[rownames(vsd1) %in% sigTrmt,],intgroup=c('treatment','stage')) # treatment genes still separate by life stage
plotPCA(vsd1[rownames(vsd1) %in% sigTrmtOnly,],intgroup=c('treatment','stage'))

pca_dat=plotPCA(vsd1[rownames(vsd1) %in% sigTrmtOnly,],intgroup=c('stage','treatment'),returnData=T)
ggplot(pca_dat,aes(x=group,y=PC1,fill=group))+
  geom_boxplot(alpha=0.8)+theme_classic()+
  scale_fill_manual(values=c('dark blue','light blue','seagreen','darkseagreen1'))+
  theme(legend.position = "none")+
  xlab('') + ylab('PC1')+
  theme(axis.text.x=element_text(size=15,face='bold'),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16))

##### writing out results
norm_df=assay(vsd1)
out=as.data.frame(cbind(norm_df,res$padj,res_stage$padj,res_origin$padj,res$log2FoldChange,res_stage$log2FoldChange,res_origin$log2FoldChange))
names(out)[45:50]=c('trmt_padj','stage_padj','origin_padj','trmt_LFC','stage_LFC','origin_LFC')

write.csv(out,file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_SYM_vsd_pvals_July2021.csv')

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

write.csv(trmt_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankSYM_ARtrmt.csv",quote=F,row.names=F)
write.csv(origin_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankSYM_ARorigin.csv",quote=F,row.names=F)
write.csv(stage_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankSYM_ARstage.csv",quote=F,row.names=F)


# let's take a look at some heatmaps
# scale values
means=apply(norm_df,1,mean) # means of rows
norm_df_scaled=norm_df-means #rescale expression data so it's up and down
df=as.data.frame(colData(dds)[,c("treatment","origin",'stage')])
display.brewer.pal(n = 6, name = 'Set2')
brewer.pal(n = 6, name = "Set2")
ann_colors=list(
  treatment=c(Control='#66C2A5', Heat='orange'),
  origin=c(off='#A6D854','in'='#FFD92F'),
  stage=c(Adult='#8DA0CB',Recruit='#E78AC3')
)

# treatment DEGs
col=color= colorRampPalette(rev(c('turquoise','turquoise','black','yellow','yellow')),bias=0.8)(30)
pheatmap(norm_df_scaled[sigTrmt,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment then stage for all sig trmt DEGs
#treatment only
col=color= colorRampPalette(rev(c('turquoise','turquoise','black','yellow','yellow')),bias=1.0)(30)
pheatmap(norm_df_scaled[sigTrmtOnly,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment then stage for all sig trmt DEGs

# origin DEGs
col=color= colorRampPalette(rev(c('turquoise','turquoise','black','yellow','yellow')),bias=1.1)(30)
pheatmap(norm_df_scaled[sigOrigin,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment then stage for all sig trmt DEGs

# stage DEGs
col=color= colorRampPalette(rev(c('turquoise','turquoise','black','yellow','yellow')),bias=1.1)(30)
pheatmap(norm_df_scaled[sigStage,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment then stage for all sig trmt DEGs

#### now group stage and treatment to see which genes are DE by life stage due to treatment ####
design(dds) <- ~ origin + stage + treatment + stage:treatment
dds_int <- DESeq(dds)
res_int <- results(dds_int, name='stageRecruit.treatmentHeat') # no genes sig

# combine the factors of interest into a single factor with all combinations of the original factors
dds$group <- factor(paste0(dds$stage, dds$treatment))
# change the design to include just this factor, e.g. ~ group
design(dds) <- ~ origin + group
dds <- DESeq(dds)
res_group=results(dds)
summary(res_group)
resultsNames(dds)

res_Adult <- results(dds, contrast=c("group", "AdultHeat", "AdultControl")) # 56 up, 36 down
res_Recruit <- results(dds, contrast=c("group", "RecruitHeat", "RecruitControl")) # 100 up, 55 down
res_stage_control <- results(dds, contrast=c('group', 'AdultControl', 'RecruitControl')) #1075 up, 4379 down
res_stage_heat <- results(dds, contrast=c('group','AdultHeat','RecruitHeat')) #819 up, 2530 down
res_origin <- results(dds,contrast=c('origin','in','off')) #664 up, 1016 down

write.csv(as.data.frame(res_Adult),file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_symAdult_trmt_response.csv')
write.csv(as.data.frame(res_Recruit),file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_symRecruit_trmt_response.csv')
write.csv(as.data.frame(res_stage_control),file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_sym_stage_control_AR.csv') # adults / recruits

save(dds, res_Adult,res_Recruit, res_stage_control, res_stage_heat, res_origin, file = 'DESeq2_sym_by_stage.RData')
load(file = 'downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_sym_by_stage.RData')

###### Generating directional GO output using pvals
#Adult response
res_Adult$direction=ifelse(res_Adult$log2FoldChange>0,1,-1) #red=higher expression in heat
res_Adult$logP<-(-log((res_Adult$padj+0.0000000001),10))*(res_Adult$direction)
#Recruit response
res_Recruit$direction=ifelse(res_Recruit$log2FoldChange>0,1,-1)
res_Recruit$logP<-(-log((res_Recruit$padj+0.0000000001),10))*(res_Recruit$direction)
#output files
Adult_out<-as.data.frame(cbind("gene"=row.names(res_Adult),"logP"=res_Adult$logP))
Recruit_out<-as.data.frame(cbind("gene"=row.names(res_Recruit),"logP"=res_Recruit$logP))
write.csv(Adult_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankSYM_AdultResponseToTreat.csv",quote=F,row.names=F)
write.csv(Recruit_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankSYM_RecruitResponseToTreat.csv",quote=F,row.names=F)

#### generating GO output using LFCs
Adult_out<-as.data.frame(cbind("gene"=row.names(res_Adult),"LFC"=res_Adult$log2FoldChange))
Recruit_out<-as.data.frame(cbind("gene"=row.names(res_Recruit),"LFC"=res_Recruit$log2FoldChange))
write.csv(Adult_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_LFC/GO_LFC_SYM_AdultResponseToTreat.csv",quote=F,row.names=F)
write.csv(Recruit_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_LFC/GO_LFC_SYM_RecruitResponseToTreat.csv",quote=F,row.names=F)



##### moving on
Adult_trmtGenes=rownames(subset(res_Adult,padj < 0.1))
Recruit_trmtGenes=rownames(subset(res_Recruit,padj < 0.1))

candidates=list('Adult'=Adult_trmtGenes,'Recruit'=Recruit_trmtGenes)
venn(candidates)

#Adult_Control='#762A83', Adult_Heat='#AF8DC3',Recruit_Control="#1B7837",Recruit_Heat="#7FBF7B"
library(VennDiagram)
venn.diagram(
  x=candidates,
  category.names = c('Adult','Recruit'),
  filename='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/plots/SYM_TrmtResponse_byLifeStageVenn.png',
  col=c('#762A83',"#1B7837"),
  lwd=5,
  margin= 0.2,
  fill = c(alpha("#762A83",0.3), alpha('#1B7837',0.3)),
  cex = c(2.5,2,3),
  cat.col = 'black',
  cat.cex = 2.5,
  cat.default.pos = "outer",
  cat.pos = c(-180, 0),
  cat.dist = c(0.07, 0.07),
  rotation.degree=90
)


# make output of annotated genes overlapping and their LFC
overlap=Adult_trmtGenes[Adult_trmtGenes %in% Recruit_trmtGenes]
iso2gene=read.table('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/kb8_blast_annot_July2021/kb8_iso2gene.tab',sep = '\t', colClasses = "character")
overlap_genes=iso2gene[iso2gene$V1 %in% overlap,]
AdultSYM_LFC = subset(res_Adult, subset = rownames(res_Adult) %in% overlap_genes$V1, select=log2FoldChange)
RecruitSYM_LFC = subset(res_Recruit, subset = rownames(res_Recruit) %in% overlap_genes$V1, select=log2FoldChange)
all_order=overlap_genes[match(rownames(AdultSYM_LFC),overlap_genes$V1),]
merge_df=as.data.frame(cbind(all_order,AdultSYM_LFC,RecruitSYM_LFC))
merge_df$gene=as.character(lapply(strsplit(merge_df$V2, split='OS='),'[',1))
merge_df$GN=as.character(lapply(strsplit(merge_df$V2, split='GN='),'[',2))
merge_df$GN=as.character(lapply(strsplit(merge_df$GN, split='PE='),'[',1))
merge_df$V2=NULL

head(merge_df)
names(merge_df)=c('isogroup','adult sym (LFC)', 'recruit sym (LFC)', 'gene','GN')
#write.csv(x=merge_df,file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/interesting_genes/SYM_stage_CoreGenes.csv')


# how many of these overlap with full dataset treatment genes?
AandR_genes=c(Adult_trmtGenes,Recruit_trmtGenes)
AandR_UniqGenes=unique(AandR_genes)
AandR_nodups=AandR_genes[!(duplicated(AandR_genes) | duplicated(AandR_genes, fromLast = TRUE))]

candidates=list('Adult and recruit responses'=AandR_UniqGenes,'Treatment'=sigTrmt)
venn(candidates)

candidates=list('Adult'=Adult_trmtGenes,'Recruit'=Recruit_trmtGenes, 'Treatment'=sigTrmt)
venn(candidates)

# try plotting heatmap of LFCs between treatment and control
adultLFC=as.data.frame(res_Adult$log2FoldChange,row.names = rownames(res_Adult))
recruitLFC=as.data.frame(res_Recruit$log2FoldChange,row.names = rownames(res_Recruit))
master_df=as.data.frame(merge(adultLFC,recruitLFC,by=0,all=TRUE))
master_sub=subset(master_df,Row.names %in% AandR_nodups)
master_order=master_sub[match(AandR_nodups,master_sub$Row.names),]
rownames(master_order)=master_order$Row.names
master_order$Row.names=NULL

#col=color=colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.1)(100)
col=color = colorRampPalette(c('turquoise','turquoise','black','yellow','yellow'),bias=1.0)(50)
pheatmap(na.omit(master_order),clustering_distance_rows="correlation",cluster_rows = T,cluster_cols = F,show_rownames = F,labels_col = c('Adult','Recruit'),color=col)

pheatmap(na.omit(master_order),cluster_rows = F,cluster_cols = F,show_rownames = F,labels_col = c('Adult','Recruit'),color=col)

# let's look at eigengene
vsd2=vst(dds,blind = F)
pca_dat=plotPCA(vsd2[rownames(vsd2) %in% Recruit_trmtGenes,],intgroup=c('stage','treatment'),returnData=T)
ggplot(pca_dat,aes(x=group,y=PC1,fill=group))+
  geom_boxplot(alpha=0.8)+theme_classic()+
  scale_fill_manual(values=c('dark blue','light blue','seagreen','darkseagreen1'))+
  theme(legend.position = "none")+
  xlab('') + ylab('PC1')+
  theme(axis.text.x=element_text(size=15,face='bold'),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16))

# heatmap with all samples
norm_df=assay(vsd2)
means=apply(norm_df,1,mean) # means of rows
norm_df_scaled=norm_df-means #rescale expression data so it's up and down
colnames(norm_df_scaled)=gsub('CL','RC',colnames(norm_df_scaled))
colnames(norm_df_scaled)=gsub('HL','RH',colnames(norm_df_scaled))

df=as.data.frame(colData(dds)[,c("treatment","origin",'stage')])
df$stage_treatment=paste(df$stage, df$treatment, sep='_')
df$treatment=NULL
df$stage=NULL

bleach_dat=read.csv('downstream_analyses/CG_Larv_WGCNA_traits.csv')
bleach_dat=subset(bleach_dat,select=c('colony','bleachstatus'))
bleach_dat$colony=gsub('X','',bleach_dat$colony)
bleach_order=bleach_dat[match(rownames(df),bleach_dat$colony),]
df$bleach_score=bleach_order$bleachstatus
rownames(df)=gsub('CL','RC',rownames(df))
rownames(df)=gsub('HL','RH',rownames(df))

display.brewer.pal(n = 6, name = 'Set2')
brewer.pal(n = 6, name = "Set2")
display.brewer.pal(n=6,name='YlOrBr')
bleach=brewer.pal(n = 6, name = "YlOrBr")
brewer.pal(n=6,name='PRGn')

# ann_colors=list(
#   treatment=c(Control='#66C2A5', Heat='orange'),
#   origin=c(off='#A6D854','in'='#E78AC3'),
#   stage=c(Adult='#8DA0CB',Recruit='#FFD92F'),
#   bleach_score=bleach
# )

df$origin=NULL
ann_colors=list(
  bleach_score=bleach,
  stage_treatment=c(Adult_Control='#762A83', Adult_Heat='#AF8DC3',Recruit_Control="#1B7837",Recruit_Heat="#7FBF7B")
  #origin=c(off='aquamarine3','in'='orange'),
)

# treatment DEGs
col=color= colorRampPalette(rev(c('yellow','yellow','black','turquoise','turquoise')),bias=0.8)(30)
pheatmap(norm_df_scaled[AandR_nodups,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) 

#put in order of stage then treatment then origin then bleach score
norm_df_scaled<-norm_df_scaled[,c(1,5,8,12,16,20,23,27,31,34,38,42,2,9,13,17,24,28,32,35,39,3,6,10,14,18,21,25,29,33,36,40,43,4,7,11,15,19,22,26,30,37,41,44)]
norm_df_scaled<-norm_df_scaled[,c(1,2,3,4,10,11,12,5,6,7,8,9,13,14,15,20,21,16,17,18,19,22,23,24,25,31,32,33,26,27,28,29,30,34,35,36,37,42,43,44,38,39,40,41)]
norm_df_scaled<-norm_df_scaled[,c(1,3,4,5,6,7,9,10,11,12,8,2,17,14,16,18,19,13,15,20,21,22:44)]

pheatmap(norm_df_scaled[AandR_nodups,],color=col,cluster_cols=F,cluster_rows = T,show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors, gaps_col = c(12,21,33)) 
pheatmap(norm_df_scaled[AandR_UniqGenes,],
         color=col,cluster_cols=F,
         cluster_rows = T,
         show_rownames=FALSE,
         annotation_col=df,
         annotation_colors = ann_colors, 
         fontsize = 12,
         gaps_col = c(12,21,33)) 

#pheatmap(norm_df_scaled[photo_genes,],color=col,cluster_cols=F,cluster_rows = T,show_rownames=T,annotation_col=df,annotation_colors = ann_colors, gaps_col = c(12,21,33)) 

pheatmap(norm_df_scaled[AandR_nodups,],color=col,cluster_cols=F,cluster_rows = F,show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors, gaps_col = c(12,21,33)) 

col=color= colorRampPalette(rev(c('yellow','yellow','black','turquoise','turquoise')),bias=0.8)(30)
pheatmap(norm_df_scaled[Adult_trmtGenes,],color=col,cluster_cols=F,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) 

col=color= colorRampPalette(rev(c('yellow','yellow','black','turquoise','turquoise')),bias=0.8)(30)
pheatmap(norm_df_scaled[Recruit_trmtGenes,],color=col,cluster_cols=F,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) 

stage_origin_trmt=read.csv('downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_SYM_vsd_pvals_July2021.csv',row.names = 1)
sigTrmtOnly=rownames(subset(stage_origin_trmt,trmt_padj < 0.1 & !stage_padj < 0.1 & !origin_padj < 0.1))
col=color= colorRampPalette(rev(c('yellow','yellow','black','turquoise','turquoise')),bias=1)(30)
pheatmap(norm_df_scaled[sigTrmtOnly,],color=col,cluster_cols=F,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors,gaps_col = c(12,21,33)) 

col=color= colorRampPalette(rev(c('yellow','yellow','black','turquoise','turquoise')),bias=0.7)(30)
sigTrmtAll=rownames(subset(stage_origin_trmt,trmt_padj < 0.1))
pheatmap(norm_df_scaled[sigTrmtAll,],color=col,cluster_cols=F,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors,gaps_col = c(12,21,33)) 
pheatmap(norm_df_scaled[sigTrmtAll,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors,gaps_col = c(12,21,33)) 

stage_heat_genes=rownames(subset(res_stage_heat, padj < 0.1))
col=color= colorRampPalette(rev(c('turquoise','turquoise','black','yellow','yellow')),bias=1.3)(30)
pheatmap(norm_df_scaled[stage_heat_genes,],color=col,cluster_cols=F,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) 

candidates=list('Adult'=Adult_trmtGenes,'Recruit'=Recruit_trmtGenes,'Treat DEGs Only'=sigTrmtOnly)
venn(candidates)

candidates=list('Adult'=Adult_trmtGenes,'Recruit'=Recruit_trmtGenes,'Treat DEGs All'=sigTrmtAll)
venn(candidates)

################ Barshis et al., 2013 chi2 test #################
# now try Barhis plots -- my thoughts are this is a detection thing, and more treatment responsive genes actually overlap between adults and recruits
#up_Rgenes=rownames(res_Recruit[which(res_Recruit$log2FoldChange > 0 & res_Recruit$padj < 0.1),])
RgenesOnly=Recruit_trmtGenes[! Recruit_trmtGenes %in% Adult_trmtGenes] #139

res_A_Rgenes=subset(res_Adult,rownames(res_Adult) %in% RgenesOnly)
res_R_sub=res_Recruit[rownames(res_Recruit) %in% RgenesOnly,]
ggplot()+geom_point(aes(x=res_R_sub$log2FoldChange,y=res_A_Rgenes$log2FoldChange))+
  geom_abline(intercept = 0,slope=1,linetype='dashed',color='blue')+ylab('Adult response log2(fold change)')+
  xlab('Recruit response log2(fold change)')+ggtitle('Recruit DEGs to treatment')+xlim(-3,3)+ylim(-3,3)

# same thing but for adult degs to trmt
AgenesOnly=Adult_trmtGenes[! Adult_trmtGenes %in% Recruit_trmtGenes] #76

res_R_Agenes=subset(res_Recruit,rownames(res_Recruit) %in% AgenesOnly)
res_A_sub=res_Adult[rownames(res_Adult) %in% AgenesOnly,]
ggplot()+geom_point(aes(x=res_A_sub$log2FoldChange,y=res_R_Agenes$log2FoldChange))+
  geom_abline(intercept = 0,slope=1,linetype='dashed',color='blue')+ylab('Recruit response log2(fold change)')+
  xlab('Adult response log2(fold change)')+ggtitle('Adult DEGs to treatment')+ylim(-3,3)+ylim(-3,3)


# check for frontloading of recruit DEGs in Adults
res_stage_control=results(dds,contrast = c('group','AdultControl','RecruitControl'))
RgenesOnlyUP=rownames(subset(res_Recruit,padj < 0.1 & log2FoldChange > 0  & !rownames(res_Recruit) %in% Adult_trmtGenes))

res_stage_control_rDEGs=subset(res_stage_control,rownames(res_stage_control) %in% RgenesOnlyUP)
res_stage_heat_rDEGs=subset(res_stage_heat,rownames(res_stage_heat) %in% RgenesOnlyUP)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_A_Rgenes=subset(res_Adult,rownames(res_Adult) %in% RgenesOnlyUP)
res_R_Rgenes=subset(res_Recruit,rownames(res_Recruit) %in% RgenesOnlyUP)
rownames(res_A_Agenes) == rownames(res_R_Agenes)
RvA_response = res_A_Rgenes$log2FoldChange / res_R_Rgenes$log2FoldChange
ggplot()+geom_point(aes(x=RvA_response,y=res_stage_control_rDEGs$log2FoldChange))+
  geom_hline(aes(x=RvA_response,y=res_stage_control_inDEGs$log2FoldChange),yintercept=1,linetype='dashed',color='red')+
  geom_vline(aes(x=RvA_response,y=res_stage_control_inDEGs$log2FoldChange),xintercept=1,linetype='dashed',color='red')+
  ylab('Adult vs Recruit control log2(fold change)')+
  xlab('Adult FC vs Recruit FC')+ggtitle('Recruit DEGs to treatment ONLY UP')

###### check for frontloading of Adult DEGs in Recruits
res_stage_control=results(dds,contrast = c('group','RecruitControl','AdultControl'))
AgenesOnlyUP=rownames(subset(res_Adult,padj < 0.1 & log2FoldChange > 0  & !rownames(res_Adult) %in% Recruit_trmtGenes))

res_stage_control_ADEGs=subset(res_stage_control,rownames(res_stage_control) %in% AgenesOnlyUP)
res_stage_heat_ADEGs=subset(res_stage_heat,rownames(res_stage_heat) %in% AgenesOnlyUP)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_R_Agenes=subset(res_Recruit,rownames(res_Recruit) %in% AgenesOnlyUP)
res_A_Agenes=subset(res_Adult,rownames(res_Adult) %in% AgenesOnlyUP)
rownames(res_R_Agenes) == rownames(res_A_Agenes)
AvR_response = res_R_Agenes$log2FoldChange / res_A_Agenes$log2FoldChange
ggplot()+geom_point(aes(x=AvR_response,y=res_stage_control_ADEGs$log2FoldChange))+
  geom_hline(aes(x=AvR_response,y=res_stage_control_ADEGs$log2FoldChange),yintercept=1,linetype='dashed',color='red')+
  geom_vline(aes(x=AvR_response,y=res_stage_control_ADEGs$log2FoldChange),xintercept=1,linetype='dashed',color='red')+
  ylab('Recruit vs Adult control log2(fold change)')+
  xlab('Recruit FC vs Adult FC')+ggtitle('Adult DEGs to treatment ONLY')
########################################

# let's check if there are any recruit DEGs being frontloaded
# aka look at if any of the treatment responsive genes not responding in recruits are upregulated relative to adult samples in control conditions
res_stage_control=results(dds,contrast = c('group','RecruitControl','AdultControl'))
up_control=subset(res_stage_control,padj < 0.1 & log2FoldChange > 0)
stage_control_up=rownames(up_control)

candidates=list('adult heat response up' = AgenesOnlyUP,'recruit vs adult control up' = stage_control_up)
venn(candidates)

# now see if adults are frontloading recruit responsive genes
res_stage_control=results(dds,contrast = c('group','AdultControl','RecruitControl'))
up_control=subset(res_stage_control,padj < 0.1 & log2FoldChange > 0)
stage_control_up=rownames(up_control)

candidates=list('recruit heat response up' = RgenesOnly,'adult vs recruit control up' = stage_control_up)
venn(candidates)





##### ~ stage + origin*treatment #####
design(dds) <- ~ stage + origin + treatment + origin:treatment
dds_int <- DESeq(dds)
res_int <- results(dds_int, name="originoff.treatmentHeat") # no sig genes

# combine the factors of interest into a single factor with all combinations of the original factors
dds$group <- factor(paste0(dds$origin, dds$treatment))
# change the design to include just this factor, e.g. ~ group
design(dds) <- ~ stage + group
dds <- DESeq(dds)
res_group=results(dds)
summary(res_group)
resultsNames(dds)

res_in <- results(dds, contrast=c("group", "inHeat", "inControl")) # 100 up, 42 down
res_off <- results(dds, contrast=c("group", "offHeat", "offControl")) # 30 up, 27 down
res_origin_control <- results(dds, contrast=c('group', 'inControl', 'offControl')) #351 up, 487 down
res_origin_heat <- results(dds, contrast=c('group','inHeat','offHeat')) #298 up, 310 down
res_stage <- results(dds,contrast=c('stage','Adult','Recruit')) #1351 up, 5788 down

write.csv(as.data.frame(res_in),file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_symInshore_trmt_response.csv')
write.csv(as.data.frame(res_off),file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_symOffshore_trmt_response.csv')
write.csv(as.data.frame(res_origin_control),file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_sym_origin_control_IO.csv') # inshore / offshore

save(dds, res_in,res_off, res_origin_control, res_origin_heat, res_stage, file = 'downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_sym_by_origin.RData')
load('downstream_analyses/Mar2020/DESeq/DESeq2_June2021/DESeq2_sym_by_origin.RData')

###### Generating directional GO output
# #inshore response
# res_in$direction=ifelse(res_in$log2FoldChange>0,1,-1) #red=higher expression in heat
# res_in$logP<-(-log((res_in$padj+0.0000000001),10))*(res_in$direction)
# 
# #off response
# res_off$direction=ifelse(res_off$log2FoldChange>0,1,-1)
# res_off$logP<-(-log((res_off$padj+0.0000000001),10))*(res_off$direction)
# 
# in_out<-as.data.frame(cbind("gene"=row.names(res_in),"logP"=res_in$logP))
# off_out<-as.data.frame(cbind("gene"=row.names(res_off),"logP"=res_off$logP))
# 
# write.csv(in_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankSYM_InshoreResponseToTreat.csv",quote=F,row.names=F)
# write.csv(off_out,file="downstream_analyses/Mar2020/DESeq/DESeq2_June2021/GO_MWU-master/GO_input_files/GOrankSYM_OffshoreResponseToTreat.csv",quote=F,row.names=F)

#### moving on
in_trmtGenes=rownames(subset(res_in,padj < 0.1))
off_trmtGenes=rownames(subset(res_off,padj < 0.1))

candidates=list('Inshore'=in_trmtGenes,'Offshore'=off_trmtGenes)
venn(candidates)

candidates=list('inshore response'=in_trmtGenes,'offshore response'=off_trmtGenes)
library(VennDiagram)
venn.diagram(
  x=candidates,
  category.names = c('inshore','offshore'),
  filename='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/plots/SYM_TrmtResponse_byOrigin_Venn.png',
  margin=0.15,
  lwd=5,
  col=c("orange","cyan4"),
  fill = c(alpha("orange",0.5), alpha('cyan4',0.5)),
  cex = c(3,2,2.2),
  cat.col = 'black',
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(-90, 90),
  cat.dist = c(0.12, 0.12)
)

#make output of overlapping genes
overlap=in_trmtGenes[in_trmtGenes %in% off_trmtGenes]
iso2gene=read.table('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/kb8_blast_annot_July2021/kb8_iso2gene.tab',sep = '\t', colClasses = "character")
overlap_genes=iso2gene[iso2gene$V1 %in% overlap,]
InshoreSYM_LFC = subset(res_in, subset = rownames(res_in) %in% overlap_genes$V1, select=log2FoldChange)
OffshoreSYM_LFC = subset(res_off, subset = rownames(res_off) %in% overlap_genes$V1, select=log2FoldChange)
all_order=overlap_genes[match(rownames(InshoreSYM_LFC),overlap_genes$V1),]
merge_df=as.data.frame(cbind(all_order,InshoreSYM_LFC,OffshoreSYM_LFC))
merge_df$gene=as.character(lapply(strsplit(merge_df$V2, split='OS='),'[',1))
merge_df$GN=as.character(lapply(strsplit(merge_df$V2, split='GN='),'[',2))
merge_df$GN=as.character(lapply(strsplit(merge_df$GN, split='PE='),'[',1))
merge_df$V2=NULL

head(merge_df)
names(merge_df)=c('isogroup','inshore sym (LFC)', 'offshore sym (LFC)', 'gene','GN')
write.csv(x=merge_df,file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/interesting_genes/SYM_origin_CoreGenes.csv')


# plot out LFCs
InandOff_genes=c(in_trmtGenes,off_trmtGenes)
InandOff_UniqGenes=unique(InandOff_genes)
InandOff_nodups=InandOff_genes[!(duplicated(InandOff_genes) | duplicated(InandOff_genes, fromLast = TRUE))]

# try plotting heatmap of LFCs between treatment and control
inLFC=as.data.frame(res_in$log2FoldChange,row.names = rownames(res_in))
offLFC=as.data.frame(res_off$log2FoldChange,row.names = rownames(res_off))
master_df=as.data.frame(merge(inLFC,offLFC,by=0,all=TRUE))
master_sub=subset(master_df,Row.names %in% InandOff_nodups)
master_order=master_sub[match(InandOff_nodups,master_sub$Row.names),]
rownames(master_order)=master_order$Row.names
master_order$Row.names=NULL

#col=color=colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.1)(100)
col=color = colorRampPalette(c('turquoise','turquoise','black','yellow','yellow'),bias=1.0)(50)
pheatmap(na.omit(master_order),clustering_distance_rows="correlation",cluster_rows = T,cluster_cols = F,show_rownames = F,labels_col = c('in','off'),color=col)

#ordered by sig DEGs
pheatmap(na.omit(master_order),cluster_rows = F,cluster_cols = F,show_rownames = F,labels_col = c('in','off'),color=col)

# let's look at eigengene
vsd3=vst(dds,blind = F)
#pca_dat=plotPCA(vsd3[rownames(vsd3) %in% off_trmtGenes,],intgroup=c('origin','treatment'),returnData=T)
pca_dat=plotPCA(vsd3[rownames(vsd3) %in% in_trmtGenes,],intgroup=c('origin','treatment'),returnData=T)
#pca_dat=plotPCA(vsd3[rownames(vsd3) %in% InandOff_nodups,],intgroup=c('origin','treatment'),returnData=T)
ggplot(pca_dat,aes(x=group,y=PC1,fill=group))+
  geom_boxplot(alpha=0.8)+theme_classic()+
  scale_fill_manual(values=c('dark blue','light blue','seagreen','darkseagreen1'))+
  theme(legend.position = "none")+
  xlab('') + ylab('PC1')+
  theme(axis.text.x=element_text(size=15,face='bold'),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16))

# heatmap with all samples
norm_df=assay(vsd3)
means=apply(norm_df,1,mean) # means of rows
norm_df_scaled=norm_df-means #rescale expression data so it's up and down
df=as.data.frame(colData(dds)[,c("treatment","origin",'stage')])
display.brewer.pal(n = 6, name = 'Set2')
brewer.pal(n = 6, name = "Set2")
ann_colors=list(
  treatment=c(Control='#66C2A5', Heat='orange'),
  origin=c(off='#A6D854','in'='#FFD92F'),
  stage=c(Adult='#8DA0CB',Recruit='#E78AC3')
)

# treatment DEGs
col=color= colorRampPalette(rev(c('turquoise','turquoise','black','yellow','yellow')),bias=0.8)(30)
pheatmap(norm_df_scaled[InandOff_nodups,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) 

#put in order of sig DEGs by stage
norm_df_sub=subset(norm_df_scaled,rownames(norm_df_scaled) %in% InandOff_nodups)
norm_order=norm_df_sub[match(InandOff_nodups,rownames(norm_df_sub)),]
pheatmap(norm_order,color=col,cluster_cols=T,cluster_rows = F,show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) 

col=color= colorRampPalette(rev(c('turquoise','turquoise','black','yellow','yellow')),bias=0.8)(30)
pheatmap(norm_df_scaled[in_trmtGenes,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment then stage for all sig trmt DEGs

col=color= colorRampPalette(rev(c('turquoise','turquoise','black','yellow','yellow')),bias=0.8)(30)
pheatmap(norm_df_scaled[off_trmtGenes,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment then stage for all sig trmt DEGs

origin_heat_genes=rownames(subset(res_origin_heat, padj < 0.1))
col=color= colorRampPalette(rev(c('turquoise','turquoise','black','yellow','yellow')),bias=1.3)(30)
pheatmap(norm_df_scaled[origin_heat_genes,],color=col,cluster_cols=F,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment then stage for all sig trmt DEGs








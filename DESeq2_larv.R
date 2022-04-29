library(DESeq2)
library(gplots)
library(pheatmap)
library(vsn)
library(RColorBrewer)
library(ggplot2)

setwd('/Users/Maria/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/')

cts=read.table("computational/Larv_AllCountsHost_new.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample
head(cts) 
length(cts[,1])  #31944 isogroups
names(cts)<-gsub("X","",names(cts))
names(cts)

#######################Creating table of conditions for your experiment -- origin, treatment, and life stage
origin=treatment=gt=c(1:length(names(cts)))
treatment[grep("C",names(cts))]="Control"
treatment[grep("H",names(cts))]="Heat"
gt=gsub("X|A|L|C|H","",names(cts))
gts=unique(gt)#"15" "26" "27" "28" "29" "2"  "4" 
ori=c('in','off','off','off','off','in','in')
gts_ori=cbind(gts,ori)
gt=as.data.frame(gt)
names(gt)='gts'
origin=merge(gt,gts_ori,by='gts',sort=FALSE)
names(origin)=c('genotype','origin')

colData=data.frame(cbind(origin,treatment),row.names = colnames(cts))
colData

#write.csv(colData,file='downstream_analyses/Mar2020/DESeq/input_files/Larv_colData.csv')

############################################ used original DESeq methods to look for outliers
#library(DESeq)
#real=newCountDataSet(cts,colData) 
#real=estimateSizeFactors(real)

#####first, quality control

#cds=estimateDispersions(real,method="blind")
#vsdBlind=varianceStabilizingTransformation(cds)

#v="/Users/maria/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/downstream_analyses/Mar2020/DESeq/"

#arrayQualityMetrics(vsdBlind,outdir=v,intgroup=c("genotype"),force=TRUE) #check .html output file in new folder

# 27H and 29C came up as outliers according to sample distance but no other outlier criteria
# let's keep all samples

#go ahead and remove low expression genes - those with count less than 2 in more than 90% of samples
cts$low = apply(cts[,1:14],1,function(x){sum(x<=2)})  #making new column counting number of samples with counts <=10 within each isogroup (host) <=5 within each dataset
cts<-cts[-which(cts$low>13),1:14] #13 is 90% of 14 samples
nrow(cts) #23806
head(cts)
counts=cts

###################### starting DESeq2 analysis -- exploring data with no design info
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ 1)
dds

# differential expression analysis - DESeq will calculate the log2fold and Wald test p-value for last variable in design formula
vsd_null<-vst(dds,blind=TRUE)
#col=c('blue','red','light blue','pink')
pca_dat=plotPCA(vsd_null,intgroup=c("origin",'treatment'),returnData=T)
percentVar <- round(100 * attr(pca_dat, "percentVar"))
pca_dat$origin_treatment=as.factor(paste(pca_dat$origin,pca_dat$treatment,sep=' '))
larvae=ggplot(pca_dat, aes(x=PC1, y=PC2,colour=origin_treatment,fill=origin_treatment)) +
  geom_point(size = 3,stroke=1,shape=22) +
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
  theme_light() + theme(text=element_text(size=16),legend.position = 'none')+ggtitle('Larvae')


# take a look at PC3
library(ggfortify)
pca=prcomp(t(assay(vsd_null)),scale. = TRUE)
pca=as.data.frame(pca$x)
ggplot(pca,aes(x=PC3, y=PC4, colour=colData$treatment,fill=colData$treatment)) +
  geom_point(size=3)
# does not separate by treatment along PC3

####### now moving on to actual differential expression analysis by treatment and origin ##########
design(dds) <- ~ origin + treatment
# set control treatment as reference level
dds$treatment <- relevel(dds$treatment, ref = "Control")
# calculate log2fold and p values for all genes and extract results
dds <- DESeq(dds)
plotDispEsts(dds)

resultsNames(dds)
res <- results(dds) # this is the treatment effect bc last term in model
res_treatment <- results(dds, name="treatment_Heat_vs_Control")
summary(res_treatment) #118 up, 165 down
metadata(res_treatment)$filterTheta #0.43

plot(metadata(res_treatment)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res_treatment)$lo.fit, col="red")
abline(v=metadata(res_treatment)$filterTheta)
# ok thats understandable

res_origin <- results(dds, name="origin_off_vs_in")
summary(res_origin) # 98 up, 127 down
metadata(res_origin)$filterTheta #0.21

# lots of low count genes in both 43% for treatment, 21% for origin

# subset out significant gene sets
p.adj.cutoff=0.1
# add T/F column if gene is sig
res_treatment$threshold <- as.logical(res_treatment$padj < p.adj.cutoff)
res_origin$threshold <- as.logical(res_origin$padj < p.adj.cutoff)
# Sort the results tables
res_treatment_sorted <- res_treatment[order(res_treatment$padj), ]
res_origin_sorted <- res_origin[order(res_origin$padj), ]
# Get significant genes
sigTrmt <- row.names(res_treatment_sorted)[which(res_treatment_sorted$threshold)]
sigOrigin <- row.names(res_origin_sorted)[which(res_origin_sorted$threshold)]
sigTrmtOnly <- sigTrmt[!sigTrmt %in% sigOrigin]
sigOriginOnly <- sigOrigin[!sigOrigin %in% sigTrmt]
sigTrmtOrigin <- sigTrmt[sigTrmt %in% sigOrigin]

candidates=list('Treatment'=sigTrmt,'Origin'=sigOrigin)
venn(candidates)
# not many genes overlapping between origin and treatment -- only 7
vsd1=vst(dds, blind=FALSE)
plotPCA(vsd1,intgroup=c('origin','treatment')) #looks the same as blind = TRUE

# let's take a look at some heatmaps by origin and treatment

# scale values
norm_df=assay(vsd1)
means=apply(norm_df,1,mean) # means of rows
norm_df_scaled=norm_df-means #rescale expression data so it's up and down
df=as.data.frame(colData(dds)[,c("treatment","origin")])
display.brewer.pal(n = 6, name = 'Set2')
brewer.pal(n = 6, name = "Set2")
ann_colors=list(
  treatment=c(Control='#66C2A5', Heat='#FC8D62'),
  origin=c(off='#A6D854','in'='#FFD92F')
)

# treatment DEGs
col=color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.4)(100)
pheatmap(norm_df_scaled[sigTrmt,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment then stage for all sig trmt DEGs

# origin DEGs
col=color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")),bias=1.3)(100)
pheatmap(norm_df_scaled[sigOrigin,],color=col,cluster_cols=T,clustering_distance_rows="correlation",show_rownames=FALSE,annotation_col=df,annotation_colors = ann_colors) #clusters mainly by treatment then stage for all sig trmt DEGs

# both cluster nicely automatically by origin / treatment :)

out=as.data.frame(cbind(norm_df,res_treatment$padj,res_origin$padj,res_treatment$log2FoldChange,res_origin$log2FoldChange))
names(out)[15:18]=c('trmt_padj','stage_padj','trmt_LFC','origin_LFC')

write.csv(out,file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/Larv_DESeq2_VsdPvalsLFC_19July2021.csv')


############## now make origin-treatment groups and test to see if inshore and offshore larvae are using the same genes to respond to treatment #####
# combine the factors of interest into a single factor with all combinations of the original factors
dds$group <- factor(paste0(dds$origin, dds$treatment))
# change the design to include just this factor, e.g. ~ group
design(dds) <- ~ group
dds <- DESeq(dds)
res_group=results(dds)
summary(res_group) #90 up, 110 down
resultsNames(dds)

res_in <- results(dds, contrast=c("group", "inHeat", "inControl")) # 9 up, 16 down
res_off <- results(dds, contrast=c("group", "offHeat", "offControl")) # 68 up, 62 down
res_origin_control <- results(dds, contrast=c('group','offControl','inControl')) # 18 up,39 down
res_origin_heat <- results(dds, contrast=c('group','offHeat','inHeat')) # 13 up, 14 down

in_trmtgenes = rownames(subset(res_in,padj < 0.1))
off_trmtgenes = rownames(subset(res_off, padj < 0.1))

candidates=list('inshore response' = in_trmtgenes, 'offshore response' = off_trmtgenes)
venn(candidates) # only 11 overlapping between inshore and offshore response

# how many overlap with overall trmt genes?
candidates=list('inshore response' = in_trmtgenes, 'offshore response' = off_trmtgenes, 'all response genes' = sigTrmt)
venn(candidates) # 50 offshore response genes not detected in full dataset

# interesting, a lot more genes responding to treatment in offshore (130 DEGs) compared to inshore (25 DEGs)
# that is opposite of what we saw in adults and recruits after 2.5 week heat stress
# maybe offshore have an initial response and then can't sustain it?

#save results
vsd2=vst(dds,blind=F)
norm_df2=assay(vsd2)

out=as.data.frame(cbind(norm_df2,res_in$padj,res_off$padj,res_origin_control$padj,res_origin_heat$padj,res_in$log2FoldChange,res_off$log2FoldChange,res_origin_control$log2FoldChange,res_origin_heat$log2FoldChange))
names(out)[15:22]=c('inTrmt_padj','offTrmt_padj', 'originControl_padj','originHeat_padj','inTrmt_LFC','offTrmt_LFC','originControl_LFC','originHeat_LFC')

write.csv(out,file='downstream_analyses/Mar2020/DESeq/DESeq2_June2021/LarvDESeq2_VsdPvalsLFC_byGroup_19July2021.csv')

####### look if any interaction genes ########
design(dds) <- ~ origin + treatment + origin:treatment
dds <- DESeq(dds)
resultsNames(dds)
res_int=results(dds,name = 'originoff.treatmentHeat')
summary(res_int) # no genes significant for interaction

# see how many treatment genes in larvae overlap with treatment genes in adults/recruits
AR_df=read.csv('DESeq2_vsd_pvals_30June2021.csv',row.names = 1)
AR_trmtDEGs=rownames(subset(AR_df,trmt_padj < 0.05))

# response to treatment
candidates = list('Adult and Recruit response' = AR_trmtDEGs, 'Larval response' = sigTrmt)
venn(candidates) # 54 overlapping

# origin DEGs
AR_originDEGs=rownames(subset(AR_df, origin_padj < 0.1))
candidates = list('Adult and Recruit origin DEGs' = AR_originDEGs, 'Larval origin DEGs' = sigOrigin)
venn(candidates) # 33 overlapping

###### 




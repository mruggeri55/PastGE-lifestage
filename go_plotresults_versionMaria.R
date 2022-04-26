# this is the final script in GO analysis using MWU test, as in
# Voolstra et al PLoS ONE 2011, 6(5): e20392.

setwd('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/downstream_analyses/Mar2020/DESeq/DESeq2_June2021/')

#Note: must run this script while results of GO_MWU are loaded into R
#must run these to initiate, then can ignore...
#gg=read.table("~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/past_blast_annot_July2021/past_iso2gene.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)  #the iso2gene file for your species
gg=read.table("~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/kb8_blast_annot_July2021/kb8_iso2gene.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)  #the iso2gene file for your species

#dfull=read.csv("DESeq2_vsd_pvals_30June2021.csv",row.names=1) # for host
dfull=read.csv('DESeq2_SYM_vsd_pvals_July2021.csv',row.names = 1) # for sym

efulldata=c(1:44) #the columns where your expression data are 1:47 host; 1:44 sym
#source("summarySE.R")
library(pheatmap)

###############################Recording significances of top GO terms - rename file below
# sig=goods
# sig$p=(10^(-(goods$pval)))
# sig
# 
# write.csv(sig,file=paste("GOsig",division,fileglob,sep="_"),quote=F,row.names=T)
# #************Note: to reconvert these scores into pvalues, you must recalculate: pval=10^(-(score))



###############plotting what significant GO genes are

#must change these each time
golist=read.table("GO_MWU-master/GO_output_files/BP_GOrankSYM_ARtrmt.csv",sep="	",header=T) #read in proper GO list from gomwu output - BP/MF/CC and module name
in.mwu=read.table("GO_MWU-master/GO_output_files/MWU_BP_GOrankSYM_ARtrmt.csv",sep=" ",header=T)
#d=read.csv("VSDs_GObinarySymPC1_NoOutlierGEO.csv") #read in proper VSD file for genes in module - generated in wgcna script
# inMOD=read.csv("topGO_Eonly_Sig_HighExpOrpheus.csv")
# inMODt=inMOD[inMOD$moduleColor>0,]
# rownames(inMODt)<-inMODt$gene
# d=as.data.frame(dfull[rownames(dfull) %in% rownames(inMODt),]) #OR, simply genes where pval for factor of interest is below specific expression level (same as your absValue cutoff)
d5=dfull[dfull$trmt_padj<=0.1 & !is.na(dfull$trmt_padj),]
edata=c(1:44) #only columns with your expression data

#in.mwu=paste("MWU",goDivision,input,sep="_")
#pv=read.table(in.mwu,header=T)

gene=subset(in.mwu,name=="photosynthesis") #write GO term of interest here from your sig list#
t=gene$term
t # Go terms with this annotation

is=as.character(golist$seq[grep(t,golist$term,ignore.case=T)])
length(is) #then pulls out the isogroups with this annotation

#####################first loop through genes matching GO term in module or in "good expression" subset
sel=c();gnms=c()
for ( i in is){
	if (i %in% rownames(d5)){
		sel=rbind(sel,d5[rownames(d5)==i,])
		gnms=append(gnms,substr(paste(as.character(gg$V2[gg$V1==i]),collapse="."),1,50))
	}
}

sel_isos=rownames(sel)
row.names(sel)=paste(gnms,rownames(sel),sep=".")
nrow(sel)
rownames(sel)
# [1] "Photosystem II reaction center protein L OS=Mesost.kb8_rep_c7376" 
# [2] "Photosystem II reaction center protein L OS=Emilia.kb8_rep_c34344"
# [3] "Photosystem I P700 chlorophyll a apoprotein A1 OS=.kb8_rep_c112"  
# [4] "Cytochrome b6-f complex iron-sulfur subunit OS=Syn.kb8_rep_c28708"

exp=sel[,edata]
if (length(exp[,1])==1) { exp=rbind(exp,exp) }
nrow(exp)
#should match 'good genes' tabulation in figure
rownames(exp) ## yayy this is working
#rownames(daat)

#daat<-exp#[1,]
#daat=rbind(daat,exp)

#write.csv(daat,file="Sym_Ori_MF_GOgenes.csv",quote=FALSE)

#------------------ heatmap of GOgenes stacked by symbiont density (or other trait...)
library(pheatmap)

#must read in traits data once then comment out
traits=read.csv("~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/downstream_analyses/CG_Larv_WGCNA_traits.csv")
rownames(traits)<-traits$colony #make sample names the rownames
traits=traits[traits$adult==1,]
traits=traits[,c(2,7,8,9,11)]
traits=traits[order(traits$bleachstatus),] #change trait here to reorder by different factor: "origin"   "treatment" "bleachstatus"

############# For DESeq2 analyses
# subset vsd and scale
norm_df=dfull[,1:44]
means=apply(norm_df,1,mean) # means of rows
norm_df_scaled=norm_df-means #rescale expression data so it's up and down
colnames(norm_df_scaled)=gsub('X','',colnames(norm_df_scaled))
colnames(norm_df_scaled)=gsub('CL','RC',colnames(norm_df_scaled))
colnames(norm_df_scaled)=gsub('HL','RH',colnames(norm_df_scaled))

# setting up metadata table
colData=read.csv('CG_colData.csv',row.names = 1)
colData$stage_treatment=paste(colData$stage, colData$treatment, sep='_')
rownames(traits)=gsub('X','',rownames(traits))
bleach_order=traits[match(rownames(colData),rownames(traits)),]
colData$bleach_score=bleach_order$bleachstatus

# subset just traits for annotation
df=subset(colData,select=c('bleach_score','stage_treatment'))
rownames(df)=gsub('CL','RC',rownames(df))
rownames(df)=gsub('HL','RH',rownames(df))

norm_df_scaled<-norm_df_scaled[,c(1,5,8,12,16,20,23,27,31,34,38,42,2,9,13,17,24,28,32,35,39,3,6,10,14,18,21,25,29,33,36,40,43,4,7,11,15,19,22,26,30,37,41,44)]
norm_df_scaled<-norm_df_scaled[,c(1,2,3,4,10,11,12,5,6,7,8,9,13,14,15,20,21,16,17,18,19,22,23,24,25,31,32,33,26,27,28,29,30,34,35,36,37,42,43,44,38,39,40,41)]
norm_df_scaled<-norm_df_scaled[,c(1,3,4,5,6,7,9,10,11,12,8,2,17,14,16,18,19,13,15,20,21,22:44)]

library(RColorBrewer)
bleach=brewer.pal(n = 6, name = "YlOrBr")

ann_colors=list(
  bleach_score=bleach,
  stage_treatment=c(Adult_Control='#762A83', Adult_Heat='#AF8DC3',Recruit_Control="#1B7837",Recruit_Heat="#7FBF7B")
  #origin=c(off='aquamarine3','in'='orange'),
)

library(RColorBrewer)
col=color= colorRampPalette(rev(c('yellow','yellow','black','turquoise','turquoise')),bias=0.6)(30)
# 4 sig photo genes
pheatmap(norm_df_scaled[sel_isos,],color=col,cluster_cols=F,cluster_rows = T,show_rownames=T,annotation_col=df,annotation_colors = ann_colors, gaps_col = c(12,21,33)) 
# now just kb8_rep_c7376 (psbL)
col=color= colorRampPalette(rev(c('yellow','yellow','black','black','aquamarine4','turquoise')),bias=0.45)(50)
pheatmap(norm_df_scaled[sel_isos[1],],color=col,cluster_cols=F,cluster_rows = F,show_rownames=F,annotation_col=df,annotation_colors = ann_colors, gaps_col = c(12,21,33)) 


# all photo genes
pheatmap(norm_df_scaled[is,],color=col,cluster_cols=F,cluster_rows = T,show_rownames=F,annotation_col=df,annotation_colors = ann_colors, gaps_col = c(12,21,33)) 


# look at the fold changes of these genes
sel
# 3/4 downregulated under heat stress, one upregulated
sel_isos
#"kb8_rep_c7376"  "kb8_rep_c34344" "kb8_rep_c112"   "kb8_rep_c28708"

# let's look at just adult FC
load(file = 'DESeq2_sym_by_stage.RData')
vsd2=vst(dds,blind = F)
norm_df=assay(vsd2)
means=apply(norm_df,1,mean) # means of rows
norm_df_scaled=norm_df-means #rescale expression data so it's up and down
norm_sub=norm_df_scaled[rownames(norm_df_scaled) %in% sel_isos,]
norm_sub=t(norm_sub)
rownames(norm_sub)=gsub('X','',rownames(norm_sub))
norm_merge=merge(norm_sub,traits,by='row.names')
norm_merge
treatment=c(1:length(norm_merge))
treatment[grep('C',norm_merge$Row.names)]='Control'
treatment[grep('H',norm_merge$Row.names)]='Heat'
norm_merge$treatment=treatment

ggplot(norm_merge,aes(x=kb8_rep_c112,y=bleachstatus))+geom_point()+geom_smooth(method='lm')
ggplot(norm_merge,aes(x=kb8_rep_c28708,y=bleachstatus))+geom_point()+geom_smooth(method='lm')
ggplot(norm_merge,aes(x=kb8_rep_c34344,y=bleachstatus))+geom_point()+geom_smooth(method='lm')
ggplot(norm_merge,aes(x=kb8_rep_c7376,y=bleachstatus))+geom_point(aes(color=treatment))+geom_smooth(method='lm')

bleach=brewer.pal(n = 6, name = "YlOrBr")
ggplot(na.omit(norm_merge),aes(x=bleachstatus,y=kb8_rep_c7376))+
  geom_point(aes(color=bleachstatus),size=2)+
  geom_smooth(method='lm',color='black',se=F)+
  ggtitle("Photosystem II reaction center protein")+
  ylab('expression')+
  xlab('bleaching score')+
  annotate('text',x=3,y=0,label='y = 0.98x - 5.4')+
  annotate('text',x=3,y=-0.5,label='R2 = 0.572')+
  annotate('text',x=3,y=-1,label='p = 4e-5')+
  scale_color_gradient(low="#FFFFD4",high="#993404")+
  theme_bw(base_size=15)


model=lm(kb8_rep_c7376~bleachstatus,data=norm_merge)
summary(model) 

# kb8_rep_c7376 significant p = 4.4e-5, y=0.98x-5.4, R2=0.572
# other genes not significant

library(reshape2)
melt=melt(norm_merge[,c(1:5,8,9)],id.vars = c('bleachstatus','Row.names','treatment'))
colnames(melt)[4:5]=c('gene','expression')
melt=na.omit(melt)
ggplot(melt,aes(y=expression,x=bleachstatus))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm',se=F,color='black')+
  facet_wrap(~gene)+
  theme_bw(base_size=15)


################# Carly's below this line
gene=as.data.frame(t(sorted))
par(mfrow=c(3,6))

SiteHost=c(rep("DB",times=14),rep("DC",times=14),rep("IB",times=15),rep("IC",times=15))
SiteSeep=c(rep("Seep",times=14),rep("Ctrl",times=14),rep("Seep",times=15),rep("Ctrl",times=15))

for(name in names(gene)){
boxplot(gene[,name]~SiteHost,xlab=paste(name))
#SiteSeep=c(rep("Seep",times=14),rep("Ctrl",times=14),rep("Seep",times=15),rep("Ctrl",times=15))
print(name)
print(anova(lm(gene[,name]~SiteSeep)))}

abline(lm(traits$MicrobePCo1[1:19]~gene[,6]),lty=2)
anova(lm(gene[,4]~SiteSeep))

#traits$MicrobePCo1[1:19]

####plotting individual gene-factor correlations

origin=transplant=c(1:length(names(sorted)))
origin[grep("KK",names(sorted))]="Keppels"
origin[grep("KO",names(sorted))]="Keppels"
origin[grep("OK",names(sorted))]="Orpheus"
origin[grep("OO",names(sorted))]="Orpheus"

transplant[grep("KK",names(sorted))]="Keppels"
transplant[grep("OK",names(sorted))]="Keppels"
transplant[grep("OO",names(sorted))]="Orpheus"
transplant[grep("KO",names(sorted))]="Orpheus"

survival=c(rep("Dead",12),rep("Live",42))

#Use Eigengene for a desired cluster
nrow(sorted)
moduleColors=rep("blue",26) #change rep to number of genes
MEcluster = moduleEigengenes(t(sorted[,c(1:42)]), moduleColors)$eigengenes
final=cbind(MEcluster,traits$GAIN[1:42])
par(mfrow=c(1,2))
plot(traits$GAIN[1:42]~MEblue,dat=final,main="Gain ~ MEsYellow")
abline(lm(traits$GAIN[1:42]~MEblue,dat=final))
summary(lm(traits$GAIN[1:42]~MEblue,dat=final))


texp=data.frame(t(sorted))
sexp=stack(texp)
#unique(sexp$ind) #unique returns a list of all the gene names in your GOterm for the module of choice


gene=subset(sexp,ind=="Caspase..apoptotic.cysteine.protease...caspase.inv.isogroup5161") #insert name of interest from unique output
final=cbind(gene,survival)

quartz()
par(mfrow=c(1,3))
plot(values~survival,dat=final,main="Caspase_isogroup5161") #can plot regression of trait value vs gene of interest




plot(gene$values~traits$OriByTrans,main="Trypsin_isogroup24511") #can plot regression of trait value vs gene of interest
plot(gene$values~traits$origin,main="Trypsin_isogroup24511")
plot(gene$values~traits$transplant,main="Trypsin_isogroup24511")
plot(gene$values~traits$Survival,main="Trypsin_isogroup24511")
plot(traits$GAIN~traits$Survival,main="Growth:Survival Trade-Off?")

quartz()
barplot(o$Sym_cell.cm2, main="", cex.main=2,
        ylab="Symbiont Density",xlab="sample")

#------------------- ggplot BY genename
# Stacking; setting treatment conditions
# this chunk must be edited, depends on your factors! hope it is not too unclear
library(ggplot2)
texp=data.frame(t(exp))
sexp=stack(texp)
head(sexp,50)
names(sexp)=c("expr","gene")
co=strsplit(row.names(texp),"")
ori=c();tr=c()
for(i in 1:length(co)) {
	ori=append(ori,co[[i]][1])
	tr=append(tr,co[[i]][2])
}	
sexp$ori=ori
sexp$tr=tr
str(sexp)
head(sexp)

#----------------------
summ=summarySE(data=sexp,measurevar="expr",groupvars=c("gene","ori","tr"))

pd=position_dodge(.3)
#summ$tr=factor(summ$tr, levels=c("o","i"))
ggplot(summ,aes(x=tr,y=expr,colour=ori,group=ori))+
	geom_errorbar(aes(ymin=expr-se,ymax=expr+se),lwd=0.4,width=0.3,position=pd)+
	geom_line(aes(group=ori,linetype=ori),position=pd)+
	geom_point(aes(group=ori,pch=ori),position=pd,size=2.5)+
	scale_colour_manual(values=c("coral","cyan3"))+
	facet_wrap(~gene,scales="free_y")+theme_bw()


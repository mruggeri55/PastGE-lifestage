library(flashClust)
library(gplots) 
library(ggplot2)
library(ggbiplot)
library(RColorBrewer)
#library(affycoretools)
library(genefilter)
library(plotrix)
library(reshape2)
library(vegan)
library(adegenet)

setwd('/Users/Maria/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/downstream_analyses/Mar2020/DESeq/DESeq2_June2021/')

#read in variance stabilized count data from DESeq2 analysis
#dat=read.csv('DESeq2HostAR_origin_response_to_treat_VsdPvalsLFC_July2021.csv',row.names=1) #for host
#dat=read.csv('DESeq2HostAR_origin_response_to_treat_VsdPvalsLFC_Oct2021.csv',row.names=1) # controlling for stage
#dat=read.csv('DESeq2_vsd_pvals_30June2021.csv',row.names = 1) # ~ all treatment genes
#dat=read.csv('Larv_DESeq2_VsdPvalsLFC_19July2021.csv',row.names = 1)
#dat=read.csv('LarvDESeq2_VsdPvalsLFC_byGroup_19July2021.csv',row.names = 1)
dat=read.csv('DESeq2_SYM_vsd_pvals_July2021.csv',row.names = 1)


head(dat)
colnames(dat)=gsub('X','',colnames(dat))

######subset into dfs for sig trmt genes and remove p vals
#for inshore response genes
#off_trmt=rownames(subset(dat,offTrmt_padj<0.1))
#in_trmt=rownames(subset(dat,inTrmt_padj<0.1))
#in_only=in_trmt[!in_trmt %in% off_trmt]
trmt=rownames(subset(dat,trmt_padj<0.1))

#dat_in_trmt=subset(dat, rownames(dat) %in% in_only)
#dat_in_trmt=subset(dat, rownames(dat) %in% c(in_only,off_trmt)) # to do all genes
dat_in_trmt=subset(dat,rownames(dat) %in% trmt) # all treatment genes from ~ stage + origin + treatment
#dat_in_trmt=dat

#dat_in_trmt=dat_in_trmt[,1:47] # for AR host
#dat_in_trmt=dat_in_trmt[,1:14] # for larvae
dat_in_trmt=dat_in_trmt[,1:44] # for sym

length(rownames(dat_in_trmt)) #167 inshore responsive genes / 316 when controlling for stage, 291 for in only, 466 for sym all trmt
degs10=rownames(dat_in_trmt)

#read in metadata
meta=read.csv('CG_coldata.csv',row.names = 1)

# for larvae set up metadata
# sample=colnames(dat_in_trmt)
# genotype=gsub("X|A|L|C|H","",sample)
# origin=c('in','in','off','off','off','off','off','off','off','off','in','in','in','in')
# treatment=rep(c('C','H'),times=7)
# meta=data.frame(genotype,origin,treatment)
# rownames(meta)=sample

# # to do in control vs off control
# control=c(rownames(meta)[meta$origin=='in' & meta$treatment == 'Control'],rownames(meta)[meta$origin=='off' & meta$treatment == 'Control'])
# heat=c(rownames(meta)[meta$origin=='in' & meta$treatment == 'Heat'],rownames(meta)[meta$origin=='off' & meta$treatment == 'Heat'])
# 
# a.vsd<-dat[,colnames(dat) %in% control]
# a.vsd.supp<-dat[,colnames(dat) %in% heat]  

# to do inshore vs offshore
o_in=rownames(meta)[meta$origin=='in']
o_off=rownames(meta)[meta$origin=='off']
# #
a.vsd<-dat_in_trmt[,colnames(dat_in_trmt) %in% o_in] # defines functions
a.vsd.supp<-dat_in_trmt[,colnames(dat_in_trmt) %in% o_off] 

# # to subset out adults
# adults_in=rownames(meta)[meta$origin =='in' & meta$stage == 'Adult']
# adults_off=rownames(meta)[meta$origin =='off' & meta$stage == 'Adult']
# #
# a.vsd<-dat_in_trmt[,colnames(dat_in_trmt) %in% adults_in]
# a.vsd.supp<-dat_in_trmt[,colnames(dat_in_trmt) %in% adults_off]

# # to subset out recruits
# recruits_in=rownames(meta)[meta$origin =='in' & meta$stage == 'Recruit']
# recruits_off=rownames(meta)[meta$origin =='off' & meta$stage == 'Recruit']
# #
# a.vsd<-dat_in_trmt[,colnames(dat_in_trmt) %in% recruits_in]
# a.vsd.supp<-dat_in_trmt[,colnames(dat_in_trmt) %in% recruits_off]

###################### Some genes have insufficient variance for DFA in data subsets!...find those genes using loop below
dframe=(a.vsd[degs10,])
for(col in rownames(dframe)){  
  min=min(dframe[col,])
  max=max(dframe[col,])
  if(min == max){print(col)}
}

#degs10=degs10[! degs10 %in% c("isogroupTR104649_c0_g1")] # for adult subset inshore response genes

###################### moving on to the analysis
pcp=prcomp(t(a.vsd[degs10,]), retx=TRUE, center=TRUE, scale.=TRUE) #scale.=TRUE
scores=pcp$x
screeplot(pcp,bstick=T) 

# adegenet: finding clusters (even though we know what clusters we want) - keep all PCs so 30 and 2 clusters
clus=find.clusters(t(a.vsd[degs10,]),max.n.clus=15)
# keep all so 30 PCs and 2 clusters

#clus=find.clusters(t(a.vsd[degs10,]),max.n.clus=13) # for adult and recruit subset
#clus=find.clusters(t(a.vsd[degs10,]),max.n.clus=7) # for larvae

#Use clus$grp to rename to heat and control
clus$grp 
##### tell the DF which groups you want to cluster; 1=control, 2=heat
#clus$grp=rep(c(1,2),times=14) #for inshore
#clus$grp=c(rep(c(1,2),times=9),1) # for offshore
#clus$grp=rep(c(1,2),times=7) # for adult inshore or recruit inshore
#clus$grp=c(rep(c(1,2),times=8),rep(c(3,4),times=9),3,rep(c(1,2),times=6)) # by origin
#clus$grp=c(rep(c(1,2,3,4),times=8),1,2,3,rep(c(1,2,3,4),times=3)) # by stage
#clus$grp=c(rep(1,times=8),rep(2,times=10),rep(1,times=6)) # for control, 1 = inshore, 2 = offshore
#clus$grp=c(rep(c(1,2),times=4)) # larvae offshore
#clus$grp=c(rep(c(1,2),times=2),1,1,2,rep(c(1,2),times=8),1,1,2) # for sym def inshore
clus$grp=c(1,2,1,2,1,1,2,1,2,1,2,1,2,1,2,1,2,1) # for sym def by offshore

# now lets build a discriminant function for these two groups:
dp=dapc(t(a.vsd[degs10,]),clus$grp)
# chose 10 PCs and 1 function dor AR, 3 PCs for larvae
#scatter(dp)
#scatter(dp,2,2, bg="white",scree.da=FALSE, legend=TRUE, solid=.4)


########## Now, add in offshore data and see where they fall along continuum
pred.sup<-predict.dapc(dp,newdata=(t(a.vsd.supp[degs10,])))
names(pred.sup)
names(a.vsd.supp)
pred.sup$assign
#pred.sup$assign<-c(rep(c(1,2),times=9),1) # for host off
#pred.sup$assign<-rep(c(1,2),times=5) # for host adult offshore
pred.sup$assign<-c(rep(c(1,2),times=2),1,1,2,rep(c(1,2),times=8),1,1,2) # for sym inshore


#must create another dataframe structure in order to plot these predicted values
test<-dp
test$ind.coord<-pred.sup$ind.scores
test$posterior<-pred.sup$posterior
test$assign<-pred.sup$assign

test$grp=as.factor(pred.sup$assign)

##retain DFA values for additional calculations
dpc=data.frame(rbind(dp$ind.coord,pred.sup$ind.scores))

#write.csv(dpc,"DAPC_Host_trmt_by_origin.csv",quote=F)
#write.csv(dpc,"DAPC_LarvaeHost_allTreatGenes283_byOff.csv",quote=F)
#write.csv(dpc,"DAPC_SYM_AllTreatGenes_byOff.csv",quote=F)
dpc=read.csv("DAPC_Host_trmt_by_origin.csv",row.names = 1)

# ###Testing significance of DFA differences - MCMCglmm

# setting up genotype, trmt and stage factors
#meta=read.csv('CG_colData.csv',row.names = 1)
order=rownames(dpc)
meta_order=meta[match(order,rownames(meta)),]
rownames(meta_order)==rownames(dpc)
coral=as.data.frame(cbind(dpc,meta_order))

#plotting on same axes -- make new column 'group'
coral$group=paste(coral$origin,coral$treatment)


#quartz()
ggplot(coral,aes(x=LD1,fill=group,color=group))+geom_density(size=1)+theme_classic()+
  scale_fill_manual(values=c(alpha('orange',0.8),alpha('orange',0.3),alpha('cyan4',0.8),alpha('cyan4',0.3)))+
  scale_color_manual(values=c('orange','orange','cyan4','cyan4'))+
  xlab('Discriminant function')+ylab('sample density')+
  theme(axis.text = element_text(size=16),axis.title=element_text(size=20),legend.text = element_text(size=16),
        legend.title = element_blank(),legend.position = 'bottom')+
  annotate('segment',x=-2.8,xend=2.5,y=0.75,yend=0.75,color='orange',size=2,arrow=arrow())+
  annotate('segment',x=-3.8,xend=4.2,y=0.7,yend=0.7,color='cyan4',size=2,arrow=arrow())

# edit segments to match peaks

library(MCMCglmm)

# weak inverse wishart prior with parameter expansion for random effect of genotype (this is standard in MCMCglmm, the results are actually identical with the default uniform prior)
prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))
cd=MCMCglmm(LD1~origin+treatment:origin,random=~genotype,data=coral,prior=prior,nitt=75000, thin=25, burnin=5000)
summary(cd)
#                         post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)               -3.1503  -3.7728  -2.5544     2800 <4e-04 ***
#   originoff                  1.2475   0.2404   2.1848     2800 0.0157 *  
#   originin:treatmentHeat     6.2983   5.5148   7.0814     3127 <4e-04 ***
#   originoff:treatmentHeat    3.6954   2.6770   4.5871     2800 <4e-04 ***

#control for stage
# post.mean   l-95% CI   u-95% CI eff.samp  pMCMC    
# (Intercept)             -2.6337283 -3.2944524 -2.0265495     2800 <4e-04 ***
#   originoff                0.9582221  0.0001154  1.8943862     2800 0.0443 *  
#   originin:treatmentHeat   5.2643377  4.5678452  5.9446839     2800 <4e-04 ***
#   originoff:treatmentHeat  2.7671162  1.8489735  3.5280381     2800 <4e-04 ***

# for larvae
# post.mean  l-95% CI  u-95% CI eff.samp  pMCMC    
# (Intercept)           -8.98404 -10.56943  -7.45701     2800 <4e-04 ***
# originoff              0.02555  -1.96087   2.12765     2800  0.996    
# originin:treatmentH   14.74306  12.91471  16.58236     2754 <4e-04 ***
# originoff:treatmentH  17.95675  16.35769  19.56309     2800 <4e-04 ***

# calculating difference in magnitudes of adult and recruit response to trmt using sampled sets of parameters:
origin_Delta=abs(cd$Sol[,"originoff:treatmentHeat"])-abs(cd$Sol[,"originin:treatmentHeat"])

# 95% credible interval:
HPDinterval(origin_Delta) 
# lower    upper
# var1 1.325224 3.779484

#MCMC p-value:
if (is.na(table(origin_Delta<0)[2])) {
  cat("p <",signif(1/length(origin_Delta),1))
} else { cat("p =",signif(table(origin_Delta<0)[2]/length(origin_Delta),2)) }
##### MCMC results #####
# p < 4e-04 aka inshore have a significantly larger response for these genes compared to offshore

# significant just for adult subset as well for inshore responsive genes
# when using combined inshore and offshore response, inshore greater plasticity p <4e-04

# when using all 605 trmt genes p=0.0025
# significant for adult subset all treatment genes (605) p = 0.042
# significant for recruit subset all treatment genes (605) p = 0.025

# no different in larvae for all treatment genes, but did define on inshore response
# p = 0.0089 defining by offshore response
# p = 0.00071 for c(inshore response, offshore response) defining based on offshore
# NS for c(inshore response, offshore response) defining based on inshore -- note still clustered by treatment

# significantly larger response (466 trmt genes) in inshore syms (p=0.00071) when defined by inshore
# BUT also sig (p = 0.0018) when define by offshore
# so can't really trust sym DAPC


####### does plasticity correlate to diffs in bleach score? for adults
bleach_dat=read.csv('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/downstream_analyses/CG_Larv_WGCNA_traits.csv')
bleach_dat$colony=gsub('X','',bleach_dat$colony)
bleach_dat=bleach_dat[bleach_dat$adult==1,]
#bleach_dat=bleach_dat[bleach_dat$larvae==1,]

#subset just adults from DF1
adult_DF=coral[coral$stage=='Adult',]
#adult_DF=coral
adult_DF$colony=rownames(adult_DF)
# now merge together
merge_df=merge(adult_DF,bleach_dat,by='colony')

plot(merge_df$LD1~merge_df$bleachstatus)
abline(lm(merge_df$LD1~merge_df$bleachstatus))

ggplot(merge_df,aes(x=LD1,y=bleachstatus, color=origin.x))+
  geom_point(size=2)+geom_line(aes(group=as.factor(genotype)))
  #geom_smooth(method='lm',se=F)
  #scale_color_manual()

#ggplot(merge_df,aes(x=LD1,y=chlorophyll, color=origin.x))+
#  geom_point(size=2)+geom_line(aes(group=as.factor(genotype)))

# subset control and heat to get the difference aka plasticity
control=subset(merge_df,subset = merge_df$treatment.x=='Control',select = c(colony,LD1,bleachstatus,genotype,origin.x))
heat=subset(merge_df,subset = merge_df$treatment.x=='Heat',select = c(colony,LD1,bleachstatus,genotype,origin.x))

# merge by genotype
merge=merge(control,heat,by='genotype')
# calculate magnitude difference along axis
merge$del_LD1=abs(merge$LD1.x)+abs(merge$LD1.y)
merge$del_bleachscore=merge$bleachstatus.x-merge$bleachstatus.y
#plastic=as.data.frame(cbind('genotype'=control$genotype,'origin'=as.character(control$origin.x),'del_LD1'=del_LD1,'del_bleachscore'=del_bleachscore))
#plastic=as.data.frame(cbind(del_LD1,del_bleachscore))
plastic=subset(merge,select=c('del_LD1','del_bleachscore'))

adult=ggplot(plastic,aes(x=del_LD1,y=del_bleachscore))+geom_point(size=2)+geom_smooth(method='lm',se=F)+ylab('change in bleach score')+xlab('change in expression profile')+ggtitle('Adult')
ggplot(plastic,aes(x=del_bleachscore,y=del_LD1, color=del_bleachscore))+geom_point(size=2)+geom_smooth(method='lm',se=F)

lm=lm(plastic$del_LD1~plastic$del_bleachscore)
summary(lm) # no relationship

ggplot(merge_df,aes(x=treatment.x,y=bleachstatus,color=origin.x,group=as.factor(genotype)))+geom_point(position = position_dodge(0.3),size=2)+geom_line(position = position_dodge(0.3),size=1.2)

library(dplyr)
# sum bleach score and plot out by treatment
sum=merge_df %>% group_by(origin.x,treatment.x) %>% summarise(mean=mean(bleachstatus),sd=sd(bleachstatus))
# origin.x treatment.x  mean    sd
# <fct>    <fct>       <dbl> <dbl>
# 1 in       Control      5.86 0.378
# 2 in       Heat         4.43 1.13 
# 3 off      Control      5.8  0.447
# 4 off      Heat         3.5  1.5
ggplot(sum,aes(x=treatment.x,y=mean,color=origin.x,group=origin.x))+
  geom_point(position = position_dodge(0.3),size=3)+
  geom_line(position = position_dodge(0.3),size=1.2)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),position = position_dodge(0.3),width=0.5,size=1.2)+
  ylab('bleaching score')+xlab('treatment')+theme_classic(base_size = 20)+
  scale_color_manual(values=c('orange','cyan4'))+labs(color='origin')

######### for larvae
coral=read.csv("DAPC_LarvaeHost_allTreatGenes283_byOff.csv", row.names=1)
#set up metadata
sample=rownames(coral)
genotype=gsub("X|A|L|C|H","",sample)
origin=c('in','in','off','off','off','off','off','off','off','off','in','in','in','in')
treatment=rep(c('C','H'),times=7)
meta=data.frame(genotype,origin,treatment)
rownames(meta)=sample
coral=merge(coral,meta,by='row.names')
rownames(coral)=coral$Row.names
coral$colony=coral$Row.names
coral$Row.names=NULL

bleach_dat=read.csv('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/downstream_analyses/CG_Larv_WGCNA_traits.csv')
bleach_dat$colony=gsub('X','',bleach_dat$colony)
bleach_dat=bleach_dat[bleach_dat$larvae==1,]

# now merge together
merge_df=merge(coral,bleach_dat,by='colony')

plot(merge_df$LD1~merge_df$chlorophyll)
abline(lm(merge_df$LD1~merge_df$chlorophyll))

control=subset(merge_df,subset = merge_df$treatment.x=='C',select = c(colony,LD1,chlorophyll,genotype,origin.x))
heat=subset(merge_df,subset = merge_df$treatment.x=='H',select = c(colony,LD1,chlorophyll,genotype,origin.x))

del_LD1=abs(control$LD1)+abs(heat$LD1)
del_chlorophyll=control$chlorophyll-heat$chlorophyll
#plastic=as.data.frame(cbind('genotype'=control$genotype,'origin'=as.character(control$origin.x),'del_LD1'=del_LD1,'del_bleachscore'=del_bleachscore))
plastic=as.data.frame(cbind(del_LD1,del_chlorophyll))

larvae=ggplot(plastic,aes(x=del_LD1,y=del_chlorophyll))+geom_point(size=2)+geom_smooth(method='lm',se=F)+ylab('change in chlorophyll')+xlab('change in expression profile')+ggtitle('B. Larvae')
ggplot(plastic,aes(x=del_chlorophyll,y=del_LD1, color=del_chlorophyll))+geom_point(size=2)+geom_smooth(method='lm',se=F)

library(gridExtra)
grid.arrange(adult,larvae,nrow=1)

lm=lm(plastic$del_LD1~plastic$del_chlorophyll)
summary(lm) # NS relationship

ggplot(merge_df,aes(x=treatment.x,y=chlorophyll,color=origin.x,group=as.factor(genotype)))+geom_point(position = position_dodge(0.3),size=2)+geom_line(position = position_dodge(0.3),size=1.2)
sum=merge_df %>% group_by(origin.x,treatment.x) %>% summarise(mean=mean(chlorophyll),sd=sd(chlorophyll))
# origin.x treatment.x    mean      sd
# <fct>    <fct>         <dbl>   <dbl>
#   1 in       C           0.00961 0.00233
# 2 in       H           0.00885 0.00120
# 3 off      C           0.00845 0.00446
# 4 off      H           0.00567 0.00335
ggplot(sum,aes(x=treatment.x,y=mean,color=origin.x,group=origin.x))+
  geom_point(position = position_dodge(0.3),size=3)+
  geom_line(position = position_dodge(0.3),size=1.2)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),position = position_dodge(0.3),width=0.5,size=1.2)+
  ylab('chlorophyll')+xlab('treatment')+theme_classic(base_size = 20)+
  scale_color_manual(values=c('orange','cyan4'))+labs(color='origin')



##### make venn diagram of sig genes overlapping between in and off response ####
in_genes=rownames(subset(dat,in_padj<0.1))
off_genes=rownames(subset(dat,off_padj<0.1))

candidates=list('inshore response'=in_trmt,'offshore response'=off_trmt)

library(VennDiagram)
venn.diagram(
  x=candidates,
  category.names = c('inshore','offshore'),
  filename='Larv_TrmtResponse_byOrigin_Venn.png',
  margin=0.15,
  lwd=5,
  col=c("orange","cyan4"),
  fill = c(alpha("orange",0.5), alpha('cyan4',0.5)),
  cex = c(2,2,3),
  cat.col = 'black',
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(-90, 90),
  cat.dist = c(0.12, 0.12)
)

#export table of sig isogroups with gene names
iso2gene=read.table('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/past_blast_annot_July2021/past_iso2gene.tab',sep = '\t', colClasses = "character")
InAndOff=in_genes[in_genes %in% off_genes]
inshore=iso2gene[iso2gene$V1 %in% in_genes & !iso2gene$V1 %in% InAndOff,] # 84 annotated
offshore=iso2gene[iso2gene$V1 %in% off_genes & !iso2gene$V1 %in% InAndOff,] # 5 annotated
InAndOff=iso2gene[iso2gene$V1 %in% InAndOff,] # 8 annotated

all=rbind(inshore,offshore,InAndOff)
all$origin=c(rep('in',times=84),rep('off',times=5),rep('in & off',times=8))

LFCs=subset(dat,subset=rownames(dat) %in% all$V1,select=c('in_LFC','off_LFC'))
all_order=all[match(rownames(LFCs),all$V1),]
merge_df=as.data.frame(cbind(all_order,LFCs))
merge_df$gene=as.character(lapply(strsplit(merge_df$V2, split='OS='),'[',1))
merge_df$GN=as.character(lapply(strsplit(merge_df$V2, split='GN='),'[',2))
merge_df$GN=as.character(lapply(strsplit(merge_df$GN, split='PE='),'[',1))
merge_df$V2=NULL
names(merge_df)=c('isogroup','origin','inshore (LFC)', 'offshore (LFC)', 'gene','GN')
write.csv(x=merge_df,file='GeneNames_trmt_by_origin.csv')

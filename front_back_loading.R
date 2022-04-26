setwd('/Users/Maria/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/downstream_analyses/Mar2020/DESeq/DESeq2_June2021/')
library(VennDiagram)
library(gplots)

### first by stage ####
res_Adult=read.csv("DESeq2_Adult_trmt_response.csv",row.names = 1)
res_Recruit=read.csv("DESeq2_Recruit_trmt_response.csv", row.names = 1)
res_stage_control=read.csv('DESeq2_stage_control_AR_response.csv',row.names = 1) # note this is adult / recruit control LFC
iso2gene=read.table('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/past_blast_annot_July2021/past_iso2gene.tab',sep = '\t', colClasses = "character")

Adult_trmtGenes=rownames(subset(res_Adult,padj < 0.1))
Recruit_trmtGenes=rownames(subset(res_Recruit,padj < 0.1))

# let's check if there are any recruit DEGs being frontloaded / backloaded
# aka look at if any of the treatment responsive genes not responding in recruits are upregulated relative to adult samples in control conditions
AgenesOnlyUP=rownames(subset(res_Adult,padj < 0.1 & log2FoldChange > 0  & !rownames(res_Adult) %in% Recruit_trmtGenes))
RgenesOnlyUP=rownames(subset(res_Recruit,padj < 0.1 & log2FoldChange > 0  & !rownames(res_Recruit) %in% Adult_trmtGenes))

AgenesOnlyDown=rownames(subset(res_Adult,padj < 0.1 & log2FoldChange < 0  & !rownames(res_Adult) %in% Recruit_trmtGenes))
RgenesOnlyDown=rownames(subset(res_Recruit,padj < 0.1 & log2FoldChange < 0  & !rownames(res_Recruit) %in% Adult_trmtGenes))

# get genes upreg or downreg in recruits relative to adults
up_control=subset(res_stage_control,padj < 0.1 & log2FoldChange < 0) #upreg in recruits rel to adults
down_control=subset(res_stage_control,padj < 0.1 & log2FoldChange > 0)

stage_control_up=rownames(up_control)
stage_control_down=rownames(down_control)

# frontloading in recruits
candidates=list('adult heat response up' = AgenesOnlyUP,'recruit vs adult control up' = stage_control_up)
venn(candidates) # 27 genes frontloaded by recruits, but responding in adults
recruit_front=iso2gene[iso2gene$V1 %in% AgenesOnlyUP[AgenesOnlyUP %in% stage_control_up],] # 5 annotated
# [1] "Collagen triple helix repeat-containing protein 1 OS=Homo sapiens OX=9606 GN=CTHRC1 PE=1 SV=1 E(blastx)=2e-23" -- Wnt signaling                                     
# [2] "Peroxidasin homolog OS=Homo sapiens OX=9606 GN=PXDN PE=1 SV=2 E(blastx)=1e-56" --- interesting, may have to do with oxidative stress / immune response                                                                       
# [3] "Short-chain collagen C4 (Fragment) OS=Ephydatia muelleri OX=6052 PE=2 SV=1 E(blastx)=3e-05"                                                          
# [4] "Sushi, von Willebrand factor type A, EGF and pentraxin domain-containing protein 1 OS=Rattus norvegicus OX=10116 GN=Svep1 PE=1 SV=1 E(blastx)=4e-124" -- cell adhesion
# [5] "Low-density lipoprotein receptor-related protein 6 OS=Homo sapiens OX=9606 GN=LRP6 PE=1 SV=2 E(blastx)=5e-47"  --- Wnt signaling

# GN=CTHRC1 and GN=LRP6 have to do with Wnt/beta-catenin signaling (plays a role in bone formation) so maybe forming CaCO3 skeleton?
# recruits could be upreg these genes in control bc growing faster but unable to sustain during heat stress whereas
# adult metabolism and skeleton formation is increasing under heat stress? -- need to look at these pathways in coral lit
# Wnt pathway also involved in development but interesting that it would be upred in adults under heat stress

# backloading in recruits
candidates=list('adult heat response down' = AgenesOnlyDown,'recruit vs adult control down' = stage_control_down)
venn(candidates) # 29 genes backloaded by recruits, but responding in adults
recruit_back=iso2gene[iso2gene$V1 %in% AgenesOnlyDown[AgenesOnlyDown %in% stage_control_down],]
# [1] "Polycystic kidney disease and receptor for egg jelly-related protein OS=Mus musculus OX=10090 GN=Pkdrej PE=2 SV=1 E(blastx)=8e-79"
# [2] "Anoctamin-4 OS=Bos taurus OX=9913 GN=ANO4 PE=2 SV=1 E(blastx)=6e-171"                                                             
# [3] "Basic phospholipase A2 S6-45 OS=Austrelaps superbus OX=29156 PE=2 SV=1 E(blastx)=2e-12"                                           
# [4] "Cryptochrome-1 OS=Gallus gallus OX=9031 GN=CRY1 PE=2 SV=1 E(blastx)=2e-19"                                                        
# [5] "Pancreatic lipase-related protein 2 OS=Cavia porcellus OX=10141 GN=PNLIPRP2 PE=1 SV=1 E(blastx)=1e-58"                

####################################
# Now for adults
up_control=subset(res_stage_control,padj < 0.1 & log2FoldChange > 0)
down_control=subset(res_stage_control,padj < 0.1 & log2FoldChange < 0)

stage_control_up=rownames(up_control)
stage_control_down=rownames(down_control)

# frontloading in adults
candidates=list('recruit heat response up' = RgenesOnlyUP,'adult vs recruit control up' = stage_control_up)
venn(candidates) # 44 genes frontloaded by adults, but responding in recruits
adult_front=iso2gene[iso2gene$V1 %in% RgenesOnlyUP[RgenesOnlyUP %in% stage_control_up],] # 17 annotated
# oooh interesting, TNF receptor, and universal stress protein

# backloading in adults
candidates=list('recruit heat response down' = RgenesOnlyDown,'adult vs recruit control down' = stage_control_down)
venn(candidates) # 54 genes frontloaded by adults, but responding in recruits
adult_back=iso2gene[iso2gene$V1 %in% RgenesOnlyDown[RgenesOnlyDown %in% stage_control_down],] # 16 annotated
# maybe ubiquitin-protein ligase interesting?

#export df with all and LFC
all=rbind(recruit_front,recruit_back,adult_front,adult_back)
all$protein=as.character(lapply(strsplit(all$V2, split='OS='),'[',1))
all$GN=as.character(lapply(strsplit(all$V2, split='GN='),'[',2))
all$GN=as.character(lapply(strsplit(all$GN, split='PE='),'[',1))
all$V2=NULL
all$stage=c(rep('recruit',times=10),rep('adult',times=33))
all$loading=c(rep('front',times=5),rep('back',times=5),rep('front',times=17),rep('back',times=16))

control_LFC=subset(res_stage_control,rownames(res_stage_control) %in% all$V1,select=log2FoldChange)
all_order=all[match(rownames(control_LFC),all$V1),]
merge_df=as.data.frame(cbind(all_order,control_LFC))
colnames(merge_df)[6]='Adult/recruit LFC control'
merge_df$'Adult LFC Heat'=res_Adult$log2FoldChange[rownames(res_Adult) %in% all$V1]
merge_df$'Recruit LFC Heat'=res_Recruit$log2FoldChange[rownames(res_Recruit) %in% all$V1]
colnames(merge_df)[1]='isogroup'

write.csv(x=merge_df,file='interesting_genes/GeneNames_frontbackloading_by_stage.csv')


#####################################
# correlation plot
RgenesOnly=Recruit_trmtGenes[! Recruit_trmtGenes %in% Adult_trmtGenes] #267

res_A_Rgenes=subset(res_Adult,rownames(res_Adult) %in% RgenesOnly)
res_R_sub=res_Recruit[rownames(res_Recruit) %in% RgenesOnly,]
ggplot()+geom_point(aes(x=res_R_sub$log2FoldChange,y=res_A_Rgenes$log2FoldChange))+
  geom_abline(intercept = 0,slope=1,linetype='dashed',color='blue')+ylab('Adult response log2(fold change)')+
  xlab('Recruit response log2(fold change)')+ggtitle('Recruit DEGs to treatment')+xlim(-6,6)+ylim(-6,6)

# same thing but for adult degs to trmt
AgenesOnly=Adult_trmtGenes[! Adult_trmtGenes %in% Recruit_trmtGenes] #150

res_R_Agenes=subset(res_Recruit,rownames(res_Recruit) %in% AgenesOnly)
res_A_sub=res_Adult[rownames(res_Adult) %in% AgenesOnly,]
ggplot()+geom_point(aes(x=res_A_sub$log2FoldChange,y=res_R_Agenes$log2FoldChange))+
  geom_abline(intercept = 0,slope=1,linetype='dashed',color='blue')+ylab('Recruit response log2(fold change)')+
  xlab('Adult response log2(fold change)')+ggtitle('Adult DEGs to treatment')+xlim(-5,7)+ylim(-5,7)

# check for frontloading of recruit DEGs in Adults
res_stage_control_rDEGs=subset(res_stage_control,rownames(res_stage_control) %in% RgenesOnlyUP)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_A_Rgenes=subset(res_Adult,rownames(res_Adult) %in% RgenesOnlyUP)
res_R_Rgenes=subset(res_Recruit,rownames(res_Recruit) %in% RgenesOnlyUP)
rownames(res_A_Rgenes) == rownames(res_R_Rgenes)
RvA_response = res_A_Rgenes$log2FoldChange / res_R_Rgenes$log2FoldChange

master_df=as.data.frame(cbind(RvA_response,res_stage_control_rDEGs$log2FoldChange),row.names = rownames(res_stage_control_rDEGs))
colnames(master_df)=c('RvA_response','res_stage_control_rDEGs')

# adult frontloading
ggplot(master_df,aes(x=RvA_response,y=res_stage_control_rDEGs))+
  geom_point()+
  #geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=0,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('Adult vs Recruit control log2(fold change)')+
  xlab('Adult LFC vs Recruit LFC')+ggtitle('Recruit DEGs to treatment ONLY UP')

ggplot(master_df,aes(x=RvA_response,y=2^(res_stage_control_rDEGs)))+
  geom_point()+
  #geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=0,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('mean Adult / Recruit control')+
  xlab('Adult LFC vs Recruit LFC')+ggtitle('Recruit DEGs to treatment ONLY UP')

############## adult backloading
res_stage_control_rDEGs=subset(res_stage_control,rownames(res_stage_control) %in% RgenesOnlyDown)

res_A_Rgenes=subset(res_Adult,rownames(res_Adult) %in% RgenesOnlyDown)
res_R_Rgenes=subset(res_Recruit,rownames(res_Recruit) %in% RgenesOnlyDown)
rownames(res_A_Rgenes) == rownames(res_R_Rgenes)
RvA_response = res_A_Rgenes$log2FoldChange / res_R_Rgenes$log2FoldChange

master_df=as.data.frame(cbind(RvA_response,res_stage_control_rDEGs$log2FoldChange),row.names = rownames(res_stage_control_rDEGs))
colnames(master_df)=c('RvA_response','res_stage_control_rDEGs')

ggplot(master_df,aes(x=RvA_response,y=res_stage_control_rDEGs))+
  geom_point()+
  #geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=0,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('Adult vs Recruit control log2(fold change)')+
  xlab('Adult LFC vs Recruit LFC')+ggtitle('Recruit DEGs to treatment ONLY Down')



###### check for frontloading of Adult DEGs in Recruits
res_stage_control$LFC_RA=res_stage_control$log2FoldChange*-1 # multiply by -1 to make LFC Recruit / Adult

res_stage_control_ADEGs=subset(res_stage_control,rownames(res_stage_control) %in% AgenesOnlyUP)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_R_Agenes=subset(res_Recruit,rownames(res_Recruit) %in% AgenesOnlyUP)
res_A_Agenes=subset(res_Adult,rownames(res_Adult) %in% AgenesOnlyUP)
rownames(res_R_Agenes) == rownames(res_A_Agenes)
AvR_response = res_R_Agenes$log2FoldChange / res_A_Agenes$log2FoldChange

master_df=as.data.frame(cbind(AvR_response,res_stage_control_ADEGs$LFC_RA),row.names = rownames(res_stage_control_ADEGs))
colnames(master_df) = c('AvR_response','res_stage_control_ADEGs')

ggplot(master_df,aes(x=AvR_response,y=res_stage_control_ADEGs))+geom_point()+
  #geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=0,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('Recruit vs Adult control log2(fold change)')+
  xlab('Recruit FC vs Adult FC')+ggtitle('Adult DEGs to treatment ONLY UP')

# recruit backloading
res_stage_control_ADEGs=subset(res_stage_control,rownames(res_stage_control) %in% AgenesOnlyDown)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_R_Agenes=subset(res_Recruit,rownames(res_Recruit) %in% AgenesOnlyDown)
res_A_Agenes=subset(res_Adult,rownames(res_Adult) %in% AgenesOnlyDown)
rownames(res_R_Agenes) == rownames(res_A_Agenes)
AvR_response = res_R_Agenes$log2FoldChange / res_A_Agenes$log2FoldChange

master_df=as.data.frame(cbind(AvR_response,res_stage_control_ADEGs$LFC_RA),row.names = rownames(res_stage_control_ADEGs))
colnames(master_df) = c('AvR_response','res_stage_control_ADEGs')

ggplot(master_df,aes(x=AvR_response,y=res_stage_control_ADEGs))+geom_point()+
  #geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=0,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('Recruit vs Adult control log2(fold change)')+
  xlab('Recruit FC vs Adult FC')+ggtitle('Adult DEGs to treatment ONLY Down')






######### BY ORIGIN #########

# read in data
res=read.csv('DESeq2HostAR_origin_response_to_treat_VsdPvalsLFC_Oct2021.csv',row.names=1)
iso2gene=read.table('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/past_blast_annot_July2021/past_iso2gene.tab',sep = '\t', colClasses = "character")

res_in=subset(res,select=c('in_padj','in_LFC'))
colnames(res_in)=c('padj','log2FoldChange')
res_off=subset(res,select=c('off_padj','off_LFC'))
colnames(res_off)=c('padj','log2FoldChange')
res_origin_control=subset(res,select=c('origin_control_padj','origin_control_LFC')) # LFC offshore / inshore
colnames(res_origin_control)=c('padj','log2FoldChange')

# subset sig treatment responsive genes
in_trmtGenes=rownames(subset(res_in,padj < 0.1)) # for inshore
off_trmtGenes=rownames(subset(res_off,padj < 0.1)) # for offshore

# let's check if there are any inshore DEGs being frontloaded / backloaded by offshore
# aka look at if any of the treatment responsive genes not responding in offshore are upregulated relative to inshore samples in control conditions
IgenesOnlyUP=rownames(subset(res_in,padj < 0.1 & log2FoldChange > 0  & !rownames(res_in) %in% off_trmtGenes)) # sig upreg in inshore
OgenesOnlyUP=rownames(subset(res_off,padj < 0.1 & log2FoldChange > 0  & !rownames(res_off) %in% in_trmtGenes)) # sig upreg in offshore

IgenesOnlyDown=rownames(subset(res_in,padj < 0.1 & log2FoldChange < 0  & !rownames(res_in) %in% off_trmtGenes)) # sig downreg in inshore
OgenesOnlyDown=rownames(subset(res_off,padj < 0.1 & log2FoldChange < 0  & !rownames(res_off) %in% in_trmtGenes)) # sig upreg in offshore

# get genes upreg or downreg in inshore relative to offshore
up_control=subset(res_origin_control,padj < 0.1 & log2FoldChange > 0) # upreg in offshore relative to inshore
down_control=subset(res_origin_control,padj < 0.1 & log2FoldChange < 0) # downreg in offshore relative to inshore

origin_control_up=rownames(up_control) # get gene names
origin_control_down=rownames(down_control)

# frontloading
candidates=list('inshore heat response up' = IgenesOnlyUP,'offshore vs inshore control up' = origin_control_up)
venn(candidates) # 1 genes frontloaded by offshore, but responding in inshore
offshore_front=IgenesOnlyUP[IgenesOnlyUP %in% origin_control_up]
offshore_front_annot=iso2gene[iso2gene$V1 %in% IgenesOnlyUP[IgenesOnlyUP %in% origin_control_up],]
# not annotated

# backloading
candidates=list('inshore heat response down' = IgenesOnlyDown,'offshore vs inshore control down' = origin_control_down)
venn(candidates) # 14 genes backloaded
offshore_back=IgenesOnlyDown[IgenesOnlyDown %in% origin_control_down]
offshore_back_annot=iso2gene[iso2gene$V1 %in% IgenesOnlyDown[IgenesOnlyDown %in% origin_control_down],]
# [1] "Hydrogenase maturation factor HoxX OS=Bradyrhizobium diazoefficiens (strain JCM 10833 / BCRC 13528 / IAM 13628 / NBRC 14792 / USDA 110) OX=224911 GN=hoxX PE=4 SV=2 E(blastx)=2e-96"
# [2] "Phenolphthiocerol/phthiocerol polyketide synthase subunit C OS=Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv) OX=83332 GN=ppsC PE=1 SV=2 E(blastx)=2e-32" -- lipid metabolism                   
# [3] "Collagen triple helix repeat-containing protein 1 OS=Homo sapiens OX=9606 GN=CTHRC1 PE=1 SV=1 E(blastx)=5e-42" -- negative regulator of collagen matrix deposition                                                                      
# [4] "Galanin receptor 2b OS=Danio rerio OX=7955 GN=galr2b PE=2 SV=1 E(blastx)=3e-22" -- Receptor for the hormone galanin                                                                                                     
# [5] "Lactadherin OS=Homo sapiens OX=9606 GN=MFGE8 PE=1 SV=3 E(blastx)=5e-20" -- Angiogenesis, Cell adhesion, Fertilization -- Contributes to phagocytic removal of apoptotic cells in many tissues
off_back_allisos=origin_control_down[origin_control_down %in% IgenesOnlyDown]
off_back_sub=subset(res_origin_control,rownames(res_origin_control) %in% off_back_allisos)
#                                 padj log2FoldChange
# isogroupTR111168_c0_g1 9.094872e-02     -3.8245473
# isogroupTR112556_c1_g2 7.887978e-02     -3.2465272  #HoxX
# isogroupTR121003_c2_g1 3.914427e-02     -1.6661571
# isogroupTR122006_c6_g5 7.555047e-02     -2.6595441
# isogroupTR134046_c3_g6 1.229414e-05     -3.1327488
# isogroupTR139817_c5_g1 5.471869e-02     -0.9021508  #MFGE8
# isogroupTR140353_c0_g1 2.738216e-02     -2.0924526  #galr2b
# isogroupTR152550_c0_g2 3.252937e-02     -2.4633740
# isogroupTR152841_c0_g1 2.047565e-02     -3.3447243  #ppsC
# isogroupTR163925_c1_g1 1.105902e-02     -1.4748941  #CTHRC1
# isogroupTR166699_c0_g1 8.151004e-02     -4.6567614
# isogroupTR172581_c0_g1 7.050413e-02     -0.6201164
# isogroupTR85815_c0_g1  2.595896e-02     -0.3341365
# isogroupTR98695_c7_g1  3.005534e-02     -1.5211170
mean(off_back_sub$log2FoldChange) #-2.28
median(off_back_sub$log2FoldChange) #-2.28

###### plot it out
res_stage_control_inDEGsOnlyUP=subset(res_origin_control,rownames(res_origin_control) %in% IgenesOnlyUP)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_in_ingenesOnly=subset(res_in,rownames(res_in) %in% IgenesOnlyUP)
res_off_ingenesOnly=subset(res_off,rownames(res_off) %in% IgenesOnlyUP)
rownames(res_in_ingenesOnly) == rownames(res_off_ingenesOnly)
OvI_response = res_off_ingenesOnly$log2FoldChange / res_in_ingenesOnly$log2FoldChange

master_df=as.data.frame(cbind(OvI_response,res_stage_control_inDEGsOnlyUP$log2FoldChange),row.names = rownames(res_stage_control_inDEGsOnlyUP))
colnames(master_df) = c('OvI_response','res_stage_control_inDEGsOnlyUP')

ggplot(master_df,aes(x=OvI_response,y=res_stage_control_inDEGsOnlyUP))+geom_point()+
  geom_hline(yintercept=1,linetype='dashed',color='orange',size=1.5)+
  geom_vline(xintercept=1,linetype='dashed',color='orange',size=1.5)+
  ylab('Offshore vs Inshore control (LFC)')+
  xlab('Offshore LFC / Inshore LFC (heat vs control)')+
  theme_bw(base_size = 18)+ggtitle('inshore DEGs up')

# now backloading
###### plot it out
res_stage_control_inDEGsOnlyDown=subset(res_origin_control,rownames(res_origin_control) %in% IgenesOnlyDown)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_in_ingenesOnly=subset(res_in,rownames(res_in) %in% IgenesOnlyDown)
res_off_ingenesOnly=subset(res_off,rownames(res_off) %in% IgenesOnlyDown)
rownames(res_in_ingenesOnly) == rownames(res_off_ingenesOnly)
OvI_response = res_off_ingenesOnly$log2FoldChange / res_in_ingenesOnly$log2FoldChange

master_df=as.data.frame(cbind(OvI_response,res_stage_control_inDEGsOnlyDown$log2FoldChange),row.names = rownames(res_stage_control_inDEGsOnlyDown))
colnames(master_df) = c('OvI_response','res_stage_control_inDEGsOnlyDown')

ggplot(master_df,aes(x=OvI_response,y=res_stage_control_inDEGsOnlyDown))+geom_point()+
  geom_hline(yintercept=-1,linetype='dashed',color='orange',size=1.5)+
  geom_vline(xintercept=1,linetype='dashed',color='orange',size=1.5)+
  ylab('Offshore vs Inshore control (LFC)')+
  xlab('Offshore LFC / Inshore LFC (heat vs control)')+
  theme_bw(base_size = 18)+ggtitle('inshore DEGs down')


#### now frontloading in inshore
# get genes upreg or downreg in inshore relative to offshore
up_control=subset(res_origin_control,padj < 0.1 & log2FoldChange < 0) # upreg in inshore relative to offshore
down_control=subset(res_origin_control,padj < 0.1 & log2FoldChange > 0) # downreg in inshore relative to offshore

origin_control_up=rownames(up_control) # get gene names
origin_control_down=rownames(down_control)

# frontloading
candidates=list('offshore heat response up' = OgenesOnlyUP,'offshore vs inshore control up' = origin_control_up)
venn(candidates) # 0 genes

# backloading
candidates=list('offshore heat response down' = OgenesOnlyDown,'offshore vs inshore control down' = origin_control_down)
venn(candidates) # 1 genes
in_back=OgenesOnlyDown[OgenesOnlyDown %in% origin_control_down]
in_back_annot=iso2gene[iso2gene$V1 %in% OgenesOnlyDown[OgenesOnlyDown %in% origin_control_down],] #not annotated

####### export df with all and LFC
#in_front=as.data.frame(cbind('isogroup'=in_front,'origin'=rep('in',times=length(in_front)),'loading'=rep('front',times=length(in_front))))
in_back=as.data.frame(cbind('isogroup'=in_back,'origin'=rep('in',times=length(in_back)),'loading'=rep('back',times=length(in_back))))
offshore_front=as.data.frame(cbind('isogroup'=offshore_front,'origin'=rep('off',times=length(offshore_front)),'loading'=rep('front',times=length(offshore_front))))
offshore_back=as.data.frame(cbind('isogroup'=offshore_back,'origin'=rep('off',times=length(offshore_back)),'loading'=rep('back',times=length(offshore_back))))

all=as.data.frame(rbind(in_back,offshore_front,offshore_back))
colnames(iso2gene)=c('isogroup','annot')
all=merge(all,iso2gene,by='isogroup',all.x=T)
all$protein=as.character(lapply(strsplit(all$annot, split='OS='),'[',1))
all$GN=as.character(lapply(strsplit(all$annot, split='GN='),'[',2))
all$GN=as.character(lapply(strsplit(all$GN, split='PE='),'[',1))
all$annot=NULL

control_LFC=subset(res_origin_control,rownames(res_origin_control) %in% all$isogroup,select=log2FoldChange)
all_order=all[match(rownames(control_LFC),all$isogroup),]
merge_df=as.data.frame(cbind(all_order,control_LFC))
colnames(merge_df)[6]='Off/In LFC control'
merge_df$'In LFC Heat'=res_in$log2FoldChange[rownames(res_in) %in% all$isogroup]
merge_df$'Off LFC Heat'=res_off$log2FoldChange[rownames(res_off) %in% all$isogroup]

write.csv(x=merge_df,file='interesting_genes/HostGeneNames_frontbackloading_by_origin.csv')





#
#
#
#
#
#
#
##### SYMBIONTS STAGE ######
setwd('/Users/Maria/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/downstream_analyses/Mar2020/DESeq/DESeq2_June2021/')

res_Adult=read.csv("DESeq2_symAdult_trmt_response.csv",row.names = 1)
res_Recruit=read.csv("DESeq2_symRecruit_trmt_response.csv", row.names = 1)
res_stage_control=read.csv('DESeq2_sym_stage_control_AR.csv',row.names = 1) # note this is adult / recruit control LFC

Adult_trmtGenes=rownames(subset(res_Adult,padj < 0.1))
Recruit_trmtGenes=rownames(subset(res_Recruit,padj < 0.1))

# let's check if there are any recruit DEGs being frontloaded / backloaded
# aka look at if any of the treatment responsive genes not responding in recruits are upregulated relative to adult samples in control conditions
AgenesOnlyUP=rownames(subset(res_Adult,padj < 0.1 & log2FoldChange > 0  & !rownames(res_Adult) %in% Recruit_trmtGenes))
RgenesOnlyUP=rownames(subset(res_Recruit,padj < 0.1 & log2FoldChange > 0  & !rownames(res_Recruit) %in% Adult_trmtGenes))

AgenesOnlyDown=rownames(subset(res_Adult,padj < 0.1 & log2FoldChange < 0  & !rownames(res_Adult) %in% Recruit_trmtGenes))
RgenesOnlyDown=rownames(subset(res_Recruit,padj < 0.1 & log2FoldChange < 0  & !rownames(res_Recruit) %in% Adult_trmtGenes))

# get genes upreg or downreg in recruits relative to adults
up_control=subset(res_stage_control,padj < 0.1 & log2FoldChange < 0) # negative = upreg in recruits
down_control=subset(res_stage_control,padj < 0.1 & log2FoldChange > 0)

stage_control_up=rownames(up_control)
stage_control_down=rownames(down_control)

# frontloading
candidates=list('adult heat response up' = AgenesOnlyUP,'recruit vs adult control up' = stage_control_up)
venn(candidates) # 22 genes frontloaded by recruits, but responding in adults
recruit_front=AgenesOnlyUP[AgenesOnlyUP %in% stage_control_up]

# backloading
candidates=list('adult heat response down' = AgenesOnlyDown,'recruit vs adult control down' = stage_control_down)
venn(candidates) # 13 genes backloaded by recruits, but responding in adults
recruit_back=AgenesOnlyDown[AgenesOnlyDown %in% stage_control_down]

iso2gene=read.table('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/kb8_blast_annot_July2021/kb8_iso2gene.tab',sep = '\t', colClasses = "character")
recruit_frontannot=iso2gene[iso2gene$V1 %in% recruit_front,] 
# [1] "Probable steroid-binding protein 3 OS=Arabidopsis thaliana OX=3702 GN=MP3 PE=1 SV=1 E(blastx)=3e-09"                       
# [2] "Regulatory protein FlaEY OS=Caulobacter vibrioides (strain ATCC 19089 / CB15) OX=190650 GN=flaEY PE=4 SV=2 E(blastx)=4e-07"
# [3] "Dextranase OS=Arthrobacter globiformis OX=1665 PE=3 SV=1 E(blastx)=2e-11"
recruit_backannot=iso2gene[iso2gene$V1 %in% recruit_back,]
# [1] "Cytochrome c-550 OS=Phaeodactylum tricornutum (strain CCAP 1055/1) OX=556484 GN=psbV PE=3 SV=1 E(blastx)=2e-52"
# [2] "Carnosine synthase 1 OS=Gallus gallus OX=9031 GN=CARNS1 PE=1 SV=1 E(blastx)=5e-16"                             
# [3] "Major basic nuclear protein 2 OS=Crypthecodinium cohnii OX=2866 GN=HCc2 PE=4 SV=1 E(blastx)=2e-08"             
# [4] "Photosystem I reaction center subunit II OS=Trieres chinensis OX=1514140 GN=psaD PE=3 SV=1 E(blastx)=2e-44"
# recruits backloading photosynthesis gene, but is responsive in adults

####################################
# Now for adults
# up in adults, down in recruits
up_control=subset(res_stage_control,padj < 0.1 & log2FoldChange > 0)
down_control=subset(res_stage_control,padj < 0.1 & log2FoldChange < 0)

#up_control=subset(res_stage_control,padj < 0.1 & log2FoldChange > 2)
#down_control=subset(res_stage_control,padj < 0.1 & log2FoldChange < -2)

stage_control_up=rownames(up_control)
stage_control_down=rownames(down_control)

# frontloading
candidates=list('recruit heat response up' = RgenesOnlyUP,'adult vs recruit control up' = stage_control_up)
venn(candidates) # 14 genes frontloaded by adults, but responding in recruits
# 0 with LFC > 2 in control threshold
adult_front=RgenesOnlyUP[RgenesOnlyUP %in% stage_control_up]
adult_frontannot=iso2gene[iso2gene$V1 %in% adult_front,]
# [1] "Regulator of rDNA transcription protein 15 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=RRT15 PE=3 SV=1 E(blastx)=4e-13" -- transcriptional regulation
# [2] "CDGSH iron-sulfur domain-containing protein 1 OS=Bos taurus OX=9913 GN=CISD1 PE=1 SV=1 E(blastx)=6e-18" -- regulates cellulcar respiration                                         
# [3] "Cytochrome c oxidase subunit 1 OS=Plasmodium falciparum OX=5833 GN=MT-CO1 PE=3 SV=1 E(blastx)=3e-75" -- also involved in cellular respiration, electron transport chain                                            
# [4] "Glutathione S-transferase A2 OS=Mus musculus OX=10090 GN=Gsta2 PE=1 SV=3 E(blastx)=5e-08"                                                        
# [5] "Thiol protease SEN102 OS=Hemerocallis sp. OX=29711 GN=SEN102 PE=2 SV=1 E(blastx)=9e-49"

# backloading
candidates=list('recruit heat response down' = RgenesOnlyDown,'adult vs recruit control down' = stage_control_down)
venn(candidates) # 27 genes frontloaded by adults, but responding in recruits
# 0 with LFC < -2 in control
adult_back=RgenesOnlyDown[RgenesOnlyDown %in% stage_control_down]
adult_backannot=iso2gene[iso2gene$V1 %in% adult_back,] # 4 annotated
# [1] "Cytochrome P450 704C1 OS=Pinus taeda OX=3352 GN=CYP704C1 PE=2 SV=1 E(blastx)=2e-57"                                                               
# [2] "3-isopropylmalate dehydrogenase OS=Haemophilus influenzae (strain ATCC 51907 / DSM 11121 / KW20 / Rd) OX=71421 GN=leuB PE=1 SV=1 E(blastx)=3e-113" -- leucine biosynthesis
# [3] "ATP synthase subunit b OS=Nostoc punctiforme (strain ATCC 29133 / PCC 73102) OX=63737 GN=atpF PE=3 SV=1 E(blastx)=8e-09"                          
# [4] "NADP-dependent glyceraldehyde-3-phosphate dehydrogenase OS=Streptococcus equinus OX=1335 GN=gapN PE=1 SV=1 E(blastx)=4e-35"

#export df with all and LFC
recruit_front=as.data.frame(cbind('isogroup'=recruit_front,'stage'=rep('recruit',times=length(recruit_front)),'loading'=rep('front',times=length(recruit_front))))
recruit_back=as.data.frame(cbind('isogroup'=recruit_back,'stage'=rep('recruit',times=length(recruit_back)),'loading'=rep('back',times=length(recruit_back))))
adult_front=as.data.frame(cbind('isogroup'=adult_front,'stage'=rep('adult',times=length(adult_front)),'loading'=rep('front',times=length(adult_front))))
adult_back=as.data.frame(cbind('isogroup'=adult_back,'stage'=rep('adult',times=length(adult_back)),'loading'=rep('back',times=length(adult_back))))

all=as.data.frame(rbind(recruit_front,recruit_back,adult_front,adult_back))
colnames(iso2gene)=c('isogroup','annot')
all=merge(all,iso2gene,by='isogroup',all.x=T)
all$protein=as.character(lapply(strsplit(all$annot, split='OS='),'[',1))
all$GN=as.character(lapply(strsplit(all$annot, split='GN='),'[',2))
all$GN=as.character(lapply(strsplit(all$GN, split='PE='),'[',1))
all$annot=NULL

control_LFC=subset(res_stage_control,rownames(res_stage_control) %in% all$isogroup,select=log2FoldChange)
all_order=all[match(rownames(control_LFC),all$isogroup),]
merge_df=as.data.frame(cbind(all_order,control_LFC))
colnames(merge_df)[6]='symAdult/symrecruit LFC control'
merge_df$'symAdult LFC Heat'=res_Adult$log2FoldChange[rownames(res_Adult) %in% all$isogroup]
merge_df$'symRecruit LFC Heat'=res_Recruit$log2FoldChange[rownames(res_Recruit) %in% all$isogroup]

write.csv(x=merge_df,file='interesting_genes/SymGeneNames_frontbackloading_by_stage.csv')

#####################################
# correlation plot
RgenesOnly=Recruit_trmtGenes[! Recruit_trmtGenes %in% Adult_trmtGenes] #139

res_A_Rgenes=subset(res_Adult,rownames(res_Adult) %in% RgenesOnly)
res_R_sub=res_Recruit[rownames(res_Recruit) %in% RgenesOnly,]
ggplot()+geom_point(aes(x=res_R_sub$log2FoldChange,y=res_A_Rgenes$log2FoldChange))+
  geom_abline(intercept = 0,slope=1,linetype='dashed',color='blue')+ylab('Adult response log2(fold change)')+
  xlab('Recruit response log2(fold change)')+ggtitle('Recruit DEGs to treatment')#+xlim(-6,6)+ylim(-6,6)

# same thing but for adult degs to trmt
AgenesOnly=Adult_trmtGenes[! Adult_trmtGenes %in% Recruit_trmtGenes] #76

res_R_Agenes=subset(res_Recruit,rownames(res_Recruit) %in% AgenesOnly)
res_A_sub=res_Adult[rownames(res_Adult) %in% AgenesOnly,]
ggplot()+geom_point(aes(x=res_A_sub$log2FoldChange,y=res_R_Agenes$log2FoldChange))+
  geom_abline(intercept = 0,slope=1,linetype='dashed',color='blue')+ylab('Recruit response log2(fold change)')+
  xlab('Adult response log2(fold change)')+ggtitle('Adult DEGs to treatment')+xlim(-3,3)+ylim(-2,2)

# check for frontloading of recruit DEGs in Adults
res_stage_control_rDEGs=subset(res_stage_control,rownames(res_stage_control) %in% RgenesOnlyUP)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_A_Rgenes=subset(res_Adult,rownames(res_Adult) %in% RgenesOnlyUP)
res_R_Rgenes=subset(res_Recruit,rownames(res_Recruit) %in% RgenesOnlyUP)
rownames(res_A_Rgenes) == rownames(res_R_Rgenes)
RvA_response = res_A_Rgenes$log2FoldChange / res_R_Rgenes$log2FoldChange

master_df=as.data.frame(cbind(RvA_response,res_stage_control_rDEGs$log2FoldChange),row.names = rownames(res_stage_control_rDEGs))
colnames(master_df)=c('RvA_response','res_stage_control_rDEGs')

ggplot(master_df,aes(x=RvA_response,y=res_stage_control_rDEGs))+
  geom_point()+
  geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=1,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('Adult vs Recruit control log2(fold change)')+
  xlab('Adult FC vs Recruit FC')+ggtitle('Recruit DEGs to treatment ONLY UP')


###### check for frontloading of Adult DEGs in Recruits
res_stage_control$LFC_RA=res_stage_control$log2FoldChange*-1 # multiply by -1 to make LFC Recruit / Adult

res_stage_control_ADEGs=subset(res_stage_control,rownames(res_stage_control) %in% AgenesOnlyUP)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_R_Agenes=subset(res_Recruit,rownames(res_Recruit) %in% AgenesOnlyUP)
res_A_Agenes=subset(res_Adult,rownames(res_Adult) %in% AgenesOnlyUP)
rownames(res_R_Agenes) == rownames(res_A_Agenes)
AvR_response = res_R_Agenes$log2FoldChange / res_A_Agenes$log2FoldChange

master_df=as.data.frame(cbind(AvR_response,res_stage_control_ADEGs$LFC_RA),row.names = rownames(res_stage_control_ADEGs))
colnames(master_df) = c('AvR_response','res_stage_control_ADEGs')

ggplot(master_df,aes(x=AvR_response,y=res_stage_control_ADEGs))+geom_point()+
  geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=1,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('Recruit vs Adult control log2(fold change)')+
  xlab('Recruit FC vs Adult FC')+ggtitle('Adult DEGs to treatment ONLY')


######### SYMBIONTS BY ORIGIN #########

res_in=read.csv("DESeq2_symInshore_trmt_response.csv",row.names = 1)
res_off=read.csv("DESeq2_symOffshore_trmt_response.csv", row.names = 1)
res_origin_control=read.csv('DESeq2_sym_origin_control_IO.csv',row.names = 1) # note this is inshore / offshore control LFC

in_trmtGenes=rownames(subset(res_in,padj < 0.1))
off_trmtGenes=rownames(subset(res_off,padj < 0.1))

# let's check if there are any off DEGs being frontloaded / backloaded
# aka look at if any of the treatment responsive genes not responding in offs are upregulated relative to in samples in control conditions
IgenesOnlyUP=rownames(subset(res_in,padj < 0.1 & log2FoldChange > 0  & !rownames(res_in) %in% off_trmtGenes))
OgenesOnlyUP=rownames(subset(res_off,padj < 0.1 & log2FoldChange > 0  & !rownames(res_off) %in% in_trmtGenes))

IgenesOnlyDown=rownames(subset(res_in,padj < 0.1 & log2FoldChange < 0  & !rownames(res_in) %in% off_trmtGenes))
OgenesOnlyDown=rownames(subset(res_off,padj < 0.1 & log2FoldChange < 0  & !rownames(res_off) %in% in_trmtGenes))

# get genes upreg or downreg in offs relative to ins
up_control=subset(res_origin_control,padj < 0.1 & log2FoldChange < 0)
down_control=subset(res_origin_control,padj < 0.1 & log2FoldChange > 0)

origin_control_up=rownames(up_control)
origin_control_down=rownames(down_control)

# frontloading
candidates=list('in heat response up' = IgenesOnlyUP,'off vs in control up' = origin_control_up)
venn(candidates) # 11 genes frontloaded by offs, but responding in ins
off_front=IgenesOnlyUP[IgenesOnlyUP %in% origin_control_up]
off_annots=iso2gene[iso2gene$V1 %in% IgenesOnlyUP[IgenesOnlyUP %in% origin_control_up],]
#[1] "Dapdiamide synthesis protein DdaC OS=Enterobacter agglomerans OX=549 GN=ddaC PE=4 SV=1 E(blastx)=2e-07"

# backloading
candidates=list('in heat response down' = IgenesOnlyDown,'off vs in control down' = origin_control_down)
venn(candidates) # 1 genes backloaded by offs, but responding in ins
off_back=IgenesOnlyDown[IgenesOnlyDown %in% origin_control_down]
off_annot=iso2gene[iso2gene$V1 %in% IgenesOnlyDown[IgenesOnlyDown %in% origin_control_down],]
# gene not annotated

####################################
# Now for ins
up_control=subset(res_origin_control,padj < 0.1 & log2FoldChange > 0)
down_control=subset(res_origin_control,padj < 0.1 & log2FoldChange < 0)

origin_control_up=rownames(up_control)
origin_control_down=rownames(down_control)

# frontloading
candidates=list('off heat response up' = OgenesOnlyUP,'in vs off control up' = origin_control_up)
venn(candidates) # 1 genes frontloaded by ins, but responding in offs
in_front=OgenesOnlyUP[OgenesOnlyUP %in% origin_control_up]
in_front_annot=iso2gene[iso2gene$V1 %in% OgenesOnlyUP[OgenesOnlyUP %in% origin_control_up],]
# not annotated

# backloading
candidates=list('off heat response down' = OgenesOnlyDown,'in vs off control down' = origin_control_down)
venn(candidates) # 2 genes frontloaded by ins, but responding in offs
in_back=OgenesOnlyDown[OgenesOnlyDown %in% origin_control_down]
in_back_annot=iso2gene[iso2gene$V1 %in% OgenesOnlyDown[OgenesOnlyDown %in% origin_control_down],]
# not annotated

#export df with all and LFC
in_front=as.data.frame(cbind('isogroup'=in_front,'origin'=rep('in',times=length(in_front)),'loading'=rep('front',times=length(in_front))))
in_back=as.data.frame(cbind('isogroup'=in_back,'origin'=rep('in',times=length(in_back)),'loading'=rep('back',times=length(in_back))))
off_front=as.data.frame(cbind('isogroup'=off_front,'origin'=rep('off',times=length(off_front)),'loading'=rep('front',times=length(off_front))))
off_back=as.data.frame(cbind('isogroup'=off_back,'origin'=rep('off',times=length(off_back)),'loading'=rep('back',times=length(off_back))))

all=as.data.frame(rbind(in_front,in_back,off_front,off_back))
colnames(iso2gene)=c('isogroup','annot')
all=merge(all,iso2gene,by='isogroup',all.x=T)
all$protein=as.character(lapply(strsplit(all$annot, split='OS='),'[',1))
all$GN=as.character(lapply(strsplit(all$annot, split='GN='),'[',2))
all$GN=as.character(lapply(strsplit(all$GN, split='PE='),'[',1))
all$annot=NULL

control_LFC=subset(res_origin_control,rownames(res_origin_control) %in% all$isogroup,select=log2FoldChange)
all_order=all[match(rownames(control_LFC),all$isogroup),]
merge_df=as.data.frame(cbind(all_order,control_LFC))
colnames(merge_df)[6]='symInshore/symOffshore LFC control'
merge_df$'symInshore LFC Heat'=res_in$log2FoldChange[rownames(res_in) %in% all$isogroup]
merge_df$'symOffshore LFC Heat'=res_off$log2FoldChange[rownames(res_off) %in% all$isogroup]

write.csv(x=merge_df,file='interesting_genes/SymGeneNames_frontbackloading_by_origin.csv')


#####################################
# correlation plot
OgenesOnly=off_trmtGenes[! off_trmtGenes %in% in_trmtGenes] #267

res_I_Ogenes=subset(res_in,rownames(res_in) %in% OgenesOnly)
res_O_sub=res_off[rownames(res_off) %in% OgenesOnly,]
ggplot()+geom_point(aes(x=res_O_sub$log2FoldChange,y=res_I_Ogenes$log2FoldChange))+
  geom_abline(intercept = 0,slope=1,linetype='dashed',color='blue')+ylab('in response log2(fold change)')+
  xlab('off response log2(fold change)')+ggtitle('off DEGs to treatment')+xlim(-2,2)+ylim(-2,2)

# same thing but for in degs to trmt
IgenesOnly=in_trmtGenes[! in_trmtGenes %in% off_trmtGenes] #150

res_O_Igenes=subset(res_off,rownames(res_off) %in% IgenesOnly)
res_I_sub=res_in[rownames(res_in) %in% IgenesOnly,]
ggplot()+geom_point(aes(x=res_I_sub$log2FoldChange,y=res_O_Igenes$log2FoldChange))+
  geom_abline(intercept = 0,slope=1,linetype='dashed',color='blue')+ylab('off response log2(fold change)')+
  xlab('in response log2(fold change)')+ggtitle('in DEGs to treatment')+xlim(-3,3)+ylim(-3,3)

###### check for frontloading of off DEGs in ins
res_origin_control_offDEGsUP=subset(res_origin_control,rownames(res_origin_control) %in% OgenesOnlyUP)

# compare off foldchange to heat / in fold change to heat vs foldchange of off vs in in control
res_I_Ogenes=subset(res_in,rownames(res_in) %in% OgenesOnlyUP)
res_O_Ogenes=subset(res_off,rownames(res_off) %in% OgenesOnlyUP)
rownames(res_I_Ogenes) == rownames(res_O_Ogenes)
OvI_response = res_I_Ogenes$log2FoldChange / res_O_Ogenes$log2FoldChange

master_df=as.data.frame(cbind(OvI_response,res_origin_control_offDEGsUP$log2FoldChange),row.names = rownames(res_origin_control_offDEGsUP))
colnames(master_df)=c('OvI_response','res_origin_control_offDEGs')

ggplot(master_df,aes(x=OvI_response,y=res_origin_control_offDEGs))+
  geom_point()+
  #geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=1,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('in vs off control log2(fold change)')+
  xlab('in FC vs off FC')+ggtitle('off DEGs to treatment ONLY UP')

# now backloading
res_origin_control_offDEGsDown=subset(res_origin_control,rownames(res_origin_control) %in% OgenesOnlyDown)
res_I_Ogenes=subset(res_in,rownames(res_in) %in% OgenesOnlyDown)
res_O_Ogenes=subset(res_off,rownames(res_off) %in% OgenesOnlyDown)
rownames(res_I_Ogenes) == rownames(res_O_Ogenes)
OvI_response = res_I_Ogenes$log2FoldChange / res_O_Ogenes$log2FoldChange

master_df=as.data.frame(cbind(OvI_response,res_origin_control_offDEGsDown$log2FoldChange),row.names = rownames(res_origin_control_offDEGsDown))
colnames(master_df)=c('OvI_response','res_origin_control_offDEGsDown')

ggplot(master_df,aes(x=OvI_response,y=res_origin_control_offDEGsDown))+
  geom_point()+
  #geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=-1,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('in vs off control log2(fold change)')+
  xlab('in FC vs off FC')+ggtitle('off DEGs to treatment ONLY Down')


###### check for frontloading of in DEGs in offs
res_origin_control$LFC_OI=res_origin_control$log2FoldChange*-1 # multiply by -1 to make LFC off / in

res_origin_control_IDEGsUP=subset(res_origin_control,rownames(res_origin_control) %in% IgenesOnlyUP)

# compare off foldchange to heat / in fold change to heat vs foldchange of off vs in in control
res_O_Igenes=subset(res_off,rownames(res_off) %in% IgenesOnlyUP)
res_I_Igenes=subset(res_in,rownames(res_in) %in% IgenesOnlyUP)
rownames(res_O_Igenes) == rownames(res_I_Igenes)
IvO_response = res_O_Igenes$log2FoldChange / res_I_Igenes$log2FoldChange

master_df=as.data.frame(cbind(IvO_response,res_origin_control_IDEGsUP$LFC_OI),row.names = rownames(res_origin_control_IDEGsUP))
colnames(master_df) = c('IvO_response','res_origin_control_IDEGsUP')

ggplot(master_df,aes(x=IvO_response,y=res_origin_control_IDEGsUP))+geom_point()+
  #geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=1,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('off vs in control log2(fold change)')+
  xlab('off FC vs in FC')+ggtitle('in DEGs to treatment ONLY UP')

# now backloading
res_origin_control_IDEGsDown=subset(res_origin_control,rownames(res_origin_control) %in% IgenesOnlyDown)

# compare off foldchange to heat / in fold change to heat vs foldchange of off vs in in control
res_O_Igenes=subset(res_off,rownames(res_off) %in% IgenesOnlyDown)
res_I_Igenes=subset(res_in,rownames(res_in) %in% IgenesOnlyDown)
rownames(res_O_Igenes) == rownames(res_I_Igenes)
IvO_response = res_O_Igenes$log2FoldChange / res_I_Igenes$log2FoldChange

master_df=as.data.frame(cbind(IvO_response,res_origin_control_IDEGsDown$LFC_OI),row.names = rownames(res_origin_control_IDEGsDown))
colnames(master_df) = c('IvO_response','res_origin_control_IDEGsDown')

ggplot(master_df,aes(x=IvO_response,y=res_origin_control_IDEGsDown))+geom_point()+
  #geom_smooth(method='lm',se=F)+
  geom_hline(yintercept=-1,linetype='dashed',color='red')+
  geom_vline(xintercept=1,linetype='dashed',color='red')+
  ylab('off vs in control log2(fold change)')+
  xlab('off FC vs in FC')+ggtitle('in DEGs to treatment ONLY Down')



######## For larvae origin specific responses ######

# read in data
res=read.csv('LarvDESeq2_VsdPvalsLFC_byGroup_19July2021.csv',row.names = 1)
iso2gene=read.table('~/Desktop/Kenkel_lab/bioinformatics/bioinformatics_class_proj/computational/past_blast_annot_July2021/past_iso2gene.tab',sep = '\t', colClasses = "character")

res_in=subset(res,select=c('inTrmt_padj','inTrmt_LFC'))
colnames(res_in)=c('padj','log2FoldChange')
res_off=subset(res,select=c('offTrmt_padj','offTrmt_LFC'))
colnames(res_off)=c('padj','log2FoldChange')
res_origin_control=subset(res,select=c('originControl_padj','originControl_LFC')) # LFC offshore / inshore
colnames(res_origin_control)=c('padj','log2FoldChange')

# subset sig treatment responsive genes
in_trmtGenes=rownames(subset(res_in,padj < 0.1)) # for inshore
off_trmtGenes=rownames(subset(res_off,padj < 0.1)) # for offshore

# let's check if there are any inshore DEGs being frontloaded / backloaded by offshore
# aka look at if any of the treatment responsive genes not responding in offshore are upregulated relative to inshore samples in control conditions
IgenesOnlyUP=rownames(subset(res_in,padj < 0.1 & log2FoldChange > 0  & !rownames(res_in) %in% off_trmtGenes)) # sig upreg in inshore
OgenesOnlyUP=rownames(subset(res_off,padj < 0.1 & log2FoldChange > 0  & !rownames(res_off) %in% in_trmtGenes)) # sig upreg in offshore

IgenesOnlyDown=rownames(subset(res_in,padj < 0.1 & log2FoldChange < 0  & !rownames(res_in) %in% off_trmtGenes)) # sig downreg in inshore
OgenesOnlyDown=rownames(subset(res_off,padj < 0.1 & log2FoldChange < 0  & !rownames(res_off) %in% in_trmtGenes)) # sig upreg in offshore

# get genes upreg or downreg in inshore relative to offshore
up_control=subset(res_origin_control,padj < 0.1 & log2FoldChange > 0) # upreg in offshore relative to inshore
down_control=subset(res_origin_control,padj < 0.1 & log2FoldChange < 0) # downreg in offshore relative to inshore

origin_control_up=rownames(up_control) # get gene names
origin_control_down=rownames(down_control)

# frontloading
candidates=list('inshore heat response up' = IgenesOnlyUP,'offshore vs inshore control up' = origin_control_up)
venn(candidates) # 0 genes frontloaded by offshore, but responding in inshore

# backloading
candidates=list('inshore heat response down' = IgenesOnlyDown,'offshore vs inshore control down' = origin_control_down)
venn(candidates) # 0 genes backloaded

###### plot it out
res_stage_control_inDEGsOnlyUP=subset(res_origin_control,rownames(res_origin_control) %in% IgenesOnlyUP)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_in_ingenesOnly=subset(res_in,rownames(res_in) %in% IgenesOnlyUP)
res_off_ingenesOnly=subset(res_off,rownames(res_off) %in% IgenesOnlyUP)
rownames(res_in_ingenesOnly) == rownames(res_off_ingenesOnly)
OvI_response = res_off_ingenesOnly$log2FoldChange / res_in_ingenesOnly$log2FoldChange

master_df=as.data.frame(cbind(OvI_response,res_stage_control_inDEGsOnlyUP$log2FoldChange),row.names = rownames(res_stage_control_inDEGsOnlyUP))
colnames(master_df) = c('OvI_response','res_stage_control_inDEGsOnlyUP')

ggplot(master_df,aes(x=OvI_response,y=res_stage_control_inDEGsOnlyUP))+geom_point()+
  geom_hline(yintercept=1,linetype='dashed',color='orange',size=1.5)+
  geom_vline(xintercept=1,linetype='dashed',color='orange',size=1.5)+
  ylab('Offshore vs Inshore control (LFC)')+
  xlab('Offshore LFC / Inshore LFC (heat vs control)')+
  theme_bw(base_size = 18)+ggtitle('inshore DEGs up')

# now backloading
###### plot it out
res_stage_control_inDEGsOnlyDown=subset(res_origin_control,rownames(res_origin_control) %in% IgenesOnlyDown)

# compare recruit foldchange to heat / adult fold change to heat vs foldchange of Recruit vs Adult in control
res_in_ingenesOnly=subset(res_in,rownames(res_in) %in% IgenesOnlyDown)
res_off_ingenesOnly=subset(res_off,rownames(res_off) %in% IgenesOnlyDown)
rownames(res_in_ingenesOnly) == rownames(res_off_ingenesOnly)
OvI_response = res_off_ingenesOnly$log2FoldChange / res_in_ingenesOnly$log2FoldChange

master_df=as.data.frame(cbind(OvI_response,res_stage_control_inDEGsOnlyDown$log2FoldChange),row.names = rownames(res_stage_control_inDEGsOnlyDown))
colnames(master_df) = c('OvI_response','res_stage_control_inDEGsOnlyDown')

ggplot(master_df,aes(x=OvI_response,y=res_stage_control_inDEGsOnlyDown))+geom_point()+
  geom_hline(yintercept=-1,linetype='dashed',color='orange',size=1.5)+
  geom_vline(xintercept=1,linetype='dashed',color='orange',size=1.5)+
  ylab('Offshore vs Inshore control (LFC)')+
  xlab('Offshore LFC / Inshore LFC (heat vs control)')+
  theme_bw(base_size = 18)+ggtitle('inshore DEGs down')


#### now frontloading in inshore
# get genes upreg or downreg in inshore relative to offshore
up_control=subset(res_origin_control,padj < 0.1 & log2FoldChange < 0) # upreg in inshore relative to offshore
down_control=subset(res_origin_control,padj < 0.1 & log2FoldChange > 0) # downreg in inshore relative to offshore

origin_control_up=rownames(up_control) # get gene names
origin_control_down=rownames(down_control)

# frontloading
candidates=list('offshore heat response up' = OgenesOnlyUP,'offshore vs inshore control up' = origin_control_up)
venn(candidates) # 0 genes

# backloading
candidates=list('offshore heat response down' = OgenesOnlyDown,'offshore vs inshore control down' = origin_control_down)
venn(candidates) # 0 genes

# check for diffs in variance
nrow(in_genes_only[in_genes_only$log2FoldChange > 0 & in_genes_only$log2FoldChange > res_off_ingenes$log2FoldChange,]) # 42 genes where LC is higher in inshore
nrow(in_genes_only[in_genes_only$log2FoldChange > 0 & in_genes_only$log2FoldChange < res_off_ingenes$log2FoldChange,]) # 6 genes where LC is smaller in inshore
up=c(42,6)
chisq.test(up, p=c(0.5,0.5))







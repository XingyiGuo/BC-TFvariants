
library(data.table)
library(lme4)

setwd('/home/ec2-user/project1/BCTF/xingyiProc/BCACsnpAnnotation')

bcac <- fread("/home/ec2-user/project1/BCTF/bcac_matrix.txt",header=T)
names(bcac)

# bcac <- fread("/home/ec2-user/project/BCTF/bcac_matrix.txt",
#  select=c(1:24,32,37,38,44,56,57,58,59,60,70,71,72,73,74,75,76,77,78,79,82,83,
#           84,85,87,90,91,99,100,104,106,108,111,112,113,114,117,122,131),header=T)
#names(bcac)
dim(bcac) # 11792542, 137
#bcac<-data.frame(bcac)

names(bcac)[23]<-"AHR_ChIP_Seq_peaks"  

bcac <- within (bcac, {
ff  <-as.numeric(as.character(eaf))
ff.p<-as.numeric(as.character(eaf.p))
ff.n<-as.numeric(as.character(eaf.n))
bb  <-as.numeric(as.character(beta))
bb.p<-as.numeric(as.character(beta.p))
bb.n<-as.numeric(as.character(beta.n))
ss.p<-as.numeric(as.character(se.p))
ss.n<-as.numeric(as.character(se.n))
pp.p<-as.numeric(as.character(p.p))
pp.n<-as.numeric(as.character(p.n))
tv<-bb/se
tv.p<-bb.p/ss.p
tv.n<-bb.n/ss.n
})

bcac<-bcac[!duplicated(bcac[,3:4]),-c(7:8,11:20)]
names(bcac)
dim(bcac)  # 11712034 SNPs
bcac<-bcac[!is.na(bcac$tv) & bcac$chrZ>=1 & bcac$chrZ<=22,]
dim(bcac)  # 11337849 SNPs, 138 vars

bcac0<-bcac[bcac$tv==0,]
n0<-dim(bcac0)[1]
rm(bcac0)
bcac<-bcac[!is.na(bcac$tv) & bcac$chr>=1 & bcac$chr<=22,]
set.seed(123)
bcac$tv[bcac$tv==0]<-runif(n0,-0.001,0.001)

dim(bcac) # 11337849 SNPs, 61 vars
range(bcac$tv)  
range(bcac$tv.p,na.rm=T)    
range(bcac$tv.n,na.rm=T) #  -11.42478  13.26699

bcac<-bcac[order(bcac$chrZ,bcac$position_b37),]
#rm(CHR,gene_pos,gene_posx,i,bp.bp,bp.pp,KB,bp.x,bp.xa,bp.y,bp.ya)   

bcac$ID<-paste(bcac$chrZ,'-',bcac$position_b37,sep='')  
bcac$pc<-cut(bcac$p,breaks=c(0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,
                   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),label=1:15)
table(bcac$pc)

nx<-ceiling(10139060*c(1e-6,1e-5-1e-6,1e-4-1e-5,1e-3-1e-4,1e-2-1e-3,
                      1e-1-1e-2,rep(0.1,9)))
idv<-vector('list',15)
set.seed(123)
for (i in 1:15) {
bdx<-bcac[bcac$pc==i,] 
idv[[i]]<-sample(bdx$ID,nx[i])
}
table(duplicated(unlist(idv))) 

bcac$chosen<-ifelse(bcac$ID %in% unlist(idv),1,0)
rm(idv)
table(bcac$chosen) #1: 10139064 SNPs chosen for deflated genome, 0:  1198785 SNPs (enriched with highly significant SNPs)

gene_pos<-read.table('/home/ec2-user/project1/BCTF/Genecode.v27.gene.list.txt',head=F)
gene_pos<-gene_pos[!duplicated(gene_pos[,c(1,4,5,8)]) & 
                    gene_pos$V7=='protein_coding',]
dim(gene_pos)
head(gene_pos)
gene_pos<-gene_pos[gene_pos$V1 %in% 1:22 ,c(1,4,5,8)]
names(gene_pos)<-c('chr','a','b','gene')
head(gene_pos)
dim(gene_pos)

names(bcac)[3:4]<-c("chr","position")

KB<-100000
bcac$loci <- paste0(bcac$chr,'_',floor(bcac$position/KB))

length(unique(bcac$loci)) # 26864 loci

rm(bdx,i,KB,n0,nx) 

ccv_hpp<-read.csv('Table_S2_CCV_HPP.csv',head=T,stringsAsFactors=F)
head(ccv_hpp)
dim(ccv_hpp)
ccv_hpp<-ccv_hpp[!duplicated(ccv_hpp[,5:6]),]
head(ccv_hpp)
dim(ccv_hpp)
length(unique(ccv_hpp$Variant_Id)) 
length(unique(ccv_hpp$finemapping_region))  # 149 regions 
 
tab6<-read.csv('Table_S6_178genes.csv',stringsAsFactors=F,header=T)
dim(tab6)
head(tab6)
table(tab6$CCV)
length(unique(tab6$region)) #  84 regions
length(unique(tab6$gene))   # 178 genes
intersect(unique(ccv_hpp$finemapping_region),unique(tab6$region))
tab6$gene[!(tab6$gene %in% gene_pos$gene)]
ccv_hpp$ccv_region<-ifelse(ccv_hpp$finemapping_region %in% tab6$region,1,0)
head(ccv_hpp)

intersect(names(bcac),names(ccv_hpp))

bcac<-merge(bcac,ccv_hpp,all=T)
bcac$snpID<-paste('S_',bcac$chr,'.',bcac$position,sep='')
dim(bcac)  
bcac<-bcac[!is.na(bcac$tv),]

risk_SNPs<-read.csv('risk_SNPs_from_WuLong.csv',stringsAsFactors=F,header=T)
dim(risk_SNPs)
head(risk_SNPs)

rpt_snps<-read.table('/home/ec2-user/project1/BCTF/GWAS_risk_SNPs.txt',head=T)
names(rpt_snps)[2:3]<-c('CHR','POS')


anno<-fread('BCACsnp.annotated.snp',
             select=c(1,2,103,114,116,118,126,136,145,150,210,290,291), head=F)
dim(anno) # 10446448, 13
anno<-as.data.frame(anno[!duplicated(anno[,1:2]),])
dim(anno) # 10429719, 13
names(anno)<-c('chr','position',"RegulomeDB_score","CADD_phred","DANN_score",
 "MKL_non_coding_score","XF_noncoding_score","ENCODE_Dnase_score",
 "Enhancer_permissive","CAGE_peak_permissive",
 "M_CAP_score","E028_15_coreMarks","E027_15_coreMarks")   

head(anno)

anno<-within(anno, {
RegulomeDB_score<-as.numeric(substr(RegulomeDB_score,1,1))
CADD_phred<-as.numeric(CADD_phred)
DANN_score<-as.numeric(DANN_score)
MKL_non_coding_score<-as.numeric(MKL_non_coding_score)
XF_noncoding_score<-as.numeric(XF_noncoding_score)
ENCODE_Dnase_score<-as.numeric(ENCODE_Dnase_score)  # too many missing
M_CAP_score<-as.numeric(M_CAP_score) # too many missing
E028_15_coreMarksn<-ifelse(E028_15_coreMarks %in% c('Active TSS','Bivalent/Poised TSS',
                   'Flanking Active TSS','Flanking Bivalent TSS/Enh'), 1, 
  ifelse(E028_15_coreMarks %in% c('Bivalent Enhancer','Enhancers','Genic enhancers'),2, 
  ifelse(E028_15_coreMarks %in% c('Strong transcription','Weak transcription'),3,
  ifelse(E028_15_coreMarks %in% c('Repressed PolyComb','Weak Repressed PolyComb'),4,0))))
E027_15_coreMarksn<-ifelse(E027_15_coreMarks %in% c('Active TSS','Bivalent/Poised TSS',
                   'Flanking Active TSS','Flanking Bivalent TSS/Enh'), 1, 
  ifelse(E027_15_coreMarks %in% c('Bivalent Enhancer','Enhancers','Genic enhancers'),2, 
  ifelse(E027_15_coreMarks %in% c('Strong transcription','Weak transcription'),3,
  ifelse(E027_15_coreMarks %in% c('Repressed PolyComb','Weak Repressed PolyComb'),4,0))))
 } )
anno<-anno[anno$chr %in% 1:22,]
anno$chr<-as.integer(anno$chr)
intersect(names(bcac),names(anno))

bcac<-merge(bcac,anno,all=T)
bcac$Pcc<-ifelse(bcac$p<5e-8,1,0)
bcac$Pcc.p<-ifelse(bcac$pp.p<5e-8,1,0)
bcac$Pcc.n<-ifelse(bcac$pp.n<5e-8,1,0)

dim(bcac)  # 11337856
bcac<-bcac[!is.na(bcac$tv),]
table(is.na(bcac$tv)) 
dim(bcac)  # 11337849, 99 vars
table(bcac$Pcc)
table(bcac$chosen)
# table(bcac$Pcc.p)
# table(bcac$Pcc.n)

####  Association of single TF with Chi-square for BC risk
names(bcac)[c(12,14,20,25,26,39,47,57,62,65,72,75,78,87,88,92:96,98:100,105,110,111,117,119,121,125)]

outr1<-apply(bcac[,c(12,14,20,25,26,39,47,57,62,65,72,75,78,87,88,92:96,98:100,105,110,111,117,119,121,125)],2,
   function (x) summary(lmer(I(bcac$tv^2)~x+(1|bcac$loci),control = lmerControl(calc.derivs = FALSE)))$coef[2,])
t(outr1)
  
bcacx<-bcac[bcac$chosen==1,]	
outr2<-apply(bcacx[,c(12,14,20,25,26,39,47,57,62,65,72,75,78,87,88,92:96,98:100,105,110,111,117,119,121,125)],2,
   function (x) summary(lmer(I(bcacx$tv^2)~x+(1|bcacx$loci)))$coef[2,])
t(outr2)

outr3<-apply(bcacx[,11:125],2,function (x)  table(x))

#### Number and percentage of TF-occupied SNPs
gwas_loci<-read.csv('GWAS_loci.csv',head=T)
head(gwas_loci) # 169 loci
out<-matrix(NA,dim(gwas_loci)[1],63)
for (i in 1:dim(gwas_loci)[1]) {
#i<-1
bcacx<-bcac[bcac$chr==gwas_loci[i,1] & bcac$position>=gwas_loci[i,2] & bcac$position<=gwas_loci[i,3],]
out[i,]<-c(i,table(bcacx$Pcc)[1:2],
table(bcacx$AR.16h_peaks ,bcacx$Pcc)[c(2,4)],                             
table(bcacx$BRD4.E2.rep1_peaks  ,bcacx$Pcc)[c(2,4)],                      
table(bcacx$Cebpbsc150V0422111_peaks ,bcacx$Pcc)[c(2,4)],                 
table(bcacx$COT2_REP1_peaks ,bcacx$Pcc)[c(2,4)],                          
table(bcacx$CtBP_ChIPSeq_peaks,bcacx$Pcc)[c(2,4)],                        
table(bcacx$ER_REP2_peaks ,bcacx$Pcc)[c(2,4)],                            
table(bcacx$E2F1_peaks,bcacx$Pcc)[c(2,4)],                                
table(bcacx$Fosl2V0422111_peaks,bcacx$Pcc)[c(2,4)],                       
table(bcacx$FoxA1_T47D_fullMedia_peaks ,bcacx$Pcc)[c(2,4)],              
table(bcacx$Foxm1sc502V0422111_peaks,bcacx$Pcc)[c(2,4)],                  
table(bcacx$GabpV0422111_peaks,bcacx$Pcc)[c(2,4)],                        
table(bcacx$Gata3V0422111_peaks,bcacx$Pcc)[c(2,4)],                       
table(bcacx$GREB1_REP1_peaks,bcacx$Pcc)[c(2,4)],                          
table(bcacx$Hae2f1Ucd_peaks,bcacx$Pcc)[c(2,4)],                           
table(bcacx$Hdac2sc6296V0422111_peaks,bcacx$Pcc)[c(2,4)],                 
table(bcacx$JundV0422111_peaks,bcacx$Pcc)[c(2,4)],                        
table(bcacx$KLF4_REP1_peaks,bcacx$Pcc)[c(2,4)],                           
table(bcacx$MaxV0422111_peaks ,bcacx$Pcc)[c(2,4)],                        
table(bcacx$MYC.6h_peaks ,bcacx$Pcc)[c(2,4)],                             
table(bcacx$Nr2f2sc271940V0422111_peaks,bcacx$Pcc)[c(2,4)],              
table(bcacx$P300V0422111_peaks,bcacx$Pcc)[c(2,4)],                        
table(bcacx$Pmlsc71910V0422111_peaks,bcacx$Pcc)[c(2,4)],                 
table(bcacx$Pol2_peaks,bcacx$Pcc)[c(2,4)],                                
table(bcacx$Sin3ak20V0422111_peaks ,bcacx$Pcc)[c(2,4)],                   
table(bcacx$SrfV0422111_peaks,bcacx$Pcc)[c(2,4)],                         
table(bcacx$Taf1V0422111_peaks ,bcacx$Pcc)[c(2,4)],                       
table(bcacx$Tcf12V0422111_peaks,bcacx$Pcc)[c(2,4)],                       
table(bcacx$Tcf7l2Ucd_peaks ,bcacx$Pcc)[c(2,4)],                          
table(bcacx$TLE3_REP1_peaks,bcacx$Pcc)[c(2,4)],               
table(bcacx$Znf217Ucd_peaks  ,bcacx$Pcc)[c(2,4)])                       
if (i %% 5 ==1)  cat(i,"/",dim(gwas_loci)[1],"\n")
 }
head(out) 
write.csv(out,'out_freq2.csv')

#### Association of co-occupancy of two TFs with Chi-square for breast cancer risk
ww<-c(12,20,25,26,39,47,62,105)
N0<-8
N1<-2
outx<-matrix(NA,choose(N0,N1),16)
ij<-0
for (i in 1:N0) {
 for (i in 1:(N0-1)) {
 for (j in (i+1):N0) {
 # i<-1
 # j<-i+1
 wi<-ww[i]
 wj<-ww[j]
 ij<-ij+1
  co_occ<-ifelse(bcac[,wi]==1 & bcac[,wj]==1,1,ifelse(bcac[,wi]==1 & bcac[,wj]==0,2,
            ifelse(bcac[,wi]==0 & bcac[,wj]==1,3,0)))  
  tab<-table(co_occ)
  m1<-lmer(I(bcac$tv^2)~bcac[,wi]+bcac[,wj]+(1|bcac$loci),control = lmerControl(calc.derivs = FALSE),REML=F)
  m2<-lmer(I(bcac$tv^2)~as.factor(co_occ)+(1|bcac$loci),control = lmerControl(calc.derivs = FALSE),REML=F) 
  cf<-summary(m2)$coef 
  pint<-anova(m1,m2)
  outx[ij,]<-c(cf[2,1:3],cf[3,1:3],cf[4,1:3],pint$Pr[2],tab,i,j)
 }
  #}
colnames(outx)<-c('beta1','se1','t1','beta2','se2','t3','beta3','se3','t3','p_int','n0','n1','n2','n3','i','j')
outx


# deflated genome
bcacx<-bcac[bcac$chosen==1,]
dim(bcacx)
ww<-c(12,20,25,26,39,47,62,75,105)
N0<-9
outx1<-matrix(NA,9,14)
ij<-0
for (i in 1:N0) {
  wi<-ww[i]
  wj<-118  # TCF7
 ij<-ij+1
   co_occ<-ifelse(bcacx[,wi]==1 & bcacx[,wj]==1,1,ifelse(bcacx[,wi]==1 & bcacx[,wj]==0,2,
             ifelse(bcacx[,wi]==0 & bcacx[,wj]==1,3,0)))  
   tab<-table(co_occ)
   m1<-lmer(I(bcacx$tv^2)~bcacx[,wi]+bcacx[,wj]+(1|bcacx$loci),control = lmerControl(calc.derivs = FALSE),REML=F)
   m2<-lmer(I(bcacx$tv^2)~as.factor(co_occ)+(1|bcacx$loci),control = lmerControl(calc.derivs = FALSE),REML=F) 
   cf<-summary(m2)$coef 
   pint<-anova(m1,m2)
   outx1[ij,]<-c(cf[2,1:3],cf[3,1:3],cf[4,1:3],pint$Pr[2],tab)
  }
outx1


###Association of co-occupancy of TFs with breast cancer risk, stratified by FOXA1
# for ESR1/E2F1, ESR1/TCF12V, TCF12V/TLE3, SIN3AK/TLE3

table(bcac$FoxA1_T47D_fullMedia_peaks, bcac$TF_coocc)

mm1<- lmer(I(tv^2)~FoxA1_T47D_fullMedia_peaks*TF_coocc+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F,data=bcac)
mm2<- lmer(I(tv^2)~FoxA1_T47D_fullMedia_peaks+TF_coocc+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F,data=bcac)
anova(mm1,mm2) #interacton test
 # stratified by FOXA1 
summary(lmer(I(tv^2)~TF_coocc+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F, data=bcac[bcac$FoxA1_T47D_fullMedia_peaks==0,]))
summary(lmer(I(tv^2)~TF_coocc+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F, data=bcac[bcac$FoxA1_T47D_fullMedia_peaks==1,]))


### TF score from 22 selected TFs
bcac$TF_22score<-
bcac$AR.16h_peaks   +                           
bcac$Cebpbsc150V0422111_peaks+                  
bcac$CtBP_ChIPSeq_peaks      +                  
bcac$E2F1_peaks              +                  
bcac$ER_REP2_peaks           +                  
bcac$Fosl2V0422111_peaks     +                  
bcac$FoxA1_T47D_fullMedia_peaks+                
bcac$Foxm1sc502V0422111_peaks+                  
bcac$Gata3V0422111_peaks     +                  
bcac$GREB1_REP1_peaks        +                  
bcac$Hae2f1Ucd_peaks         +                  
bcac$Hdac2sc6296V0422111_peaks+                 
bcac$JundV0422111_peaks      +                  
bcac$Nr2f2sc271940V0422111_peaks+               
bcac$P300V0422111_peaks      +                  
bcac$Pmlsc71910V0422111_peaks+                  
bcac$Sin3ak20V0422111_peaks  +                  
bcac$SrfV0422111_peaks       +                  
bcac$Tcf12V0422111_peaks     +                  
bcac$Tcf7l2Ucd_peaks         +                  
bcac$TLE3_REP1_peaks         +                  
bcac$Znf217Ucd_peaks        ;   

 # remove FOXA1
bcac$TF_21score<-
bcac$AR.16h_peaks   +                           
bcac$Cebpbsc150V0422111_peaks+                  
bcac$CtBP_ChIPSeq_peaks      +                  
bcac$E2F1_peaks              +                  
bcac$ER_REP2_peaks           +                  
bcac$Fosl2V0422111_peaks     +                  
bcac$Foxm1sc502V0422111_peaks+                  
bcac$Gata3V0422111_peaks     +                  
bcac$GREB1_REP1_peaks        +                  
bcac$Hae2f1Ucd_peaks         +                  
bcac$Hdac2sc6296V0422111_peaks+                 
bcac$JundV0422111_peaks      +                  
bcac$Nr2f2sc271940V0422111_peaks+               
bcac$P300V0422111_peaks      +                  
bcac$Pmlsc71910V0422111_peaks+                  
bcac$Sin3ak20V0422111_peaks  +                  
bcac$SrfV0422111_peaks       +                  
bcac$Tcf12V0422111_peaks     +                  
bcac$Tcf7l2Ucd_peaks         +                  
bcac$TLE3_REP1_peaks         +                  
bcac$Znf217Ucd_peaks        ;   

bcac$TF_22scorec <- as.factor(ifelse (bcac$TF_22score==0,0,
                    ifelse (bcac$TF_22score<=2,1,
                    ifelse (bcac$TF_22score<=5,2,
                    ifelse (bcac$TF_22score<=10,3,4)))))
bcac$TF_21scorec <- as.factor(ifelse (bcac$TF_21score==0,0,
                    ifelse (bcac$TF_21score<=2,1,
                    ifelse (bcac$TF_21score<=5,2,
                    ifelse (bcac$TF_21score<=10,3,4)))))
               
			   
summary(lmer(I(tv^2)~ TF_22scorec +(1|loci),data=bcac,control = lmerControl(calc.derivs = FALSE))

mm1<- lmer(I(tv^2)~FoxA1_T47D_fullMedia_peaks*TF_21scorec+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F,data=bcac)
mm2<- lmer(I(tv^2)~FoxA1_T47D_fullMedia_peaks+TF_21scorec+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F,data=bcac)
anova(mm1,mm2) #interacton test
 # stratified by FOXA1 
summary(lmer(I(tv^2)~TF_21scorec+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F, data=bcac[bcac$FoxA1_T47D_fullMedia_peaks==0,]))
summary(lmer(I(tv^2)~TF_21scorec+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F, data=bcac[bcac$FoxA1_T47D_fullMedia_peaks==1,]))

### Association of chromatin states and stratified by TF score
bcac$TF_scorec3 <- as.factor(ifelse (bcac$TF_22score==0,0,
                    ifelse (bcac$TF_22score<=5,1,2)))

table(bcac$E027_15_coreMarks)
table(bcac$E028_15_coreMarks)
summary(lmer(I(tv^2)~E028_15_coreMarks+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F, data=bcac))

mm1<- lmer(I(tv^2)~E028_15_coreMarks*TF_scorec3+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F,data=bcac)
mm2<- lmer(I(tv^2)~E028_15_coreMarks+TF_scorec3+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F,data=bcac)
anova(mm1,mm2) #interacton test

summary(lmer(I(tv^2)~E028_15_coreMarks+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F, data=bcac[bcac$TF_scorec3==0,]))
summary(lmer(I(tv^2)~E028_15_coreMarks+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F, data=bcac[bcac$TF_scorec3==1,]))
summary(lmer(I(tv^2)~E028_15_coreMarks+(1|loci),control = lmerControl(calc.derivs = FALSE),REML=F, data=bcac[bcac$TF_scorec3==2,]))


### Genetic variations of TF-DNA binding of TFs, stratified by motifs
motifs<-read.table('motifs.txt',head=T)
dim(motifs)
names(motifs)<-c('chr','position','TCF7L2x','AHRs','ARs','CEBPBs','CTCFVEHs','E2F1s','ERs','FOSL2s','FOXA1s',
  'FOXM1s','GABPs','GATA3s','JUNDs','KLF4s','MAXs','MYCs','COT2s','SRFs','TCF12s','TCF7L2s')  
head(motifs)
motifs<-motifs[!duplicated(motifs[,1:2]),-c(3,4,7,13,16,17,18)] 
dim(motifs)  # 11712034 15
intersect(names(bcac),names(motifs))
bcac<-merge(bcac,motifs)
dim(bcac)  # 11337849      194

bcac$ARs[bcac$ARs==-1]<-200
bcac$CEBPBs[bcac$CEBPBs==-1]<-200
bcac$E2F1s[bcac$E2F1s==-1]<-200
bcac$ERs[bcac$ERs==-1]<-200
bcac$FOSL2s[bcac$FOSL2s==-1]<-200
bcac$FOXA1s[bcac$FOXA1s==-1]<-200
bcac$FOXM1s[bcac$FOXM1s==-1]<-200
bcac$GATA3s[bcac$GATA3s==-1]<-200
bcac$JUNDs[bcac$JUNDs==-1]<-200
bcac$COT2s[bcac$COT2s==-1]<-200
bcac$SRFs[bcac$SRFs==-1]<-200
bcac$TCF12s[bcac$TCF12s==-1]<-200
bcac$TCF7L2s[bcac$TCF7L2s==-1]<-200

bcac$Bmotifs<-pmin(bcac$ARs,bcac$CEBPBs,bcac$E2F1s,bcac$ERs,bcac$FOSL2s,
         bcac$FOXA1s,bcac$FOXM1s,bcac$GATA3s,bcac$JUNDs,bcac$COT2s,
		 bcac$SRFs,bcac$TCF12s,bcac$TCF7L2s,na.rm=T)
 

ww1<-c(12,20,39,47,57,62,65,75,92,25,110,117,119)
names(bcacx)[ww1]
ww2<-182:194
names(bcacx)[ww2]

out<-matrix(NA,13,42)
for (i in 1:13) {
i1<-ww1[i]
i2<-ww2[i]
Apeaks<-bcacx[,i1]
Amotifs<-bcacx[,i2]
tmp<-ifelse(Apeaks==1 & Amotifs<200,1,
     ifelse(Apeaks==1 & bcacx$Bmotifs<200,2,
     ifelse(Apeaks==1,3,0)))
tmpx<-ifelse(Apeaks==1 & Amotifs==0,1,
      ifelse(Apeaks==1 & Amotifs>0 & Amotifs<=20,2,
      ifelse(Apeaks==1 & Amotifs>20 & Amotifs<200,3,
      ifelse(Apeaks==1 & bcacx$Bmotifs==0,4,
      ifelse(Apeaks==1 & bcacx$Bmotifs>0 & bcacx$Bmotifs<=20,5,
      ifelse(Apeaks==1 & bcacx$Bmotifs>20 & bcacx$Bmotifs<200,6,
      ifelse(Apeaks==1,7,0)))))))
tab<-table(tmp)
tabx<-table(tmpx)
system.time(mod1a<-lmer(I(tv^2)~ as.factor(tmp)+(1|loci),data=bcac,REML=F,
 control = lmerControl(calc.derivs = FALSE)))  
system.time(mod1b<-lmer(I(tv^2)~ as.factor(tmpx)+(1|loci),data=bcac,REML=F,
  control = lmerControl(calc.derivs = FALSE)))  
cfa<-summary(mod1a)$coef
cfb<-summary(mod1b)$coef
out[i,]<-c(tab,cfa[2,],cfa[3,],cfa[4,],tabx,cfb[2,],cfb[3,],cfb[4,],cfb[5,],cfb[6,],cfb[7,],cfb[8,])
 }
colnames(out)<-c('N0','N1','N2','N3','beta1','se1','t1','beta2','se2','t2','beta3','se3','t3',
                 'Nx0','Nx1','Nx2','Nx3','Nx4','Nx5','Nx6','Nx7','betax1','sex1','tx1','betax2','sex2','tx2',
				 'betax3','sex3','tx3','betax4','sex4','tx4','betax5','sex5','tx5','betax6','sex6','tx6','betax7','sex7','tx7')
rownames(out)<-names(bcacx)[ww1]


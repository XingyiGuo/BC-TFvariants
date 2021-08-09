
###############################################
####  Computer system requirements         ####
####  R package:  data.table, lme4         ####
####  Unix operating system environment    ####
####  15 GB memory                         #### 
###############################################

library(data.table) # read in large data matrix fast
library(lme4)       # for mixed model
library(parallel)

setwd('/home/ec2-user/project1/chenzs/crctf/TF_only')

# Input TF matrix
```
### example data matrix file see "CRC_GWAS_TFs_matrix_final_include.matrix" under the same directory. In pratice, you need to adjust number of TFs in the mastrix. 
crc <- fread("CRC_GWAS_TFs_matrix_final_include.matrix",header=T)
names(crc)
dim(crc)
 ```
# cauculate t values for GWAS associations
```
crc <- within (crc, {
freq  <-as.numeric(as.character(Freq1))
beta  <-as.numeric(as.character(Effect))
se  <-as.numeric(as.character(StdErr))
tv<-beta/se
 })
```
 # remove duplicates and those with missing t-values
 ```
crc <- crc[!duplicated(crc[,236:237]),-c(5:7,11:15)]
crc<-crc[!is.na(crc$tv) & crc$chrZ>=1 & crc$chrZ<=22,]
dim(crc)  
```
 # replace t-value ==0 with random values from -0.001 to 0.001 
```
crc0<-crc[crc$tv==0,]
n0<-dim(crc0)[1]
rm(crc0)
set.seed(123)
crc$tv[crc$tv==0]<-runif(n0,-0.001,0.001)

range(crc$tv)  
```
##### create a random genome based on the random distribution of t-values for GWAS association #####
  # sort data by chromosome and positions
  ```
crc<-crc[order(crc$chrZ,crc$pos),]
  # create SNPs ID using Chr-position
crc$ID<-paste(crc$chrZ,'-',crc$pos,sep='')  
crc$pc<-cut(crc$p,breaks=c(0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,
                   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),label=1:15)
table(crc$pc) # the smallest freq of 15 crc$pc categories = 815915
 # estimate the expected number of SNPs in each p-values category
nx<-ceiling(10*815915*c(1e-6,1e-5-1e-6,1e-4-1e-5,1e-3-1e-4,1e-2-1e-3,1e-1-1e-2,rep(0.1,9)))
idv<-vector('list',15)
set.seed(123)
for (i in 1:15) {
bdx<-crc[crc$pc==i,] 
idv[[i]]<-sample(bdx$ID,nx[i])
}
table(duplicated(unlist(idv))) 
crc$chosen<-ifelse(crc$ID %in% unlist(idv),1,0)
rm(idv)
table(crc$chosen) #1: 8159194 SNPs chosen for deflated genome, 0: 952306 SNPs (enriched with highly significant SNPs)
```
# input gene position data
```
gene_pos<-read.table('/home/ec2-user/project1/chenzs/crctf/TF_only/gencode.v19.annotation.gtf.Gene',head=F)
gene_pos$V3<-gsub("chr(\\w+)","\\1",gene_pos$V3, perl=T)
gene_pos <- gene_pos[!duplicated(gene_pos[,c(3,4,5,2)]) & gene_pos$V7=='protein_coding',]
gene_pos<-gene_pos[gene_pos$V3 %in% 1:22 ,c(3,4,5,2)]
names(gene_pos)<-c('chr','a','b','gene')
dim(gene_pos)
```
#### define LD block based on 100kb
```
names(crc)[228,229]<-c("chr","position")
KB<-100000
crc$loci <- paste0(crc$chr,'_',floor(crc$position/KB))

length(unique(crc$loci)) # 26578 loci
rm(bdx,i,KB,n0,nx) 

# crc$Pcc<-ifelse(crc$p<5e-8,1,0) # define a binary outcome variable for GWAS fignificance, which have lower statistical power than models with continuous Chi-square as the outome variable  

crc <- as.data.frame (crc)
crc <- crc[, - which(colnames(crc) %in% "TF_219")] # remove this TF due to error
```
####  Association of single TF with continuous Chi-square (==tv^2) for CRC risk using genealized mixed models
# computing time for each model: about 5 minutes
```
names(crc)[9:227]
outr1<-mclapply(crc[,9:227], function (x) {summary(lmer(I(crc$tv^2)~x+(1|crc$loci),control = lmerControl(calc.derivs = FALSE)))$coef[2,]}, mc.cores = 6)

# check errors in parallel
# bad <- sapply (outr1, inherits, what = "try-error")
# a <- names(outr1[bad])
# a <- a[-1] # remove TF_219
# for (i in 1:length(a)) {
#  df <- crc[,a[i]]
#  outr1 [[a[i]]] <- summary(lmer(I(crcx$tv^2)~df+(1|crcx$loci), control = lmerControl(calc.derivs = FALSE)))$coef[2,]
#  }

# output results - 
output <- data.frame(matrix(unlist(outr1), ncol=3, byrow=TRUE)
#rownames (output) <- names(crc)[9:227][-grep("TF_219", names(crc)[9:227])]
rownames (output) <- names(crc)[9:226]
colnames (output) <- c('bata','se','tv')
write.table(output,'output.txt', row.names=TRUE, sep="\t", quote=FALSE)
```
# Analysis for the deflated genome
```
crcx <- crc [crc$chosen==1,]
outr2<-mclapply(crcx[,9:226], function (x) {summary(lmer(I(crcx$tv^2)~x+(1|crcx$loci),control = lmerControl(calc.derivs = FALSE)))$coef[2,]}, mc.cores = 6)
```
# Output results
```
outputx <- data.frame(matrix(unlist(outr2), ncol=3, byrow=TRUE)
rownames (outputx) <- names(crcx)[9:226]
colnames (outputx) <- c('bata','se','tv')
write.table(outputx,'outputx.txt', row.names=TRUE, sep="\t", quote=FALSE)
```
# TF frequency
```
outfreq <- apply(crc[,9:226],2,function (x) table(x))
outfreqx <- apply(crcx[,9:226],2,function (x) table(x))
write.table(t(outfreq),'outfreq.txt', row.names=TRUE, sep="\t", quote=FALSE)
write.table(t(outfreqx),'outfreqx.txt', row.names=TRUE, sep="\t", quote=FALSE)
```

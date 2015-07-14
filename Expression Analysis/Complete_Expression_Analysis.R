# Script to analyze expression data
# Sadhika Malladi, Jonathan Ma, Andrew Beck

# The following script includes:
# 1. Expression and clinical data consolidation
# 2. Wilcoxon Rank-Sum test with respect to sex in each tissue
# 3. Benjamini-Hochberg method to adjust p-values
# 4. Fisher's Method to combine cancer-specific p-values
# 5. Fold Change Computation
# 6. GWAS-Like Analysis
# 7. Gene Set Analysis


# Set your working directory to wherever the RData file is stored
setwd('~/googledrive/ExpressionSexProject/Expression_Signatures/')

load("TCGA.Clin.RNAseq.v2.RData")

# Load relevant libraries
library(gdata)
library(MADAM)
# Uncomment to install GWASTools package: 
# source("http://bioconductor.org/biocLite.R")
# biocLite("GWASTools")
library(GWASTools)
library(gap)
library(NCBI2R)
library(GSA)


# Create a data frame with clinical information ---------------------------
# Includes several overall graphs to describe trends
expression.data <- dat
clinical.data <- clin

expression.data <- t(expression.data)

# Create a data frame with tissue type, barcode, and sample type (tumor or normal)
expression.dataframe <- data.frame(t(data.frame(strsplit(rownames(expression.data),"-|\\."),stringsAsFactors=F)))
expression.dataframe$Barcode <- apply(expression.dataframe,MARGIN=1,function(x)paste(x[c(2:4)],collapse="-"))
expression.dataframe <- expression.dataframe[,c("X1","X5","Barcode")]
colnames(expression.dataframe) <- c("TissueType","SampleType","Barcode")
expression.dataframe$SampleType=substr(expression.dataframe$SampleType,1,2)

# Patterns for tumor and normal sample types
# Normal: Number 11 followed by any capital letter
# Tumor: Number 0 followed by 1 or 3 followed by any capital letter
Tumor.expression <- t(expression.data[ expression.dataframe$SampleType=="01" | expression.dataframe$SampleType=="03",])
Normal.expression <- t(expression.data[ expression.dataframe$SampleType=="11",])
Tumor.dataframe <- expression.dataframe[ expression.dataframe$SampleType=="01"|expression.dataframe$SampleType=="03",]
Normal.dataframe <- expression.dataframe[ expression.dataframe$SampleType=="11",]

# Rownames of clinical data are barcodes. 
# To access rows of clinical data, you need a vector of barcodes.
rownames(clinical.data) <- toupper(gsub("[A-Z]+.(tcga).([a-z0-9]+).([a-z0-9]+)","\\1-\\2-\\3",rownames(clinical.data)))
colnames(Tumor.expression) <- Tumor.dataframe$Barcode
colnames(Normal.expression) <- Normal.dataframe$Barcode

# We only want to do computations on normals and tumors for which we have clinical data 
# (in this case, gender)
Normal.dataframe <- Normal.dataframe[Normal.dataframe$Barcode%in% rownames(clinical.data),]
Normal.expression <- Normal.expression[ ,Normal.dataframe$Barcode]
Tumor.dataframe <- Tumor.dataframe[ Tumor.dataframe$Barcode %in% rownames(clinical.data),]
Tumor.expression <- Tumor.expression[ ,Tumor.dataframe$Barcode]

# Include age for future analyses
Normal.dataframe$yearstobirth <- clinical.data[ Normal.dataframe$Barcode,"yearstobirth"]
Normal.dataframe$sex <- clinical.data[ Normal.dataframe$Barcode,"gender"]
Tumor.dataframe$yearstobirth <- clinical.data[ Tumor.dataframe$Barcode,"yearstobirth"]
Tumor.dataframe$sex <- clinical.data[ Tumor.dataframe$Barcode,"gender"]

# Note: sort returns reordered values, order returns reordered indices
# Order dataframe by ordering the ages and returning reordered indices
Tumor.dataframe <- Tumor.dataframe[ order(Tumor.dataframe$yearstobirth),]
# Plotting age, colored by tissue type -- not used in later analyses, just for fun
plot(sort(Tumor.dataframe$yearstobirth),col=as.factor(Tumor.dataframe$TissueType), main="Age colored by tissue type",xlab="Tissue Types (as factors)",ylab="Age")
# Age distribution by tissue type
boxplot(as.numeric(Tumor.dataframe$yearstobirth)~Tumor.dataframe$TissueType, main="Age vs Tissue Type (Tumor)")
boxplot(as.numeric(Normal.dataframe$yearstobirth)~Normal.dataframe$TissueType, main="Age vs Tissue Type (Normal)")

Tumor.TypeBySex <- table(Tumor.dataframe$TissueType,Tumor.dataframe$sex)
Normal.TypeBySex <- table(Normal.dataframe$TissueType,Normal.dataframe$sex)

# Subset tumor and normal dataframes so that they only include the mixed tissue types
# Mixed tissue types have more than 10 cases in male and more than 10 cases in female
Tumor.dataframe <- Tumor.dataframe[ Tumor.dataframe$TissueType %in% rownames(Tumor.TypeBySex)[apply(Tumor.TypeBySex,1,function(x)sum(x>10)==2)],]
Normal.dataframe <- Normal.dataframe[ Normal.dataframe$TissueType %in% rownames(Normal.TypeBySex)[apply(Normal.TypeBySex,1,function(x)sum(x>10)==2)],]
Tumor.expression <- Tumor.expression[,Tumor.dataframe$Barcode]
Normal.expression <- Normal.expression[,Normal.dataframe$Barcode]

# Removes cases that have "0" (i.e. less than 10) males or females 
Tumor.dataframe$TissueType=drop.levels(Tumor.dataframe$TissueType)
# Plot tissue type by sex
tissueBySex<-table(Tumor.dataframe$TissueType,Tumor.dataframe$sex)
tissueBySex<-as.table(tissueBySex[order(apply(tissueBySex,1,sum)),])
plot(tissueBySex,col=c("hotpink","lightblue"),las=3,main="Gender Distribution of 404 Normal Samples of 10 Types",ylab="",xlab="")

# p value matrices: no data is inputted here, only dimensions are set
# row = genes, col = tissue type (unique)
Tumor.WilcP <- matrix(NA,nrow=nrow(Tumor.expression),ncol=length(unique(Tumor.dataframe$TissueType)),dimnames=list(rownames(Tumor.expression),unique(Tumor.dataframe$TissueType)))
Normal.WilcP <- matrix(NA,nrow=nrow(Normal.expression),ncol=length(unique(Normal.dataframe$TissueType)),dimnames=list(rownames(Normal.expression),unique(Normal.dataframe$TissueType)))
Tumor.WilcEst <- Tumor.WilcP
Normal.WilcEst <- Normal.WilcP

# Remember: these will be saved to your working directory!
saveRDS(Tumor.expression,file="TumorExpression.RDS")
saveRDS(Normal.expression,file="NormalExpression.RDS")
saveRDS(Tumor.dataframe,file="TumorDataframe.RDS")
saveRDS(Normal.dataframe,file="NormalDataframe.RDS")


# Calculate p-values ------------------------------------------------------
for(i in colnames(Tumor.WilcP)) {
  # Get all barcodes of the current tissue type
  Tumor.tissueCases<-Tumor.dataframe$Barcode[Tumor.dataframe$TissueType==i]
  # Expression values of all genes for current tissue type (genes=rows, barcodes=columns)
  Tumor.tissueExpression<-Tumor.expression[,Tumor.tissueCases]
  # A list of genders of all barcodes (cases) of the current tissue type
  Tumor.tissueSex<-Tumor.dataframe$sex[Tumor.dataframe$TissueType==i]
  
  # Compute the wilcoxon rank-sum test with respect to sex using the expression for the
  # current tissue.
  
  # In this step, we are saving the median differences
  Tumor.WilcEst[,i]=apply(Tumor.tissueExpression,1,function(x)(
    if(sum(x>0)>(length(x)/10)){ 
      wilcox.test(x~Tumor.tissueSex,conf.int=T)$est
    }else{
      NA
    })
  )
  # Extract the p-value
  Tumor.WilcP[,i]=apply(Tumor.tissueExpression,1,function(x)((wilcox.test(x~Tumor.tissueSex)$p.value)))
  
  cat("Tissue ", i, " Tumor")
  cat("\n")
}
saveRDS(Tumor.WilcEst,"MedianDiffEstimate_Tumors_6.23.RDS")
saveRDS(Tumor.WilcP,"WilcoxonRawP_Tumors_6.23.RDS")
 
# Same procedure as above, but for normal tissues
 for(i in colnames(Normal.WilcP)) {
   Normal.tissueCases<-Normal.dataframe$Barcode[Normal.dataframe$TissueType==i]
   Normal.tissueExpression<-Normal.expression[,Normal.tissueCases]
   Normal.tissueSex<-Normal.dataframe$sex[Normal.dataframe$TissueType==i]
   Normal.WilcEst[,i]=apply(Normal.tissueExpression,1,function(x)(
     if(sum(x>0)>(length(x)/10)){ 
       wilcox.test(x~Normal.tissueSex,conf.int=T)$est
     }else{
       NA
     })
   )
   Normal.WilcP[,i]=apply(Normal.tissueExpression,1,function(x)((wilcox.test(x~Normal.tissueSex)$p.value)))
   
   cat("Tissue ", i, " Normal")
   cat("\n")
 }
 saveRDS(Normal.WilcEst,"MedianDiffEstimate_Normals_6.23.RDS")
 saveRDS(Normal.WilcP,"WilcoxonRawP_Normals_6.23.RDS")

# Adjust p-values ---------------------------------------------------------
 
# Manually adjust p values in order to find the most significant genes in each cancer type separately
adj.Tumor.WilcP.vector=c(Tumor.WilcP)
adj.Tumor.WilcP.vector=p.adjust(adj.Tumor.WilcP.vector,method='fdr')
adj.Tumor.WilcP=matrix(adj.Tumor.WilcP.vector,ncol=ncol(Tumor.WilcP),byrow=F)
rownames(adj.Tumor.WilcP)=rownames(Tumor.WilcP)
colnames(adj.Tumor.WilcP)=colnames(Tumor.WilcP)

adj.Normal.WilcP.vector=c(Normal.WilcP)
adj.Normal.WilcP.vector=p.adjust(adj.Normal.WilcP.vector,method='fdr')
adj.Normal.WilcP=matrix(adj.Normal.WilcP.vector,ncol=ncol(Normal.WilcP),byrow=F)
rownames(adj.Normal.WilcP)=rownames(Normal.WilcP)
colnames(adj.Normal.WilcP)=colnames(Normal.WilcP)

saveRDS(adj.Tumor.WilcP,file='WilcoxonAdjustedP_Tumors_6.24.RDS')
saveRDS(adj.Normal.WilcP,file='WilcoxonAdjustedP_Normals_6.24.RDS')


# Fisher's Method -----------------------------------------------------------

# Combine the p values generated by the wilcoxon test to form pan-cancer p values
Tumor.MetaPs=fisher.method(Tumor.WilcP,na.rm=T)
write.table(Tumor.MetaPs,"Meta.Ps_Tumor_6.23.14.txt",col.names=NA,sep="\t")
Normal.MetaPs=fisher.method(Normal.WilcP,na.rm=T)
write.table(Normal.MetaPs,"Meta.Ps_Normals_6.23.14.txt",col.names=NA,sep="\t")


# Fold Change Computation -------------------------------------------------
# Fold Change - mean, median, sum
# male/female
Tumor.ExpressionMale=Tumor.expression[,Tumor.dataframe$sex=="male"]
Tumor.ExpressionFemale=Tumor.expression[,Tumor.dataframe$sex=="female"]
Normal.ExpressionMale=Normal.expression[,Normal.dataframe$sex=="male"]
Normal.ExpressionFemale=Normal.expression[,Normal.dataframe$sex=="female"]

Tumor.FC=rowMeans(Tumor.ExpressionMale)/rowMeans(Tumor.ExpressionFemale)
Normal.FC=rowMeans(Normal.ExpressionMale)/rowMeans(Normal.ExpressionFemale)

Tumor.FCMed=apply(Tumor.ExpressionMale,1,median)/apply(Tumor.ExpressionFemale,1,median)
Normal.FCMed=apply(Normal.ExpressionMale,1,median)/apply(Normal.ExpressionFemale,1,median)

Tumor.onRatio=apply(Tumor.ExpressionMale,1,sum)/apply(Tumor.ExpressionFemale,1,sum)
Normal.onRatio=apply(Normal.ExpressionMale,1,sum)/apply(Normal.ExpressionFemale,1,sum)

# Correlation between normal and tumor fold changes
cor.test(Tumor.FC,Normal.FC,method="sp")
cor.test(Tumor.FCMed,Normal.FCMed,method="sp")
cor.test(Tumor.onRatio,Normal.onRatio,method="sp")
plot(Tumor.FC,Tumor.onRatio,main="Fold Change vs Number of nonzero (Tumor)")

# Aggregate and save all fold change analyses
FCs=cbind(Tumor.FC,Normal.FC,Tumor.FCMed,Normal.FCMed,Tumor.onRatio,Normal.onRatio)
rownames(FCs)=colnames(expression.data)
write.table(FCs,"FoldChanges_Gender.txt",col.names=NA,sep="\t")


# GWAS Like Analysis ------------------------------------------------------
# Note: Meta.Ps.combined.txt combines the p-values of male and female cases through meta-analysis
Full.pVals=read.table("Meta.Ps.combined.txt",sep="\t",header=T,row.names=1)
Full.pVals=Full.pVals[complete.cases(Full.pVals),]
# Preliminary Manhattan Plot
mhtplot(data.frame(Full.pVals))

Full.pVals=Full.pVals[order(Full.pVals[,"chr"],Full.pVals[,"pos"]),]
# pdf("ManhattanPlot_GenderPs_TumorsAndNormals.pdf") # uncomment to save to PDF
par(mfrow=c(2,1))
manhattanPlot(Full.pVals[,"Tumor_P.Adjust"],Full.pVals[,"UCSC_chrom.as"],ylim=c(0,15),signif=NULL,main="Manhattan Plot of Gender-Associated Transcripts in Tumor")
manhattanPlot(Full.pVals[,"Normal_P.Adjust"],Full.pVals[,"UCSC_chrom.as"],ylim=c(0,11),signif=NULL,main="Manhattan Plot of Gender-Associated Transcripts in Normal")
# dev.off()

# pdf("SignificanceAlongXandY.pdf") # uncomment to save to PDF
par(mfrow=c(1,2))

plot(Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="Y","UCSC_pos.as"],-log10(Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="Y","Tumor_P.value"]),col="blue",main="Significance Along Y-Chromosome",ylab="-log10(p)",xlab = "Y-Chromosomal Position")
plot(Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="X","UCSC_pos.as"],-log10(Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="X","Tumor_P.value"]),col="blue",main="Significance Along X-Chromosome",ylab="-log10(p)",xlab = "X-Chromosomal Position")
# dev.off()

plot(Full.pVals[,"Tumor_S"],Full.pVals[,"Normal_S"])
# Spearman's correlation coefficient: significance of correlation between tumor and non-tumor pvalues
cor.test(Full.pVals[,"Tumor_S"],Full.pVals[,"Normal_S"],method="sp")
# correlation between p-values of tumor and non-tumor in males
cor.test(Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="Y","Tumor_S"],Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="Y","Normal_S"],method="sp") 
# correlation between p-values of tumor and non-tumor in females
cor.test(Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="X","Tumor_S"],Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="X","Normal_S"],method="sp")
# tumors and non-tumors p-value significance in non-sex chromosomes
cor.test(Full.pVals[Full.pVals[,"UCSC_chrom.as"]!="X" & Full.pVals[,"UCSC_chrom.as"]!="Y","Tumor_S"],Full.pVals[Full.pVals[,"UCSC_chrom.as"]!="X" & Full.pVals[,"UCSC_chrom.as"]!="Y","Normal_S"])

# pdf("ScatterOfS_Stats_SexAndNonSex.pdf")
# par(mfrow=c(3,1))

#plots in same order as correlation tests above
plot(Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="Y","Tumor_S"],Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="Y","Normal_S"],method="sp")
plot(Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="X","Tumor_S"],Full.pVals[Full.pVals[,"UCSC_chrom.as"]=="X","Normal_S"])
plot(Full.pVals[Full.pVals[,"UCSC_chrom.as"]!="X" & Full.pVals[,"UCSC_chrom.as"]!="Y","Tumor_S"],Full.pVals[Full.pVals[,"UCSC_chrom.as"]!="X" & Full.pVals[,"UCSC_chrom.as"]!="Y","Normal_S"])
# dev.off()

# similar analysis in Normals
plot(-log10(Full.pVals[,"Tumor_P.value"]))
#p value vs chromosome number
boxplot(-log10(Full.pVals[,"Tumor_P.value"])~Full.pVals[,"UCSC_chrom.as"])

#not really used for anything - still included.
miss=read.table("noChrom.txt",sep="\t")
miss=as.numeric(miss[,1])
genes=unlist(lapply(strsplit(row.names(Tumor.MetaPs),"|",fixed=T),function(x)(x[[2]])))


# Gene Set Analysis (GSA) -------------------------------------------------

# Gene Set Analysis -- this section was not used because analysis was too general
library(GSA)
# entrez is one of the forms used to identify genes
geneset.obj<- GSA.read.gmt("msigdb.v4.0.entrez.gmt")
# 3D array: genesets x number of attributes x number of tissue types
outArray=array(NA,dim=c(length(geneset.obj$genesets),5,length(unique(Tumor.dataframe$TissueType))),dimnames=list(genesets=geneset.obj$geneset.names,attributes=c("Gene_set","Gene_set_name","Score","p-value","FDR"),unique(Tumor.dataframe$TissueType)))
# x = all tumor expression data for one tissue type
# y = male/female matrix with male = 1 and female = 2 (a vector)
# gene names are row names of expression data
# multiple cases of one tissue type so only loop over unique tissue types
for(i in 1:length(unique(Tumor.dataframe$TissueType))){
  x=Tumor.expression[,Tumor.dataframe$TissueType==unique(Tumor.dataframe$TissueType)[i]]
  y=Tumor.dataframe$sex[Tumor.dataframe$TissueType==unique(Tumor.dataframe$TissueType)[i]]
  y[y=="male"]=1
  y[y=="female"]=2
  # The second string of the split string is the numeric value of the gene
  # genenames is a vector of all numeric values of genes
  genenames=unlist(lapply(strsplit(rownames(Tumor.expression),"|",fixed=T),function(xx)(xx[2])))
  genesets=geneset.obj$genesets
  geneset.names=geneset.obj$geneset.names
  # Note: 'genenames' and 'genesets' are the names of the parameters of the GSA function
  # returns significance, measured in p-values, of sets of genes in genesets as applied to genenames
  GSA.obj<-GSA(x,y, genenames=genenames, genesets=genesets, resp.type="Two class unpaired", nperms=100)
  # list out the gene sets with false discovery rates (p values) lower than the false discovery rate cutpoint (FDRcut).
  # FDR=false discovery rate (i.e. rate of false positives), which is the same as the significance level
  out=GSA.listsets(GSA.obj, geneset.names=geneset.names,FDRcut=1,maxchar=9999)
  # combine the positive and negative into a table by rows
  outTable=rbind(out$positive,out$negative)
  write.table(outTable,file=paste('GSA_',unique(Tumor.dataframe$TissueType)[i]),sep='\t')
}

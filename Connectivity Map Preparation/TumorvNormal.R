#Separates data into males and females, separates them by tissue type, and generates probes for cmap
#This runs the wilcoxon rank-sum test with respect to tissue type (tumor vs normal)
#This script essentially constructs the resistance and sensitivity signatures for each cancer
#Jonathan Ma and Sadhika Malladi
#7.4.2014

setwd("C:/Users/Jonathan Ma/Google Drive/ExpressionSexProject/Expression_Signatures")
load("TCGA.Clin.RNAseq.v2.RData")
library(gdata)
library(MADAM)
#to install GWASTools package: 
# source("http://bioconductor.org/biocLite.R")
# biocLite("GWASTools")
library(GWASTools)
library(gap)
library(NCBI2R)
library(GSA)
library(jetset)
expression.data <- dat
clinical.data <- clin

expression.data=t(expression.data)

expression.dataframe <- data.frame(t(data.frame(strsplit(rownames(expression.data),"-|\\."),stringsAsFactors=F)))
expression.dataframe$Barcode <- apply(expression.dataframe,MARGIN=1,function(x)paste(x[c(2:4)],collapse="-"))
expression.dataframe <- expression.dataframe[,c("X1","X5","Barcode")]
colnames(expression.dataframe) <- c("TissueType","SampleType","Barcode")
expression.dataframe$SampleType=substr(expression.dataframe$SampleType,1,2)

rownames(clinical.data) <- toupper(gsub("[A-Z]+.(tcga).([a-z0-9]+).([a-z0-9]+)","\\1-\\2-\\3",rownames(clinical.data)))
clinical.data=clinical.data[rownames(clinical.data)%in%expression.dataframe$Barcode,]
expression.dataframe=expression.dataframe[expression.dataframe$Barcode%in%rownames(clinical.data),]

expression.dataframe$sex <- clinical.data[ expression.dataframe$Barcode,"gender"]
male.dataframe=expression.dataframe[expression.dataframe$sex=="male",]
male.dataframe=male.dataframe[male.dataframe$Barcode %in% rownames(clinical.data),]
female.dataframe=expression.dataframe[expression.dataframe$sex=="female",]
female.dataframe=female.dataframe[female.dataframe$Barcode %in% rownames(clinical.data),]

rownames(expression.data)=gsub("[A-Z]+-TCGA.([A-Z0-9]+).([A-Z0-9]+).*","TCGA-\\1-\\2",rownames(expression.data))
male.expression=expression.data[male.dataframe$Barcode,]
female.expression=expression.data[female.dataframe$Barcode,]

male.dataframe$Category=NA
temp=male.dataframe$Category
temp[male.dataframe$SampleType=="01" | male.dataframe$SampleType=="03"]="Tumor"
temp[male.dataframe$SampleType=="11"]="Normal"
male.dataframe$Category=temp

female.dataframe$Category=NA
temp=female.dataframe$Category
temp[female.dataframe$SampleType=="01" | female.dataframe$SampleType=="03"]="Tumor"
temp[female.dataframe$SampleType=="11"]="Normal"
female.dataframe$Category=temp

male.expression=t(male.expression)
female.expression=t(female.expression)

#only do calculations on cancers for which we have clinical data, which are those in male.dataframe and female.dataframe
#also sorts data so that they match the order in male.dataframe and female.dataframe respectively
male.expression=male.expression[,male.dataframe$Barcode]
female.expression=female.expression[,female.dataframe$Barcode]

#remove all barcodes for which category is NA
male.dataframe=male.dataframe[!is.na(male.dataframe$Category),]
female.dataframe=female.dataframe[!is.na(female.dataframe$Category),]
male.expression=male.expression[,male.dataframe$Barcode]
female.expression=female.expression[,female.dataframe$Barcode]

#only choose cancers for which there are at least 10 male and 10 female cases in normals
normal.tissuetypes=unique(expression.dataframe$TissueType)
normal.tissuetypes=normal.tissuetypes[!is.na(normal.tissuetypes)]
indices=integer(0)
for(i in 1:length(normal.tissuetypes))
{
  if(length(male.dataframe$TissueType[male.dataframe$TissueType==normal.tissuetypes[i] & male.dataframe$Category=="Normal"])>10 & length(female.dataframe$TissueType[female.dataframe$TissueType==normal.tissuetypes[i] & female.dataframe$Category=="Normal"])>10)
  {
    indices=c(indices,i)
  }
}
normal.tissuetypes=normal.tissuetypes[indices]

#only choose cancers for which there are at least 10 male and 10 female cases in tumors
tumor.tissuetypes=unique(expression.dataframe$TissueType)
tumor.tissuetypes=tumor.tissuetypes[!is.na(tumor.tissuetypes)]
indices=integer(0)
for(i in 1:length(tumor.tissuetypes))
{
  if(length(male.dataframe$TissueType[male.dataframe$TissueType==tumor.tissuetypes[i] & male.dataframe$Category=="Tumor"])>10 & length(female.dataframe$TissueType[female.dataframe$TissueType==tumor.tissuetypes[i] & female.dataframe$Category=="Tumor"])>10)
  {
    indices=c(indices,i)
  }
}
tumor.tissuetypes=tumor.tissuetypes[indices]

target.tissuetypes=intersect(tumor.tissuetypes,normal.tissuetypes)

male.dataframe=male.dataframe[male.dataframe$TissueType%in%target.tissuetypes,]
female.dataframe=female.dataframe[female.dataframe$TissueType%in%target.tissuetypes,]
male.expression=male.expression[,male.dataframe$Barcode]
female.expression=female.expression[,female.dataframe$Barcode]

#remove all genes which are 0 in all cases for at least 1 tissue
current.tissue.rows=integer(0)
overall.rows=integer(0)
for(i in 1:length(target.tissuetypes))
{
  male.tissueCases<-male.dataframe$Barcode[male.dataframe$TissueType==target.tissuetypes[i]]
  male.tissueExpression<-male.expression[,male.tissueCases]
  female.tissueCases<-female.dataframe$Barcode[female.dataframe$TissueType==target.tissuetypes[i]]
  female.tissueExpression<-female.expression[,female.tissueCases]
  male.category=male.category<-male.dataframe$Category[male.dataframe$TissueType==target.tissuetypes[i]]
  female.category=female.category<-female.dataframe$Category[female.dataframe$TissueType==target.tissuetypes[i]]
  current.tissue.rows=integer(0)
  for(j in 1:nrow(male.expression))
  {
    if(length(unique(male.tissueExpression[j,]))>1 & length(unique(female.tissueExpression[j,]))>1)
    {
      current.tissue.rows=c(current.tissue.rows,j)
    }
  }
  if(i==1)
  {
    overall.rows=current.tissue.rows
  }
  else
  {
    overall.rows=intersect(overall.rows,current.tissue.rows)
  }
}
male.expression=male.expression[overall.rows,]
female.expression=female.expression[overall.rows,]

#Wilcoxon tests on tumor vs normal in males and females
male.wilcP=matrix(NA,nrow=nrow(male.expression),ncol=length(target.tissuetypes),dimnames=list(rownames(male.expression),target.tissuetypes))
female.wilcP=matrix(NA,nrow=nrow(female.expression),ncol=length(target.tissuetypes),dimnames=list(rownames(female.expression),target.tissuetypes))
male.wilcEst=male.wilcP
female.wilcEst=female.wilcP

for(i in 1:ncol(male.wilcP)) {
  #get all barcodes of the current tissue type
  male.tissueCases<-male.dataframe$Barcode[male.dataframe$TissueType==target.tissuetypes[i]]
  #expression values of all genes for current tissue type (genes=rows, barcodes=columns)
  male.tissueExpression<-male.expression[,male.tissueCases]
  #a list of genders of all barcodes (cases) of the current tissue type
  male.category<-male.dataframe$Category[male.dataframe$TissueType==target.tissuetypes[i]]
  for(j in 1:nrow(male.tissueExpression))
  {
    if(j%%1000==0)
    {
      cat("\nGene ",j," in tissue ",colnames(male.wilcP)[i])
    }
    temp.dataframe=data.frame(male.tissueExpression[j,],male.category)
    twilc=wilcox.test(temp.dataframe[,1]~temp.dataframe[,2],conf.int=T)
    male.wilcP[j,i]=twilc$p.value
    male.wilcEst[j,i]=twilc$est
  }
  cat("\nTissue ", colnames(male.wilcP)[i], "in males completed")
}
saveRDS(male.wilcP,"male.wilcP.tvn.RDS")
saveRDS(male.wilcEst,"male.wilcEst.tvn.RDS")

for(i in 1:ncol(female.wilcP)) {
  #get all barcodes of the current tissue type
  female.tissueCases<-female.dataframe$Barcode[female.dataframe$TissueType==target.tissuetypes[i]]
  #expression values of all genes for current tissue type (genes=rows, barcodes=columns)
  female.tissueExpression<-female.expression[,female.tissueCases]
  #a list of genders of all barcodes (cases) of the current tissue type
  female.category<-female.dataframe$Category[female.dataframe$TissueType==target.tissuetypes[i]]
  for(j in 1:nrow(female.tissueExpression))
  {
    if(j%%1000==0)
    {
      cat("\nGene ",j," in tissue ",colnames(male.wilcP)[i])
    }
    temp.dataframe=data.frame(female.tissueExpression[j,],female.category)
    twilc=wilcox.test(temp.dataframe[,1]~temp.dataframe[,2],conf.int=T)
    female.wilcP[j,i]=twilc$p.value
    female.wilcEst[j,i]=twilc$est
  }
  cat("\nTissue ", colnames(female.wilcP)[i], " in females completed")
}

saveRDS(female.wilcP,"female.wilcP.tvn.RDS")
saveRDS(female.wilcEst,"female.wilcEst.tvn.RDS")

setwd("C:/Users/Jonathan Ma/Google Drive/ExpressionSexProject/Connectivity_Map_v2")
#only run this section for testing
male.wilcP=readRDS('male.wilcP.tvn.RDS')
female.wilcP=readRDS('female.wilcP.tvn.RDS')


#manually adjust p values in order to find the most significant genes in each cancer type separately
adj.male.wilcP=c(male.wilcP)
adj.male.wilcP=p.adjust(male.wilcP,method='fdr')
adj.male.wilcP=matrix(male.wilcP,ncol=ncol(male.wilcP),byrow=F)
rownames(adj.male.wilcP)=rownames(male.wilcP)
colnames(adj.male.wilcP)=colnames(male.wilcP)
saveRDS(adj.male.wilcP,'adj.male.wilcP.tvn.RDS')

adj.female.wilcP=c(female.wilcP)
adj.female.wilcP=p.adjust(female.wilcP,method='fdr')
adj.female.wilcP=matrix(female.wilcP,ncol=ncol(female.wilcP),byrow=F)
rownames(adj.female.wilcP)=rownames(female.wilcP)
colnames(adj.female.wilcP)=colnames(female.wilcP)
saveRDS(adj.female.wilcP,'adj.female.wilcP.tvn.RDS')

# Fold Change Analysis: tumor vs normal
male.fc.matrix=matrix(NA,nrow=nrow(male.expression),ncol=length(target.tissuetypes),dimnames=list(rownames(male.expression),target.tissuetypes))
female.fc.matrix=matrix(NA,nrow=nrow(female.expression),ncol=length(target.tissuetypes),dimnames=list(rownames(female.expression),target.tissuetypes))
for(i in 1:length(target.tissuetypes))
{
  current.male.tumor.dataframe=male.dataframe[male.dataframe$TissueType==target.tissuetypes[i] & male.dataframe$Category=="Tumor",]
  current.male.normal.dataframe=male.dataframe[male.dataframe$TissueType==target.tissuetypes[i] & male.dataframe$Category=="Normal",]
  current.male.tumor.expression=male.expression[,current.male.tumor.dataframe$Barcode]
  current.male.normal.expression=male.expression[,current.male.normal.dataframe$Barcode]
  current.female.tumor.dataframe=female.dataframe[female.dataframe$TissueType==target.tissuetypes[i] & female.dataframe$Category=="Tumor",]
  current.female.normal.dataframe=female.dataframe[female.dataframe$TissueType==target.tissuetypes[i] & female.dataframe$Category=="Normal",]
  current.female.tumor.expression=female.expression[,current.female.tumor.dataframe$Barcode]
  current.female.normal.expression=female.expression[,current.female.normal.dataframe$Barcode]
  
  male.fc.matrix[,i]=apply(current.male.tumor.expression,1,median)/apply(current.male.normal.expression,1,median)
  female.fc.matrix[,i]=apply(current.female.tumor.expression,1,median)/apply(current.female.normal.expression,1,median)
}
saveRDS(male.fc.matrix,"male.fc.matrix.tvn.RDS")
saveRDS(female.fc.matrix,"female.fc.matrix.tvn.RDS")


# Meta Analysis: Fisher's Method for tumor vs normal
setwd("C:/Users/Jonathan Ma/Google Drive/ExpressionSexProject/Results/Pathway Analysis/Meta")
male.metaP=fisher.method(male.wilcP,na.rm=T)
female.metaP=fisher.method(female.wilcP,na.rm=T)
saveRDS(male.metaP,'male.metaP.tvn.RDS')
saveRDS(female.metaP,'female.metaP.tvn.RDS')

overall.male.fc=apply(male.expression[,male.dataframe$Category=="Tumor"],1,median)/apply(male.expression[,male.dataframe$Category=="Normal"],1,median)
overall.female.fc=apply(female.expression[,female.dataframe$Category=="Tumor"],1,median)/apply(female.expression[,female.dataframe$Category=="Normal"],1,median)
saveRDS(overall.male.fc,"overall.male.fc.RDS")
saveRDS(overall.female.fc,"overall.female.fc.RDS")

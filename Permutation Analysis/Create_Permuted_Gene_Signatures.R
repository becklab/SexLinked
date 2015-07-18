#Code intended to validate LincsCloud results through Permutation Analysis, inspired by and modeled after the PMASE approach of Liu et al. (2011)
#Creates permuted gene expression signatures for input into LincsCloud C3 (Compute Connectivity on the Cloud)
#This approach uses the Wilcoxon option of the Significance Analysis of Microarrays method, which is much faster than the default wilcox.test function in R.
#Author: Jonathan Ma
#Created: 2/20/2015
#
#PLEASE NOTE: This code may cause the console to stop outputting. To remedy this, just run the command sink(NULL) in the console.

library(gdata)
library(MADAM)
library(GWASTools)
library(gap)
library(NCBI2R)
library(GSA)
library(jetset)
library(GSEABase)
library(Biobase)
library(BiocGenerics)
library(hgu133a2.db)
library(biomaRt)
library(siggenes)
library(multtest)
library(samr)

#set to your own working directory
setwd("C:/Users/Jonathan Ma/Google Drive/ExpressionSexProject/Expression_Signatures")

#read in local data files created from cleaned TCGA expression data.
te=readRDS(file="TumorExpression.RDS")
ne=readRDS(file="NormalExpression.RDS")
td=readRDS(file="TumorDataframe.RDS")
nd=readRDS(file="NormalDataframe.RDS")

#specify list of cancers for which we want to predict drugs
#Meta denotes a pan-cancer meta-analysis
cancers=c("HNSC","KICH","KIRC","LIHC","LUAD","LUSC","THCA","Meta")
td=td[td$TissueType%in%cancers,]
nd=nd[nd$TissueType%in%cancers,]

#removes cases that have "0" (i.e. less than 10) males or females 
#needed because as.factor keeps levels
td$TissueType=drop.levels(td$TissueType)

#specify your number of permutations here. 1000 chosen as convenience to limit LincsCloud C3 computational time
numberofpermutations=1000

allgenes=unique(c(rownames(te),rownames(ne)))

#create list of affymetrix id to hugo gene conversions
hugo.allgenes=gsub("(.*)\\|.*","\\1",allgenes)
hugo.allgenes=hugo.allgenes[!is.na(hugo.allgenes)]
hugo.allgenes=hugo.allgenes[hugo.allgenes!="?"]
hugo.allgenes=unique(hugo.allgenes)

mart=useMart('ensembl',dataset='hsapiens_gene_ensembl')
all.affy=getBM(attributes=c('affy_hg_u133_plus_2','hgnc_symbol'),filters='hgnc_symbol',values=hugo.allgenes,mart=mart)



numberofpermutations=1000

for(x in 1:numberofpermutations)
{
  #randomly reorder td$sex and nd$sex to permute gender at the sample level
  perm.td=td
  perm.nd=nd
  indices=sample(1:nrow(perm.td),nrow(perm.td),replace=F)
  perm.td$sex=perm.td$sex[indices]
  indices=sample(1:nrow(perm.nd),nrow(perm.nd),replace=F)
  perm.nd$sex=perm.nd$sex[indices]
  
  pmtd=perm.td[perm.td$sex=='male',]
  pftd=perm.td[perm.td$sex=='female',]
  pmnd=perm.nd[perm.nd$sex=='male',]
  pfnd=perm.nd[perm.nd$sex=='female',]
  
  pmte=te[,colnames(te)%in%pmtd$Barcode]
  pfte=te[,colnames(te)%in%pftd$Barcode]
  pmne=ne[,colnames(ne)%in%pmnd$Barcode]
  pfne=ne[,colnames(ne)%in%pfnd$Barcode]
  
  pme=cbind(pmte,pmne)
  pfe=cbind(pfte,pfne)
  
  pmtd=cbind(pmtd,"tumor")
  colnames(pmtd)[6]="TissueCategory"
  pftd=cbind(pftd,"tumor")
  colnames(pftd)[6]="TissueCategory"
  pmnd=cbind(pmnd,"normal")
  colnames(pmnd)[6]="TissueCategory"
  pfnd=cbind(pfnd,"normal")
  colnames(pfnd)[6]="TissueCategory"
  
  pmd=rbind(pmtd,pmnd)
  pfd=rbind(pftd,pfnd)
  
  pme=pme[,pmd$Barcode]
  pfe=pfe[,pfd$Barcode]
  
  #p value matrices: no data is inputted here, only dimensions are set
  #p values indicate significance of gene expression disparities between tumor and normal
  # row = genes, col = tissue type (unique)
  #unique discards duplicated values
  male.tvn.wilcP <- matrix(NA,nrow=nrow(pme),ncol=length(unique(pmd$TissueType)),dimnames=list(rownames(pme),unique(pmd$TissueType)))
  female.tvn.wilcP <- matrix(NA,nrow=nrow(pfe),ncol=length(unique(pfd$TissueType)),dimnames=list(rownames(pfe),unique(pfd$TissueType)))
  
  #calculate male tumor vs normal wilcoxon p values using the Wilcoxon option of SAM, based on these sex permutations
  for(i in colnames(male.tvn.wilcP)) { #iterate thru all tissue types
    #get all barcodes of the current tissue type
    male.tissueCases=pmd$Barcode[pmd$TissueType==i]
    female.tissueCases=pfd$Barcode[pfd$TissueType==i]
    #get expression values of all genes for the current tissue type
    male.tissueExpression=pme[,male.tissueCases]
    female.tissueExpression=pfe[,female.tissueCases]
    #get a list of all tissue categories (i.e. tumor or normal) for the current tissue type
    male.tissueCategories=pmd$TissueCategory[pmd$TissueType==i]
    female.tissueCategories=pfd$TissueCategory[pfd$TissueType==i]
    
    #0=normal, 1=tumor
    #suppressWarnings suppresses warnings
    #capture.output stops the sam() function from printing "currently under testing." the SAM object that is returned by sam() is stored in male.sam or female.sam
    #the message "currently under testing" is stored under the variable "uselessmessage"
    uselessmessage=capture.output({male.sam=suppressWarnings(sam(male.tissueExpression,male.tissueCategories,method=wilc.stat,use.row=T))})
    uselessmessage2=capture.output({female.sam=suppressWarnings(sam(female.tissueExpression,female.tissueCategories,method=wilc.stat,use.row=T))})
    male.tvn.wilcP[,i]=male.sam@p.value
    female.tvn.wilcP[,i]=female.sam@p.value
  }
  
  male.wilcP=male.tvn.wilcP
  female.wilcP=female.tvn.wilcP
  
  #eliminate NA's
  for(i in 1:ncol(male.wilcP))
  {
    male.wilcP=male.wilcP[!is.na(male.wilcP[,i]),]
    male.wilcP=male.wilcP[male.wilcP[,i]!=0,]
    
    female.wilcP=female.wilcP[!is.na(female.wilcP[,i]),]
    female.wilcP=female.wilcP[female.wilcP[,i]!=0,]
  }
  
  #aggregate cancer-specific p-values into pan-cancer p-values using Fisher's method
  #done so that we can include a meta-analysis
  male.metaP=fisher.method(male.wilcP,na.rm=T)
  female.metaP=fisher.method(female.wilcP,na.rm=T)
  
  
  #calculate benjamini-hochberg values by adjusting for multiple significance testing
  adj.male.wilcP.vector=c(male.wilcP)
  adj.male.wilcP.vector=p.adjust(adj.male.wilcP.vector,method='fdr')
  adj.male.wilcP=matrix(adj.male.wilcP.vector,ncol=ncol(male.wilcP),byrow=F)
  rownames(adj.male.wilcP)=rownames(male.wilcP)
  colnames(adj.male.wilcP)=colnames(male.wilcP)
  
  adj.female.wilcP.vector=c(female.wilcP)
  adj.female.wilcP.vector=p.adjust(adj.female.wilcP.vector,method='fdr')
  adj.female.wilcP=matrix(adj.female.wilcP.vector,ncol=ncol(female.wilcP),byrow=F)
  rownames(adj.female.wilcP)=rownames(female.wilcP)
  colnames(adj.female.wilcP)=colnames(female.wilcP)
  
  adj.male.wilcP=cbind(adj.male.wilcP,male.metaP$p.adj)
  adj.female.wilcP=cbind(adj.female.wilcP,female.metaP$p.adj)
  colnames(adj.male.wilcP)[ncol(adj.male.wilcP)]="Meta"
  colnames(adj.female.wilcP)[ncol(adj.female.wilcP)]="Meta"
  
  #prepare matrices to calculate gene expression fold changes of tumor vs normal
  male.fc.matrix=matrix(NA,nrow=nrow(pme),ncol=ncol(male.tvn.wilcP),dimnames=list(rownames(pme),colnames(male.tvn.wilcP)))
  female.fc.matrix=matrix(NA,nrow=nrow(pfe),ncol=ncol(female.tvn.wilcP),dimnames=list(rownames(pfe),colnames(female.tvn.wilcP)))
  
  #calculate tumor vs normal vold changes
  for(i in colnames(male.tvn.wilcP)) { #iterate thru all tissue types including Meta
    male.tumor.tissueCases=pmtd$Barcode[pmtd$TissueType==i]
    female.tumor.tissueCases=pftd$Barcode[pftd$TissueType==i]
    male.normal.tissueCases=pmnd$Barcode[pmnd$TissueType==i]
    female.normal.tissueCases=pfnd$Barcode[pfnd$TissueType==i]
    
    male.tumor.tissueExpression=pmte[,male.tumor.tissueCases]
    male.normal.tissueExpression=pmne[,male.normal.tissueCases]
    female.tumor.tissueExpression=pfte[,female.tumor.tissueCases]
    female.normal.tissueExpression=pfne[,female.normal.tissueCases]
    
    male.fc.matrix[,i]=apply(male.tumor.tissueExpression,1,median)/apply(male.normal.tissueExpression,1,median)
    female.fc.matrix[,i]=apply(female.tumor.tissueExpression,1,median)/apply(female.normal.tissueExpression,1,median)
  }
  
  #calculate pan-cancer tumor vs normal expression fold changes
  male.meta.fc.vector=apply(pmte[rownames(male.fc.matrix),],1,median)/apply(pmne[rownames(male.fc.matrix),],1,median)
  female.meta.fc.vector=apply(pfte[rownames(female.fc.matrix),],1,median)/apply(pfne[rownames(female.fc.matrix),],1,median)
  
  male.fc.matrix=cbind(male.fc.matrix,male.meta.fc.vector)
  female.fc.matrix=cbind(female.fc.matrix,female.meta.fc.vector)
  colnames(male.fc.matrix)[ncol(male.fc.matrix)]="Meta"
  colnames(female.fc.matrix)[ncol(female.fc.matrix)]="Meta"
  
  #for each cancer (as well as meta-analysis) find genes significantly differentially expressed between tumors and normals within each sex,
  #then sort them by fold change of expression. From this, identify sensitivity signatures as the 250 genes with lowest p values and fold change
  #<1, and resistance signatures as the 250 genes with lowest p values and fold change >1
  for(i in cancers)
  {
    male.up=rownames(male.fc.matrix)[rownames(male.fc.matrix)%in%rownames(adj.male.wilcP)]
    relevant.fc.matrix=male.fc.matrix[male.up,]
    male.up=male.up[relevant.fc.matrix[,i]>1 & adj.male.wilcP[,i]<0.05] #overexpressed in male tumors relative to male normals
    sorted.male.pvals=sort(adj.male.wilcP[,i],decreasing=F)
    sorted.male.pvals=cbind(sorted.male.pvals,integer(0))
    sorted.male.up=rownames(sorted.male.pvals)[rownames(sorted.male.pvals)%in%male.up]
    
    female.up=rownames(female.fc.matrix)[rownames(female.fc.matrix)%in%rownames(adj.female.wilcP)]
    relevant.fc.matrix=female.fc.matrix[female.up,]
    female.up=female.up[relevant.fc.matrix[,i]>1 & adj.female.wilcP[,i]<0.05] #overexpressed in female tumors relative to female normals
    sorted.female.pvals=sort(adj.female.wilcP[,i],decreasing=F)
    sorted.female.pvals=cbind(sorted.female.pvals,integer(0))
    sorted.female.up=rownames(sorted.female.pvals)[rownames(sorted.female.pvals)%in%female.up]
    
    male.down=rownames(male.fc.matrix)[rownames(male.fc.matrix)%in%rownames(adj.male.wilcP)]
    relevant.fc.matrix=male.fc.matrix[male.down,]
    male.down=male.down[relevant.fc.matrix[,i]<1 & adj.male.wilcP[,i]<0.05] #overexpressed in male normals relative to male tumors
    sorted.male.pvals=sort(adj.male.wilcP[,i],decreasing=F)
    sorted.male.pvals=cbind(sorted.male.pvals,integer(0))
    sorted.male.down=rownames(sorted.male.pvals)[rownames(sorted.male.pvals)%in%male.down]
    
    female.down=rownames(female.fc.matrix)[rownames(female.fc.matrix)%in%rownames(adj.female.wilcP)]
    relevant.fc.matrix=female.fc.matrix[female.down,]
    female.down=female.down[relevant.fc.matrix[,i]<1 & adj.female.wilcP[,i]<0.05] #overexpressed in female normals relative to female tumors
    sorted.female.pvals=sort(adj.female.wilcP[,i],decreasing=F)
    sorted.female.pvals=cbind(sorted.female.pvals,integer(0))
    sorted.female.down=rownames(sorted.female.pvals)[rownames(sorted.female.pvals)%in%female.down]
    
    #only use HUGO gene symbols
    male.up.hugo=gsub("(.*)\\|.*","\\1",sorted.male.up)
    female.up.hugo=gsub("(.*)\\|.*","\\1",sorted.female.up)
    male.down.hugo=gsub("(.*)\\|.*","\\1",sorted.male.down)
    female.down.hugo=gsub("(.*)\\|.*","\\1",sorted.female.down)
    
    #eliminate genes which have no hugo symbols
    male.up.hugo=male.up.hugo[!is.na(male.up.hugo)]
    male.down.hugo=male.down.hugo[!is.na(male.down.hugo)]
    female.up.hugo=female.up.hugo[!is.na(female.up.hugo)]
    female.down.hugo=female.down.hugo[!is.na(female.down.hugo)]
    
    male.up.hugo=male.up.hugo[male.up.hugo!="?"]
    female.up.hugo=female.up.hugo[female.up.hugo!="?"]
    male.down.hugo=male.down.hugo[male.down.hugo!="?"]
    female.down.hugo=female.down.hugo[female.down.hugo!="?"]
    
    unique.male.up.hugo=unique(male.up.hugo)
    unique.male.down.hugo=unique(male.down.hugo)
    unique.female.up.hugo=unique(female.up.hugo)
    unique.female.down.hugo=unique(female.down.hugo)
    
    if(length(unique.male.up.hugo)>250)
    {
      male.up.hugo.truncated=unique.male.up.hugo[1:250]
    } else
    {
      male.up.hugo.truncated=unique.male.up.hugo
    }
    
    if(length(unique.female.up.hugo)>250)
    {
      female.up.hugo.truncated=unique.female.up.hugo[1:250]
    } else
    {
      female.up.hugo.truncated=unique.female.up.hugo
    }
    
    if(length(unique.male.down.hugo)>250)
    {
      male.down.hugo.truncated=unique.male.down.hugo[1:250]
    } else
    {
      male.down.hugo.truncated=unique.male.down.hugo
    }
    
    if(length(unique.female.down.hugo)>250)
    {
      female.down.hugo.truncated=unique.female.down.hugo[1:250]
    } else
    {
      female.down.hugo.truncated=unique.female.down.hugo
    }
    
    #convert to affymetrix using the all.affy data frame from earlier
    #map HUGO gene symbols to AffyMetrix Gene ID's
    male.up.affy=all.affy[all.affy[,2]%in%male.up.hugo.truncated,1]
    female.up.affy=all.affy[all.affy[,2]%in%female.up.hugo.truncated,1]
    male.down.affy=all.affy[all.affy[,2]%in%male.down.hugo.truncated,1]
    female.down.affy=all.affy[all.affy[,2]%in%female.down.hugo.truncated,1]
    
    male.up.affy=male.up.affy[male.up.affy!=""]
    female.up.affy=female.up.affy[female.up.affy!=""]
    male.down.affy=male.down.affy[male.down.affy!=""]
    female.down.affy=female.down.affy[female.down.affy!=""]
    
    unique.male.up.affy=unique(male.up.affy)
    unique.female.up.affy=unique(female.up.affy)
    unique.male.down.affy=unique(male.down.affy)
    unique.female.down.affy=unique(female.down.affy)
    
    #convert signatures to GeneSetCollections
    male.up.gs=GeneSet(unique.male.up.affy)
    female.up.gs=GeneSet(unique.female.up.affy)
    male.down.gs=GeneSet(unique.male.down.affy)
    female.down.gs=GeneSet(unique.female.down.affy)
    male.up.collection=GeneSetCollection(male.up.gs)
    female.up.collection=GeneSetCollection(female.up.gs)
    male.down.collection=GeneSetCollection(male.down.gs)
    female.down.collection=GeneSetCollection(female.down.gs)
    
    #save the signatures (now converted to GeneSetCollections) in GMT format, required by LincsCloud C3
    setwd(paste0("C:/Users/Jonathan Ma/Documents/Research/LincsCloud_C3/Permuted_Signatures_SAM/",i))
    toGmt(male.up.collection,con=paste0(x,"_",i,"_male_permuted_upInTumors.gmt"))
    toGmt(female.up.collection,con=paste0(x,"_",i,"_female_permuted_upInTumors.gmt"))
    toGmt(male.down.collection,con=paste0(x,"_",i,"_male_permuted_downInTumors.gmt"))
    toGmt(female.down.collection,con=paste0(x,"_",i,"_female_permuted_downInTumors.gmt"))
  }
  cat("\nFinished generating iteration ",x)
  
}
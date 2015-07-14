#Analyzes copy number data downloaded from Cbioportal
#Jonathan Ma and Sadhika Malladi
#6.19.2014

# This script computes the wilcoxon rank-sum test with respect to sex on 
# DNA copy number data. The resulting p-values are then adjusted with the Benjamini
# Hochberg method and saved.
# Note: the script is complex because of the difficulty of dealing with cBioPortal data.

# Set your working directory to the appropriate path!
setwd("C:/Users/admin/Google Drive/ExpressionSexProject/DNA_Copy_Number")

#data not available on cbioportal for READ, COAD, and LGG
fileNames=dir(pattern="copy_v2.RDS")
allCancerTypes=gsub("([a-z]+)_copy_v2.RDS","\\1",fileNames)
allCancerTypes=toupper(allCancerTypes)

#since clinical data from cbioportal itself is sparse, need to use master TCGA clinical data file
clinical.data=readRDS('MasterClinical.RDS')
#convert clinical data barcodes to standard format by making them uppercase
clinical.data$barcode <- toupper(clinical.data$barcode)

#wilc.mat: rows=probes, 2 columns (1 for est, 1 for p values)
#each column of wilc.mat is a tissue type
#each row of wilc.mat is a gene in the dna copy number data file (and each data file has the same genes)
wilcP <- matrix(c(1,0),nrow=nrow(readRDS(fileNames[1])),ncol=length(fileNames),byrow=T)
wilcEst <- matrix(c(1,0),nrow=nrow(readRDS(fileNames[1])),ncol=length(fileNames),byrow=T)
all.copy.data=readRDS(fileNames[1])
colnames(all.copy.data)=gsub("TCGA.([A-Z0-9]+).([A-Z0-9]+)","TCGA-\\1-\\2",colnames(all.copy.data))
removeRowsList=c(1)
removeRowsList=removeRowsList[-1]
naValues=logical(0)
naIndexes=integer(0)

naRows=logical(0)
currentRowNAs=logical(0)
currentRowNAValues=integer(0)
naRowIndexes=integer(0)
overallNARowIndexes=integer(0)

for(i in 1:length(fileNames))
{ 
  copy.data=readRDS(fileNames[i])
  colnames(copy.data)=gsub("TCGA.([A-Z0-9]+).([A-Z0-9]+)","TCGA-\\1-\\2",colnames(copy.data))
  #filter to only include barcodes for which we have clinical data
  copy.data=copy.data[,colnames(copy.data)%in%clinical.data$barcode]
  tissue.clinical=clinical.data[clinical.data$barcode%in%colnames(copy.data),]
  
  for(currentRow in 1:nrow(copy.data))
  {
    currentRowNAs=is.na(copy.data[currentRow,])
    currentRowNAValues=c(1:ncol(copy.data))
    currentRowNAValues=currentRowNAValues[currentRowNAs]
    naRows=c(naRows,length(currentRowNAValues)==ncol(copy.data))
  }
  
  naRowIndexes=c(1:nrow(copy.data)) #rows that are all NA's for the current copy.data
  naRowIndexes=naRowIndexes[naRows]
  overallNARowIndexes=c(overallNARowIndexes,naRowIndexes)
  overallNARowIndexes=unique(overallNARowIndexes)
  overallNARowIndexes=sort(naRowIndexes,decreasing=T)
  
  
  cat("\nComputing wilcoxon test for cancer type: ",allCancerTypes[i])
  for(j in 1:nrow(copy.data)) #j is the current gene. we are looping over all genes in copy.data
  {
    numUnique=0
    numUnique=length(unique(copy.data[j,]))
    if(numUnique>2 && !(j %in% naRowIndexes))
    {
      cat("\nAnalyzing gene ",j," of ",nrow(copy.data))
      temp.dataframe=data.frame(copy.data[j,],tissue.clinical$gender)
      colnames(temp.dataframe)=c("copynumbers","gender")
      twilc=wilcox.test(temp.dataframe$copynumbers~temp.dataframe$gender,conf.int=T)
      wilcEst[j,i]=twilc$est
      wilcP[j,i]=twilc$p.value
    }
    else
    {
      removeRowsList=c(removeRowsList,j)
    }
  }
  
  saveRDS(wilcEst[,i],paste0(allCancerTypes[i],"_copy_wilcEst_Unremoved.RDS",sep=""))
  saveRDS(wilcP[,i],paste0(allCancerTypes[i],"_copy_wilcP_Unremoved.RDS",sep=""))
}

colnames(wilcEst)=allCancerTypes
colnames(wilcP)=allCancerTypes
rownames(wilcEst)=rownames(all.copy.data)
rownames(wilcP)=rownames(all.copy.data)

saveRDS(wilcEst,file="MedianDiffEstimate_Copy_Unremoved.RDS")
saveRDS(wilcP,file="WilcoxonRawP_Copy_Unremoved.RDS")

removeRowsList=c(overallNARowIndexes,removeRowsList)
removeRowsList=unique(removeRowsList) #each row can only be removed once
removeRowsList=sort(removeRowsList,decreasing=T) #must remove backwards


all.copy.data=all.copy.data[-removeRowsList,]
wilcEst=wilcEst[-removeRowsList,]
wilcP=wilcP[-removeRowsList,]

saveRDS(wilcEst,"MedianDiffEstimate_Copy_Removed.RDS")
saveRDS(wilcP,"WilcoxonRawP_Copy_Removed.RDS")

colnames(wilcEst)=allCancerTypes
colnames(wilcP)=allCancerTypes
rownames(wilcEst)=rownames(all.copy.data)
rownames(wilcP)=rownames(all.copy.data)

#remove rows that have p-values or estimates we were unable to determine
#wilcEst=wilcEst[-removeRowsList,]
#wilcP=wilcP[-removeRowsList,]

saveRDS(wilcEst,"MedianDiffEstimate_Copy.RDS")
saveRDS(wilcP,"WilcoxonRawP_Copy.RDS")

adj.wilcP.vector=c(wilcP)
adj.wilcP.vector=p.adjust(adj.wilcP.vector,method='fdr')
adj.wilcP=matrix(adj.wilcP.vector,ncol=ncol(wilcP),byrow=F)
rownames(adj.wilcP)=rownames(wilcP)
colnames(adj.wilcP)=colnames(wilcP)
saveRDS(adj.wilcP,"WilcoxonAdjustedP_Copy.RDS")

#Fisher Test - tumors & normals
# combine the p values generated by the wilcoxon test
MetaPs=fisher.method(wilcP,na.rm=T)
write.table(MetaPs,"Meta.Ps_Copy_6.19.14.txt",col.names=NA,sep="\t")

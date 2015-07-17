#Computes foldchanges for methylation with respect to sex
#Jonathan Ma
#6.19.2014

#set to your desired working directory
setwd("C:/Users/admin/Google Drive/ExpressionSexProject/DNA_Copy_Number")

#read in data downloaded from cBioPortal (see the FileDownloader)
#data not available on cbioportal for READ, COAD, and LGG
fileNames=dir(pattern="_methylation.RDS")
allCancerTypes=gsub("([a-z]+)_methylation.RDS","\\1",fileNames)
allCancerTypes=toupper(allCancerTypes)

#since clinical data from cbioportal itself is sparse, need to use master TCGA clinical data file
clinical.data=readRDS('MasterClinical.RDS')
#convert clinical data barcodes to standard format by making them uppercase
clinical.data$barcode <- toupper(clinical.data$barcode)

#wilc.mat: rows=probes, 2 columns (1 for est, 1 for p values)
#each column of wilc.mat is a tissue type
#each row of wilc.mat is a gene in the dna copy number data file (and each data file has the same genes)
all.methyl.data=readRDS(fileNames[1])
colnames(all.methyl.data)=gsub("TCGA.([A-Z0-9]+).([A-Z0-9]+)","TCGA-\\1-\\2",colnames(all.methyl.data))

#initialize matrix to contain fold changes of all cancers
fc.matrix=matrix(NA,nrow=nrow(readRDS(fileNames[1])),ncol=length(fileNames))
rownames(fc.matrix)=rownames(readRDS(fileNames[1]))
colnames(fc.matrix)=allCancerTypes

#loop thru all cancers
for(i in 1:length(fileNames))
{ 
  cat("\nComputing fold changes for ",allCancerTypes[i])
  methyl.data=readRDS(fileNames[i])
  
  #standardize barcodes to matach masterclinical.RDS
  colnames(methyl.data)=gsub("TCGA.([A-Z0-9]+).([A-Z0-9]+)","TCGA-\\1-\\2",colnames(methyl.data))
  #filter to only include barcodes for which we have clinical data
  methyl.data=methyl.data[,colnames(methyl.data)%in%clinical.data$barcode]
  tissue.clinical=clinical.data[clinical.data$barcode%in%colnames(methyl.data),]
  
  #split into male and female
  methyl.male=methyl.data[,tissue.clinical$gender=="male"]
  methyl.female=methyl.data[,tissue.clinical$gender=="female"]
  
  #calculate fold changes
  fc=apply(methyl.male,1,median)/apply(methyl.female,1,median)
  for(j in 1:nrow(methyl.male))
  {
    current.male=methyl.male[j,]
    current.male=current.male[!is.na(current.male)]
    current.female=methyl.female[j,]
    current.female=current.female[!is.na(current.female)]
    fc.matrix[j,i]=median(current.male)/median(current.female)
  }
}

#save the fold changes
saveRDS(fc.matrix,'FoldChanges_Methylation_TissueSpecific.RDS')

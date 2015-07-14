# Computes fold changes for copy number data -- tissue specific
# Jonathan Ma
# 6.19.2014

# Set your working directory appropriately!
#setwd("C:/Users/Jonathan Ma/Google Drive/ExpressionSexProject/Methylation_Cbio")

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
all.copy.data=readRDS(fileNames[1])
colnames(all.copy.data)=gsub("TCGA.([A-Z0-9]+).([A-Z0-9]+)","TCGA-\\1-\\2",colnames(all.copy.data))

fc.matrix=matrix(NA,nrow=nrow(readRDS(fileNames[1])),ncol=length(fileNames))
rownames(fc.matrix)=rownames(readRDS(fileNames[1]))
colnames(fc.matrix)=allCancerTypes
for(i in 1:length(fileNames))
{ 
  cat("\nComputing fold changes for ",allCancerTypes[i])
  copy.data=readRDS(fileNames[i])
  colnames(copy.data)=gsub("TCGA.([A-Z0-9]+).([A-Z0-9]+)","TCGA-\\1-\\2",colnames(copy.data))
  #filter to only include barcodes for which we have clinical data
  copy.data=copy.data[,colnames(copy.data)%in%clinical.data$barcode]
  tissue.clinical=clinical.data[clinical.data$barcode%in%colnames(copy.data),]
  
  copy.male=copy.data[,tissue.clinical$gender=="male"]
  copy.female=copy.data[,tissue.clinical$gender=="female"]
  fc=apply(copy.male,1,median)/apply(copy.female,1,median)
  for(j in 1:nrow(copy.male))
  {
    current.male=copy.male[j,]
    current.male=current.male[!is.na(current.male)]
    current.female=copy.female[j,]
    current.female=current.female[!is.na(current.female)]
    fc.matrix[j,i]=median(current.male)/median(current.female)
  }
}
saveRDS(fc.matrix,'FoldChanges_Copy_TissueSpecific.RDS')

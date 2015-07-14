# This script can be used to download the copy number data via cBioPortal.
# Jonathan Ma

setwd("C:/Users/Jonathan Ma/Google Drive/ExpressionSexProject/DNA_Copy_Number")
#read in a list of all human protein-coding genes
genes=read.table("human_genes.txt",stringsAsFactors=F,sep="\t",header=F,comment.char='',quote='')
genes <- genes[!duplicated(genes$V3),]
gene.synonym <- strsplit(genes$V4,split="|",fixed=T)
names(gene.synonym) <- genes$V3

geneNames = genes$V3


library(cgdsr)
library(BBmisc)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
allCancerStudies=getCancerStudies(mycgds)

#only choose cancer types in mixed
#data not available on cbioportal for READ, COAD, and LGG
cancerStudies=c(
                paad=51,
                skcm=60,
                kich=28,
                sarc=58,
                lihc=35,
                kirp=32,
                gbm=23,
                laml=2,
                blca=8,
                hnsc=26,
                lusc=41,
                luad=38,
                thca=66,
                kirc=31)


for (i in 1:length(cancerStudies))
{
  #choose which cancer study (corresponding to which tissue type) you want
  myCancerStudy=allCancerStudies[cancerStudies[i],1]
  #the category of cases you want (log2 of the copy number)
  myCaseList=getCaseLists(mycgds,myCancerStudy)[5,1]
  #myCaseList=myCaseList$case_ids
  #myCaseList=unlist(myCaseList)
  #myCaseList=strsplit(myCaseList," ")
  #myCaseList=unlist(myCaseList)
  #The name of the category of data (-log2 of the copy number)
  desiredGeneProfile=getGeneticProfiles(mycgds,myCancerStudy)[4,1]
  #the name of the type of cancer
  cancerName=names(cancerStudies)[i]
  cat("Downloading DNA copy number data for: ",cancerName)
 test.data <- t(getProfileData(mycgds,1,desiredGeneProfile,myCaseList,))
  copy.data=t(getProfileData(mycgds,geneNames[1:999],desiredGeneProfile,myCaseList))
  copy.data.matrix <- matrix(NA,nrow=length(geneNames),ncol=ncol(copy.data),dimnames=list(geneNames,colnames(copy.data)))
  good.copy.data <- copy.data[ rownames(copy.data)%in% rownames(copy.data.matrix),] 
 bad.genenames <- geneNames[1:999][ !geneNames[1:999] %in% rownames(copy.data)]
  copy.data.matrix[rownames(good.copy.data),colnames(good.copy.data)] <- good.copy.data
 for(j in 1:length(bad.genenames)){
   test.data <- t(getProfileData(mycgds,bad.genenames[j],desiredGeneProfile,myCaseList))
   if(!all(is.na(test.data[1,]))){
     copy.data.matrix[bad.genenames[j],colnames(test.data)] <- test.data
   }
 }
 
  
  genechunks <- chunk(geneNames,chunk.size=1000)
  for(j in 1:length(genechunks))
  {
    cat("\nCurrently downloading chunk ",j)
    tgeneNames <- geneNames[genechunks[[j]]]
    temp.data <- t(getProfileData(mycgds,tgeneNames,desiredGeneProfile,myCaseList))
    good.copy.data <- temp.data[ rownames(temp.data)%in% rownames(copy.data.matrix),]
    bad.genenames <- tgeneNames[ tgeneNames %in% rownames(copy.data)]
    copy.data.matrix[rownames(good.copy.data),colnames(good.copy.data)] <- good.copy.data
    if(length(bad.genenames)>0)
    {
    for(k in 1:length(bad.genenames)){
      temp.data <- t(getProfileData(mycgds,bad.genenames[k],desiredGeneProfile,myCaseList))
      if(!all(is.na(temp.data[1,]))){
        copy.data.matrix[bad.genenames[k],colnames(temp.data)] <- temp.data
      }
    }
    }
  }
  
  #get clinical data (clinical parameters, i.e. sex, age, tissue type, etc.)
  saveRDS(copy.data,file=paste0(cancerName,"_DNA_Copy_log2.RDS",sep=""))
  clinical.data=getClinicalData(mycgds,myCaseList)
  cat("Downloading clinical data for: ",cancerName)
  saveRDS(clinical.data,file=paste0(cancerName,"_DNA_Copy_Clinical.RDS",sep=""))
}

#Automatically downloads methylation data from cbioportal
#Jonathan Ma
#6.18.2014

#set to your own working directory
setwd("C:/Users/Jonathan Ma/Google Drive/ExpressionSexProject/Methylation_Cbio")

#read in a table of desired genes (HUGO id's)
#designed so that you can select a subset of cBioPortal genes if you want to. If not, just list out a list of all available HUGO human genes.
genes=read.table("human_genes.txt",stringsAsFactors=F,sep="\t",header=F,comment.char='',quote='')
genes <- genes[!duplicated(genes$V3),]
gene.synonym <- strsplit(genes$V4,split="|",fixed=T)
names(gene.synonym) <- genes$V3
geneNames = genes$V3

#load libraries as well as CGDS
library(cgdsr)
library(BBmisc)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
allCancerStudies=getCancerStudies(mycgds)

#select your desired cancer types
#cancer type -> numeric correlations can be obtained from cBioPortal documentation
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

#loop thru all cancer studies (each study = 1 cancer)
for (i in 1:length(cancerStudies))
{
  myCancerStudy=allCancerStudies[cancerStudies[i],1]
  myCaseList=getCaseLists(mycgds,myCancerStudy)[6,1]
  desiredGeneProfile=getGeneticProfiles(mycgds,myCancerStudy)[5,1]
  cancerName=names(cancerStudies)[i]
  cat("\nDownloading methylation data for: ",cancerName)
  
  #note: download in chunks to avoid overflowing RAM
  genechunks <- chunk(geneNames,chunk.size=100)
  
  #get profile data
  methyl.data=t(getProfileData(mycgds,geneNames[1:99],desiredGeneProfile,myCaseList))
  methyl.data.matrix <- matrix(NA,nrow=length(geneNames),ncol=ncol(methyl.data),dimnames=list(geneNames,colnames(methyl.data)))
  good.methyl.data <- methyl.data[ rownames(methyl.data)%in% rownames(methyl.data.matrix),] 
  bad.genenames <- geneNames[1:99][ !geneNames[1:99] %in% rownames(methyl.data)]
  methyl.data.matrix[rownames(good.methyl.data),colnames(good.methyl.data)] <- good.methyl.data
  for(j in 1:length(bad.genenames)){
    test.data <- t(getProfileData(mycgds,bad.genenames[j],desiredGeneProfile,myCaseList))
    if(!all(is.na(test.data[1,]))){
      methyl.data.matrix[bad.genenames[j],colnames(test.data)] <- test.data
    }
  }
  
  #download chunk by chunk and add to the data matrix
  for(j in 1:length(genechunks))
  {
    cat("\nCurrently downloading chunk ",j)
    tgeneNames <- unlist(genechunks[j])
    temp.data <- t(getProfileData(mycgds,tgeneNames,desiredGeneProfile,myCaseList))
    good.methyl.data <- temp.data[ rownames(temp.data)%in% rownames(methyl.data.matrix),]
    bad.genenames <- tgeneNames[ tgeneNames %in% rownames(methyl.data)]
    methyl.data.matrix[rownames(good.methyl.data),colnames(good.methyl.data)] <- good.methyl.data
    if(length(bad.genenames)>0)
    {
      for(k in 1:length(bad.genenames)){
        temp.data <- t(getProfileData(mycgds,bad.genenames[k],desiredGeneProfile,myCaseList))
        if(!all(is.na(temp.data[1,]))){
          methyl.data.matrix[tgeneNames[k],colnames(temp.data)] <- temp.data
        }
      }
    }
  }
  
  #get clinical data (clinical parameters, i.e. sex, age, tissue type, etc.)
  #you can save clinical data if you want; we did not
  clinical.data=getClinicalData(mycgds,myCaseList)
  
  #save methylation data matrix for the current cancer to your desired working directory
  saveRDS(methyl.data.matrix,file=paste0(cancerName,"_methylation.RDS",sep=""))
}

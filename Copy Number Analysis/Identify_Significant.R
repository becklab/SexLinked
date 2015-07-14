# This script is to be used after the complete analysis.
# Jonathan Ma and Sadhika Malladi

# The first section simply finds the genes that are significant (p < 0.05) and notes them.
# The second section attempts to find an overlap of genes found significant by copy number and expression analyses.
# The third section notes how many genes were found significant by copy number analysis in each cancer
# The fourth section has code converting IDs -- this code is useful throughout the project!
# The final section was used to identify potential candidates for functional studies (validation)

# Set your working directory appropriately!
setwd('~/googledrive/ExpressionSexProject/DNA_Copy_Number/')

# This library contains information for converting gene identifier types
library(org.Hs.eg.db)


# Find significant Copy Number genes --------------------------------------

pValsT<-read.table(file='Meta.Ps_Copy_6.21.14.txt',header=T,stringsAsFactors=F)
rownames(pValsT)<-pValsT[,1]
origGenes<-rownames(pValsT)
pValsT<-pValsT[,-1]
pValsT<-pValsT[,'p.adj']

sigPValsT<-pValsT<=0.05
numSigT<-sum(as.numeric(sigPValsT))

sigGenePsT<-matrix(NA,nrow=numSigT,ncol=2)
sigGenePsT[,1]<-pValsT[sigPValsT]
colnames(sigGenePsT)<-c('p-values','on a sex chromosome')
rownames(sigGenePsT)<-origGenes[sigPValsT==T]
genes<-rownames(sigGenePsT)

mapped_entrez_ids<-as.list(org.Hs.egALIAS2EG)
mapped_entrez_ids<-mapped_entrez_ids[genes]
mapped_entrez_ids <- sapply(mapped_entrez_ids,"[",1)
rownames(sigGenePsT)<-mapped_entrez_ids

mapped_genes<-mappedkeys(org.Hs.egCHR)
genes_list<-as.list(org.Hs.egCHR[mapped_genes])
genes_list<-genes_list[mapped_entrez_ids]
for (i in mapped_entrez_ids) {
  if (genes_list[i]=='X'|| genes_list[i]=='Y')
    sigGenePsT[i,2]<-T
  else
    sigGenePsT[i,2]<-F
}

sigGenePsT[,2][is.na(sigGenePsT[,2])]<-0

entrez<-rownames(sigGenePsT)
entrez<-unname(entrez)
write(entrez,'SignificantGenesCopy.txt')

saveRDS(sigGenePsT,file='SignificantEntrezIDs.RDS')

rownames(sigGenePsT)<-genes
saveRDS(sigGenePsT,file='SignificantGenes.RDS')


# Match with RNA seq ------------------------------------------------------

rnaseq<-readRDS(file='~/googledrive/ExpressionSexProject/Expression_Signatures/Tumors_Significant.RDS')
CN<-readRDS(file='SignificantEntrezIDs.RDS')

CN.genes<-rownames(CN)
rna.genes<-rownames(rnaseq)
matched.genes<-CN.genes%in%rna.genes

matchedPs<-CN[matched.genes,]

matched.entrez<-unname(rownames(matchedPs))
write(matched.entrez,'MatchedRNASeqCopy.txt')
write(matched.entrez[1:3000],'MatchedRNASeqCopy_Trunc.txt')

saveRDS(matchedPs,file='MatchedRNASeqSignificantEntrezIDs_WithNormals.RDS')

# Cancer Type Counts ----------------------------------------------------------
# don't need to run previous code chunks
# careful with the working directories that are set

library(xlsx)

# number of significant genes per cancer type
wilcP<-readRDS(file='WilcoxonRawP_Copy.RDS')
wilcP<-wilcP<=0.05
cancerTypes<-colnames(wilcP)
cancerXnumSig<-matrix(NA,nrow=length(cancerTypes),ncol=1,dimnames=list(cancerTypes,'Number of Significant Genes'))
cancerCount<-apply(wilcP,2,sum)
cancerXnumSig[,1]<-cancerCount

# number of patients per cancer type
cancerXnumPat<-matrix(NA,nrow=length(cancerTypes),ncol=1,dimnames=list(cancerTypes,'Number of Patients'))
setwd('C:/Users/Sadhika/Google Drive/ExpressionSexProject/DNA_Copy_Number/Data')
files<-dir()
setwd('C:/Users/Sadhika/Google Drive/ExpressionSexProject/DNA_Copy_Number/')
for (i in files) {
  temp<-readRDS(file=i)
  currentCancer<-toupper(gsub('([A-Za-z])_copy_v2.RDS','\\1',i))
  numPat<-length(colnames(temp))
  cancerXnumPat[currentCancer,1]<-numPat
}

#write to excel sheet
wb<-loadWorkbook('CopyNumberStats.xlsx')
createSheet(wb,sheetName='Number of Significant Genes')
createSheet(wb,sheetName='Number of Patients')
saveWorkbook(wb,'CopyNumberStats.xlsx')
write.xlsx(cancerXnumSig, 'CopyNumberStats.xlsx', sheetName="Number of Significant Genes")
write.xlsx(cancerXnumPat,'CopyNumberStats.xlsx',sheetName='Number of Patients')

# Convert Hugo Gene Symbol to EntrezID ------------------------------------

#setwd('~/googledrive/ExpressionSexProject/DNA_Copy_Number/Data')
library(org.Hs.eg.db)
mapped_entrez_ids<-as.list(org.Hs.egALIAS2EG)
files<-dir(pattern='_copy_v2.RDS')

for (i in files) {
  #setwd('~/googledrive/ExpressionSexProject/DNA_Copy_Number/Data')
  current<-readRDS(i)
  origGenes<-rownames(current)
  mapped_entrez_idsN<-mapped_entrez_ids[origGenes]
  mapped_entrez_idsN <- sapply(mapped_entrez_idsN,"[",1)
  rownames(current)<-mapped_entrez_idsN
  cancer<-toupper(gsub('([a-z])_copy_v2.RDS','\\1',i))
  #setwd('~/googledrive/ExpressionSexProject/DNA_Copy_Number/')
  fileName<-gsub(' ','',paste(cancer,'_copy_v2_entrezIDs.RDS'))
  saveRDS(current,fileName)
}

# Adjusted Ps by Cancer Type ----------------------------------------------
library(org.Hs.eg.db)
mapped_entrez_ids<-as.list(org.Hs.egALIAS2EG)

wilcT<-readRDS('WilcoxonAdjustedP_Copy.RDS')
genes<-rownames(wilcT)

for (i in colnames(wilcT)) {
  Ps<-wilcT[,i]
  Ps[is.nan(Ps)]<-1
  sortedPs<-sort(Ps,decreasing=F)
  genes<-mapped_entrez_ids[names(sortedPs)]
  names(sortedPs) <- sapply(genes,"[",1)
  tempMat<-matrix(NA,nrow=length(sortedPs),ncol=1,dimnames=list(names(sortedPs),'p-value'))
  tempMat[,1]<-sortedPs
  fileName<-gsub(' ','',paste(i,'_SortedAdjustedP_Tumors.RDS'))
  saveRDS(tempMat,fileName)
}


# Functional Analysis Prep ------------------------------------------------
setwd('~/googledrive/ExpressionSexProject/DNA_Copy_Number/')
library(org.Hs.eg.db)

matchedCN<-readRDS('MatchedRNASeqSignificantEntrezIDs.RDS')
Ps<-matchedCN[,1]
sortedPs<-sort(Ps,decreasing=F)
rownames(matchedCN)<-names(sortedPs)
matchedCN[,1]<-sortedPs

fullCN<-readRDS('WilcoxonAdjustedP_Copy.RDS')
mapped_entrez_ids<-as.list(org.Hs.egALIAS2EG)
hugos<-rownames(fullCN)
mapped_entrez_ids<-mapped_entrez_ids[hugos]
mapped_entrez_ids <- sapply(mapped_entrez_ids,"[",1)
rownames(fullCN)<-mapped_entrez_ids

fullRNA<-readRDS('~/googledrive/ExpressionSexProject/Expression_Signatures/WilcoxonAdjustedP_Tumors_6.24.RDS')
rownames(fullRNA)<-gsub(".+\\|([0-9]+)","\\1",rownames(fullRNA))

interestingGenes<-rownames(matchedCN)[1:13]

filteredCN<-fullCN[interestingGenes,]
filteredRNA<-fullRNA[interestingGenes,]

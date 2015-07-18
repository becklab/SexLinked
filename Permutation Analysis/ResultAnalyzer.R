#Code intended to validate LincsCloud results through Permutation Analysis, inspired by and modeled after the PMASE approach of Liu et al. (2011)
#Intended to analyze the raw results output of LincsCloud C3 permuted queries en masse and output the resulting empirical p-values of sex-disparate
#connectivity score.
#Author: Jonathan Ma
#Created: 12/31/2014

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
library(R.matlab)
library(calibrate)
library(colbycol)
library(xlsx)
library(base)

#Cancers used in permutation analysis
cancers=c("HNSC","KICH","KIRC","LIHC","LUAD","LUSC","THCA","Meta")

numberofpermutations=1000 #number of permutations we did
rfn="result_NA_summly_n7503.txt" #file name that LincsCloud C3 assigns to the result file

for(i in 1:length(cancers))
{
  cc=cancers[i]
  
  #set to your appropriate directory
  setwd(paste0("C:/Users/Jonathan Ma/Google Drive/ExpressionSexProject/Results/LincsCloud/Results/",cc))
  
  #read in the actual (non-permuted) connectivity scores which were obtained through LincsCloud connectivity map analysis of resistance and
  #sensitivity signatures.
  diff=read.table("EnrichmentDifferences.txt",row.names=1,stringsAsFactors=F)
  
  pvals=matrix(data=0,nrow=nrow(diff),ncol=2) #col1 = exceeded actual, col2=p-value
  
  #iterate through the raw results of each permutation, computing empirical p-values for each drug based on the number of permutations whose 
  #male-female connectivity score difference is greater in magnitude than the actual (non-permuted) computed value
  
  for(j in 1:numberofpermutations)
  {
    male.exists=T
    female.exists=T
    cat("\ncomputing ",cc," permutation ",j)
    if(file.exists(paste0("C:/Users/Jonathan Ma/Documents/Research/LincsCloud_C3/Permuted_Results_SAM/",cc,"/male",j,"/summary/1/",rfn)))
    {
      setwd(paste0("C:/Users/Jonathan Ma/Documents/Research/LincsCloud_C3/Permuted_Results_SAM/",cc,"/male",j,"/summary/1"))
      mpr=read.table(rfn,stringsAsFactors=F,fill=T,sep="\t") #ignore warnings
      
      if(j==1) firstfile=mpr
    } else {
      male.exists=F
    }
    
    if(file.exists(paste0("C:/Users/Jonathan Ma/Documents/Research/LincsCloud_C3/Permuted_Results_SAM/",cc,"/female",j,"/summary/1/",rfn)))
    {
      setwd(paste0("C:/Users/Jonathan Ma/Documents/Research/LincsCloud_C3/Permuted_Results_SAM/",cc,"/female",j,"/summary/1"))
      fpr=read.table(rfn,stringsAsFactors=F,fill=T,sep="\t") #ignore warnings
    } else {
      female.exists=F
    }
    
    if(male.exists)
    {
      rownames(mpr)=mpr[,1]
      colnames(mpr)=mpr[1,]
      mpr=mpr[2:nrow(mpr),]
      diff.mpr=mpr[rownames(diff),'mean_rankpt_4']
      diff.mpr.mod=diff.mpr
      diff.mpr.mod[is.na(diff.mpr.mod)]=0
      diff.mpr.mod=as.numeric(diff.mpr.mod)
    } else {
      diff.mpr.mod=rep.int(0,times=nrow(pvals))
    }
    
    if(female.exists)
    {
      rownames(fpr)=fpr[,1]
      colnames(fpr)=fpr[1,]
      fpr=fpr[2:nrow(fpr),]
      diff.fpr=fpr[rownames(diff),'mean_rankpt_4']
      diff.fpr.mod=diff.fpr
      diff.fpr.mod[is.na(diff.fpr.mod)]=0
      diff.fpr.mod=as.numeric(diff.fpr.mod)
    } else{
      diff.fpr.mod=rep.int(0,times=nrow(pvals))
    }
    
    pdf=data.frame(diff.mpr.mod,diff.fpr.mod,abs(diff.mpr.mod-diff.fpr.mod),row.names=rownames(diff))
    colnames(pdf)=c("mscore","fscore","perm.diff")
    
    
    pvals[abs(pdf$perm.diff)>abs(diff[,3]),1]=pvals[abs(pdf$perm.diff)>abs(diff[,3]),1]+1
  }
  
  #label the perturbagen type (i.e. upregulation or knockdown)
  rownames(firstfile)=firstfile[,1]
  colnames(firstfile)=firstfile[1,]
  firstfile=firstfile[2:nrow(fpr),]
  firstfile=firstfile[rownames(pdf),]
  pert.type=firstfile$pert_type
  
  #compute empirical p-valiue
  pvals[,2]=pvals[,1]/1000
  setwd(paste0("C:/Users/Jonathan Ma/Documents/Research/LincsCloud_C3/Permuted_Pvals_SAM"))
  pvals.df=data.frame(pvals)
  rownames(pvals.df)=rownames(diff)
  colnames(pvals.df)=c("exceeded actual","pval")
  
  #adjust p-values using fdr method, yielding Benjamini-hochberg values
  adj.pval=p.adjust(pvals.df$pval,method='fdr')
  pvals.df=cbind(pvals.df,adj.pval)
  
  #add the actual (non permuted) connectivity scores
  pvals.df=cbind(pvals.df,diff$males.comp,diff$females.comp,diff$diff,diff$Common.Drug.Name)
  colnames(pvals.df)=c("exceeded.actual","pval","adj.pval","actual.male.cs","actual.female.cs","actual.diff","common.drug.name")
  pvals.df=cbind(pvals.df,pert.type)
  
  ranked.pvals.df=pvals.df[order(pvals.df$adj.pval),]
  
  #save results
  write.table(pvals.df,file=paste0(cc,"_permutedpvals.txt"),row.names=T,col.names=T)
  write.xlsx(x=pvals.df,file=paste0(cc,"_permutedpvals.xlsx"),sheetName="sheet1",col.names=T,row.names=T)
  write.table(ranked.pvals.df,file=paste0(cc,"_ranked_permutedpvals.txt"),row.names=T,col.names=T)
  write.xlsx(x=pvals.df,file=paste0(cc,"_ranked_permutedpvals.xlsx"),sheetName="sheet1",col.names=T,row.names=T)
}

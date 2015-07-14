# Analyze LincsCloud results: generate plots seen in paper
# Plots have connectivity scores for males vs. connectivity scores for females.
# The "three corners" are of interest: effective in males only, effective in females only, and effective in both
# Sadhika Malladi

# Change working directory appropriately
setwd('~/googledrive/ExpressionSexProject/Results/Lincscloud/Results/')
library(calibrate)

cancers <- dir()

for (i in cancers) {
  setwd(paste0('~/googledrive/ExpressionSexProject/Results/Lincscloud/Results/',i))
  
  if (i =='Pan-Cancer') {
    i<-'Meta'
  }
  
  #ignore warnings
  females<-read.table(paste0('Females_',i,'.txt'),stringsAsFactors=F,row.names=1,header=T,fill=T,sep='\t')
  males<-read.table(paste0('Males_',i,'.txt'),stringsAsFactors=F,row.names=1,header=T,fill=T,sep='\t')
  
  # compare drugs with info in both males and females, isolate mean rankpt 4
  femalesComp<-females[rownames(females)%in%rownames(males),]
  malesComp<-males[rownames(femalesComp),]
  females.comp<-femalesComp[,'mean_rankpt_4']
  names(females.comp)<-rownames(femalesComp)
  males.comp<-malesComp[,'mean_rankpt_4']
  names(males.comp)<-rownames(malesComp)
  
  # handle NAs
  males.comp[is.na(males.comp)]<- 0
  females.comp[is.na(females.comp)]<-0
  
  # get common drug names
  drugnames<-females[,'pert_iname'][rownames(females)%in%rownames(males)]
  drugnames<-drugnames[!is.na(males.comp)]
  
  # calculate difference
  diff<-abs(males.comp-females.comp)
  diff<-as.numeric(diff)
  names(diff)<-names(females.comp)
  sum<-males.comp+females.comp
  sum<-as.numeric(sum)
  names(sum)<-names(diff)
  sum<-sort(sum,decreasing=T)
  
  mat<-cbind(males.comp,females.comp,diff)
  mat<-cbind(mat,drugnames)
  perts<-femalesComp[,'pert_type']
  perts<-perts[perts!='']
  mat<-cbind(mat,perts)
  colnames(mat)[4:5]<-c('Common Drug Name','Perturbation Type')
  mat<-mat[order(as.numeric(mat[,'diff']),decreasing=T),]
  
  # plot male vs female, look for deviance off diagonal
  femalesComp<-females[rownames(females)%in%rownames(males),]
  malesComp<-males[rownames(femalesComp),]
  
  femaleEnrichment<-as.numeric(mat[,'females.comp'])>0
  maleEnrichment<-as.numeric(mat[,'males.comp'])>0
  fullEnrichment<-cbind(femaleEnrichment,maleEnrichment,mat[,'diff'])
  rownames(fullEnrichment)<-rownames(mat)
  
  # find male and female enriched drugs, just in male, and then just in female
  femaleEnriched<-fullEnrichment[fullEnrichment[,1]==T & fullEnrichment[,2]==F,]
  maleEnriched<-fullEnrichment[fullEnrichment[,1]==F & fullEnrichment[,2]==T,]
  
  # highlight known drugs -- these were manually input
  if (i == 'HNSC') {
    knownDrugs<-c('methotrexate','fluorouracil','bleomycin','cetuximab','cisplatin','docetaxel')
  } else if (i=='KIRC') {
    knownDrugs<-c('everolimus','aldesleukin','bevacizumab','axitinib','sorafenib','sunitinib','temsirolimus','pazopanib hydrochloride')
  }  else if (i == 'LIHC') {
    knownDrugs<-c('sorafenib')
  }  else if (i =='LUAD' | i=='LUSC') {
    knownDrugs<-c('carboplatin-taxol','gemcitabine-cisplatin')
  }  else if (i=='THCA') {
    knownDrugs<-c('doxorubicin hydrochloride','cabozantinib-s-malate','vandetanib','sorafenib')
  }  else {
    knownDrugs<-c()
  }
  
  # check if any known drugs are found in linscloud data
  found<-sum(mat[,'Common Drug Name']%in%knownDrugs)
  if (found != 0) {
    foundDrug<-mat[mat[,'Common Drug Name']%in%knownDrugs,]
    foundDrug<-matrix(foundDrug,nrow=found,ncol=5,dimnames=list(rownames(mat)[mat[,'Common Drug Name']%in%knownDrugs],colnames(mat)))
  }
  
  if (i == 'Meta') {
    title<-'Comparison of Mean Connectivity Scores of Perturbagens in Males versus Females Across All Cancers Considered'
  } else {
    title<-paste('Comparison of Mean Connectivity Scores of Perturbagens in Males versus Females in',i)
  }
  
  # determine scale for plot
  minX<-min(as.numeric(femalesComp[,'mean_rankpt_4']),na.rm=T)
  maxX<-max(as.numeric(femalesComp[,'mean_rankpt_4']),na.rm=T)
  maxY<-max(as.numeric(malesComp[,'mean_rankpt_4']),na.rm=T)
  minY<-min(as.numeric(malesComp[,'mean_rankpt_4']),na.rm=T)
  
  # extract perturbation type (e.g., overexpression or knockdown)
  malesGenes<-malesComp[malesComp[,'pert_type']!='trt_cp',]
  malesDrugs<-malesComp[malesComp[,'pert_type']=='trt_cp',]
  femalesGenes<-femalesComp[femalesComp[,'pert_type']!='trt_cp',]
  femalesDrugs<-femalesComp[femalesComp[,'pert_type']=='trt_cp',]
  
  png('ScatterPlot_MalesVsFemales_3CornersAndKnownDrugs_NAIncluded.png',res=120,width=2500,height=1606,units='px')
  par(mar=c(10,5,5,1) + 0.1)
  if (i != 'HNSC') {
    plot(malesDrugs[,'mean_rankpt_4']~femalesDrugs[,'mean_rankpt_4'],xlab='Female Mean Connectivity Score',ylab='Male Mean Connectivity Score',main=title,cex=2,cex.lab=2, cex.axis=2, cex.main=2,col='darkorange')
    par(new=T)
    par(mar=c(10,5,5,1) + 0.1)
    plot(malesGenes[,'mean_rankpt_4']~femalesGenes[,'mean_rankpt_4'],xlab='',ylab='',main='',cex=2,cex.lab=2, cex.axis=2, cex.main=2,axes=F,col='gray23')
  } else {
    plot(malesDrugs[,'mean_rankpt_4']~femalesDrugs[,'mean_rankpt_4'],xlab='Female Mean Connectivity Score',ylab='Male Mean Connectivity Score',main=title,cex=2,cex.lab=2, cex.axis=2, cex.main=2,col='darkorange',xlim=c(45,100),ylim=c(-90,100))
    par(new=T)
    par(mar=c(10,5,5,1) + 0.1)
    plot(malesGenes[,'mean_rankpt_4']~femalesGenes[,'mean_rankpt_4'],xlab='',ylab='',main='',cex=2,cex.lab=2, cex.axis=2, cex.main=2,axes=F,col='gray23',xlim=c(45,100),ylim=c(-90,100))
  }
  abline(coef=c('0','1'),lwd=7)
  points(x=mat[rownames(femaleEnriched)[1:5],'females.comp'],y=mat[rownames(femaleEnriched)[1:5],'males.comp'],col='pink',bg='pink',pch=16,cex=5)
  points(x=mat[rownames(maleEnriched)[1:5],'females.comp'],y=mat[rownames(maleEnriched)[1:5],'males.comp'],col='blue',bg='blue',pch=16,cex=5)
  points(x=mat[names(sum)[1:5],'females.comp'],y=mat[names(sum)[1:5],'males.comp'],col='darkorchid1',bg='darkorchid1',pch=16,cex=5)
  if (found !=0) {
    points(x=foundDrug[,'females.comp'],y=foundDrug[,'males.comp'],col='green',bg='green',pch=16,cex=5) 
  }
  par(mar=c(0, 0, 0, 0))
  dev.off()
  
  # write table to describe the drugs in the 3 corners
  allDrugs<-c(rownames(femaleEnriched)[1:5],rownames(maleEnriched)[1:5],names(sum)[1:5])
  tab<-mat[allDrugs,]
  write.table(tab,'3CornersInformation_NAIncluded.txt')
  
  # write information about known drugs
  if (found != 0) {
    write.table(foundDrug,'KnownDrugsFound_NAIncluded.txt')
  }
  
  # print summary to console
  cat(i,' Known Drugs Found:',found,': ')
  if (found !=0) {
    cat(rownames(foundDrug))
  }
  
}


# Create table of stats ------------------------------------------------------------
# Creates table with number of drugs returned in each analysis and basic stats (connectivity scores and perturbation type)
# Supplementary to graphs generated above

setwd('~/googledrive/ExpressionSexProject/Results/Lincscloud/Results/')
library(calibrate)

cancers <- dir()

colN<-c('Total Number of Drugs and Genes Returned for both males and females','Number of Drugs returned','Number of genes returned','Number of Drugs with Enrichment Difference of at least 30', 'Number of drugs with min enrichment difference and opposite signs','Number of genes with Enrichment Difference of at least 30', 'Number of genes with min enrichment difference and opposite signs')
mat<-matrix(NA,nrow=length(cancers),ncol=7, dimnames=list(cancers,colN))

for (i in cancers) {
  setwd(paste0('~/googledrive/ExpressionSexProject/Results/Lincscloud/Results/',i))
  
  # ignore warnings!
  females<-read.table(paste0('Females_',i,'.txt'),stringsAsFactors=F,row.names=1,header=T,fill=T,sep='\t')
  males<-read.table(paste0('Males_',i,'.txt'),stringsAsFactors=F,row.names=1,header=T,fill=T,sep='\t')
  
  # comparison
  femalesComp<-females[rownames(females)%in%rownames(males),]
  malesComp<-males[rownames(femalesComp),]
  females.comp<-femalesComp[,'mean_rankpt_4']
  names(females.comp)<-rownames(femalesComp)
  males.comp<-malesComp[,'mean_rankpt_4']
  names(males.comp)<-rownames(malesComp)
  
  # handle NAs
  males.comp[is.na(males.comp)]<-0
  females.comp[is.na(females.comp)]<-0
  
  # get common drug names
  drugnames<-females[,'pert_iname'][rownames(females)%in%rownames(males)]
  drugnames<-drugnames[!is.na(males.comp)]
  
  # total number of drugs and genes
  mat[i,1]<-dim(malesComp)[1]
  # number of drugs
  drugs<-malesComp[malesComp[,'pert_type']=='trt_cp',]
  mat[i,2]<-dim(drugs)[1]
  # number of genes
  genes<-malesComp[malesComp[,'pert_type']!='trt_cp',]
  mat[i,3]<-dim(genes)[1]
  
  # calculate drugs diff
  drugsdiff<-abs(as.numeric(malesComp[rownames(drugs),'mean_rankpt_4'])-as.numeric(femalesComp[rownames(drugs),'mean_rankpt_4']))
  names(drugsdiff)<-rownames(drugs)
  drugsdiff<-drugsdiff[drugsdiff>30]
  mat[i,4]<-length(drugsdiff)
  drugmat<-cbind(malesComp[rownames(drugs),'mean_rankpt_4'],femalesComp[rownames(drugs),'mean_rankpt_4'])
  dimnames(drugmat)<-list(rownames(drugs),c('Male','Female'))
  drugmat<-drugmat[names(drugsdiff),]
  # find opp signs
  opp<-sum(drugmat[,1]<0 & drugmat[,2]>0)
  opp<-opp + sum(drugmat[,1]>0 & drugmat[,2]<0)
  mat[i,5]<-opp
  
  # calculate genes diff
  genesdiff<-abs(as.numeric(malesComp[rownames(genes),'mean_rankpt_4'])-as.numeric(femalesComp[rownames(genes),'mean_rankpt_4']))
  names(genesdiff)<-rownames(genes)
  genesdiff<-genesdiff[!is.na(genesdiff)]
  genesdiff<-genesdiff[genesdiff>30]
  mat[i,6]<-length(genesdiff)
  genemat<-cbind(malesComp[rownames(genes),'mean_rankpt_4'],femalesComp[rownames(genes),'mean_rankpt_4'])
  dimnames(genemat)<-list(rownames(genes),c('Male','Female'))
  genemat<-genemat[names(genesdiff),]
  # find opp signs
  opp<-sum(genemat[,1]<0 & genemat[,2]>0)
  opp<-opp + sum(genemat[,1]>0 & genemat[,2]<0)
  mat[i,7]<-opp
}

write.table(mat,'~/googledrive/ExpressionSexProject/Results/Lincscloud/Results/Stats_NAIncluded.txt',sep=',')




# This script functions to write plot comparisons of the pathways returned by ingenuity pathway analysis
# as significantly enriched in male and/or female tumor patients.
# Sadhika Malladi

library(XLConnect)
library(calibrate)

# Be sure to change the working directories appropriately!

# Create Plots ------------------------------------------------------------
# Create a plot for each cancer showing the results from ingenuity pathway analysis. Each plot
# has the pathways returned as enriched in male (blue points) and female (red points) and
# the -log10 of the corresponding p-value.

setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/MaleTvNPlots/')
files<-dir(pattern='.xls')
cancer<-gsub('(.*)_male.*','\\1',files)
for (i in cancer) {
  
  # female
  file<-paste0(i,'_female_tvn.xls')
  setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/FemaleTvNPlots/')
  tab.f<-readWorksheet(loadWorkbook(file),sheet=1)
  colnames(tab.f)<-tab.f[1,]
  tab.f<-tab.f[-1,]
  column<-as.numeric(tab.f[,2])
  tab.f<-tab.f[order(-column),]
  
  # male
  file<-paste0(i,'_male_tvn.xls')
  setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/MaleTvNPlots/')
  tab.m<-readWorksheet(loadWorkbook(file),sheet=1)
  colnames(tab.m)<-tab.m[1,]
  tab.m<-tab.m[-1,]
  column<-as.numeric(tab.m[,2])
  tab.m<-tab.m[order(-column),]
  
  # write comparison image
  # uncomment the textxy and legend functions if you want to label the five most enriched pathways on the plot
  # commented right now because pathway names became illegible (tables are written in later part of code)
  setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/MaleFemaleComparison/')
  fileName<-paste0(i,'_malefemale.png')

  png(fileName,res=120,width=2500,height=1606,units='px')
  par(mar=c(10,5,5,1) + 0.1)
  plot(as.numeric(tab.f[,2]),xlab='Pathway Number (Ranked by Significance)',ylab='-log10(p-value)',main=paste('Significance of Pathway Association with Neoplastic',i,'Expression Signature in Males and Females'),col='red',cex=2,cex.lab=2, cex.axis=2, cex.main=2)
#   textxy(tab.f[,2][1:5],tab.f[,1][1:5],tab.f[,1][1:5],col='red',cex=2)
#   legend(x=225,y=8,legend=tab.f[,1][1:5],col=rep('red',5),cex=2,pch=1)
  par(new=T)
  par(mar=c(10,5,5,1) + 0.1)
  plot(tab.m[,2],xlab='',ylab='',main='',col='blue',cex=2,axes=F)
#   textxy(tab.m[,2][1:5],tab.m[,1][1:5],tab.m[,1][1:5],col='blue',cex=2)
#   legend(x=225,y=4,legend=tab.m[,1][1:5],col=rep('blue',5),cex=2,pch=1)
  dev.off()

  # write 5 pathways for each sex
  fileName<-paste0(i,'_malefemale_fivepathways.txt')
  colnames(tab.f)[1:2]<-c('Female Pathway','-log10 p-value')
  colnames(tab.m)[1:2]<-c('Male Pathway','-log10 p-value')
  fullmat<-cbind(tab.f[1:5,1:2],tab.m[1:5,1:2])
  write.table(fullmat,file=fileName,sep='\t')

}


# Write tables for plots --------------------------------------------------
# Write tables that correspond to the plots of the previous section
# Tables include the pathways enriched most in males and females and corresponding p-values

setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/FemaleTvNPlots/')
files<-dir()
cancer<-gsub('(.*)_female.*','\\1',files)

for (i in cancer) {
  # female
  file<-paste0(i,'_female_tvn.xls')
  setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/FemaleTvNPlots/')
  tab.f<-readWorksheet(loadWorkbook(file),sheet=1)
  colnames(tab.f)<-tab.f[1,]
  tab.f<-tab.f[-1,]
  column<-as.numeric(tab.f[,2])
  tab.f<-tab.f[order(-column),]
  
  # male
  file<-paste0(i,'_male_tvn.xls')
  setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/MaleTvNPlots/')
  tab.m<-readWorksheet(loadWorkbook(file),sheet=1)
  colnames(tab.m)<-tab.m[1,]
  tab.m<-tab.m[-1,]
  column<-as.numeric(tab.m[,2])
  tab.m<-tab.m[order(-column),]
  
  # find overlap
  pathways<-unique(c(tab.f[,1][1:5],tab.m[,1][1:5]))
  female<-tab.f[tab.f[,1]%in%pathways,]
  male<-tab.m[tab.m[,1]%in%pathways,]
  names.m<-male[,1]
  names.f<-female[,1]
  male<-male[,2]
  female<-female[,2]
  names(male)<-names.m
  names(female)<-names.f
  
  # create table
  mat<-matrix(NA,nrow=length(pathways),ncol=2,dimnames=list(pathways,c('Male p-value','Female p-value')))
  for (x in rownames(mat)) {
    mat[x,1]<-male[x]
    mat[x,2]<-female[x]
  }
  
  # write table
  setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/MaleFemaleComparison/')
  fileName<-paste0(i,'pVal_Comparison.txt')
  write.table(mat,file=fileName)
}


# Create p-val Plots ------------------------------------------------------
# Plots a comparison of p-values for each pathway in males and females.
# Found to not be very indicative of any patterns. Instead, use the tables and plots generated
# above to draw comparisons manually.

library(calibrate)

setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/FemaleTvNPlots/')
files<-dir()
cancer<-gsub('(.*)_female.*','\\1',files)

for (i in cancer) {
  # female
  file<-paste0(i,'_female_tvn.xls')
  setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/FemaleTvNPlots/')
  tab.f<-readWorksheet(loadWorkbook(file),sheet=1)
  colnames(tab.f)<-tab.f[1,]
  tab.f<-tab.f[-1,]
  column<-as.numeric(tab.f[,2])
  tab.f<-tab.f[order(-column),]
  rownames(tab.f)<-tab.f[,1]
  
  # male
  file<-paste0(i,'_male_tvn.xls')
  setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/PathwayPlots/MaleTvNPlots/')
  tab.m<-readWorksheet(loadWorkbook(file),sheet=1)
  colnames(tab.m)<-tab.m[1,]
  tab.m<-tab.m[-1,]
  column<-as.numeric(tab.m[,2])
  tab.m<-tab.m[order(-column),]
  rownames(tab.m)<-tab.m[,1]
  
  # get pathways and p-values in one plot
  female<-tab.f[tab.f[,1]%in%tab.m[,1],]
  male<-tab.m[rownames(female),]
  male[,2]<-10^-as.numeric(male[,2])
  female[,2]<-10^-as.numeric(female[,2])
  mat<-cbind(female[,1],male[,2],female[,2])
  colnames(mat)<-c('Pathway','Male p-value','Female p-value')
  
  # points to label:
  # 1) greatest difference in p-value
  # 2) high in males, low in females
  # 3) high in females, low in males
  # 4) low in both 
  diff<-abs(as.numeric(mat[,2])-as.numeric(mat[,3]))
  names(diff)<-mat[,1]
  diff<-sort(diff,decreasing=T)
  enrichMat<-mat
  # 2
  rownames(enrichMat)<-enrichMat[,1]
  enrichMat<-enrichMat[names(diff),]
  maleEnrich<-enrichMat[enrichMat[,2]<enrichMat[,3],]
  malediff<-abs(as.numeric(maleEnrich[,2])-as.numeric(maleEnrich[,3]))
  maleEnrich<-cbind(maleEnrich,malediff)
  maleEnrich<-maleEnrich[order(-as.numeric(maleEnrich[,4])),]
  # 3
  femaleEnrich<-enrichMat[enrichMat[,2]>enrichMat[,3],]
  femalediff<-abs(as.numeric(femaleEnrich[,2])-as.numeric(femaleEnrich[,3]))
  femaleEnrich<-cbind(femaleEnrich,femalediff)
  femaleEnrich<-femaleEnrich[order(-as.numeric(femaleEnrich[,4])),]
  # 4
  diff<-sort(diff,decreasing=F)
  enrichMat<-enrichMat[names(diff),]
  bothEnrich<-enrichMat[enrichMat[,2]<0.001 & enrichMat[,2]<0.001,]
  
  # plot
  file<-paste0(i,'_pValPlot.png')
  title<-paste0('Comparison of P-values for Each Pathway in Males versus Females with ',i)
  png(file,res=120,width=2500,height=1606,units='px')
  par(mar=c(10,5,5,1) + 0.1)
  plot(mat[,2]~mat[,3],xlab='Female P-value',ylab='Male P-value',main=title,cex=2,cex.lab=2, cex.axis=2, cex.main=2)
  abline(coef=c('0','1'))
  par(mar=c(0, 0, 0, 0))
  dev.off()
}

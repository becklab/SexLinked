# This script formats and generates excel files that function as input to 
# ingenuity pathway analysis.
# Sadhika Malladi

# Change working directory appropriately!
setwd('~/googledrive/ExpressionSexProject/Connectivity_Map_v2/')

library(xlsx)

ps.f<-readRDS('adj.female.wilcP.RDS')
ps.m<-readRDS('adj.male.wilcP.RDS')
fc.f<-readRDS('female.fc.matrix.RDS')
fc.m<-readRDS('male.fc.matrix.RDS')
fc.f[is.na(fc.f)]<-1
fc.m[is.na(fc.m)]<-1
rownames(ps.f)<-gsub('.*\\|([0-9]+)','\\1',rownames(ps.f))
rownames(ps.m)<-rownames(ps.f)
rownames(fc.m)<-rownames(ps.f)
rownames(fc.f)<-rownames(ps.f)
# fc < 1 indicates greater enrichment in normals

for (i in colnames(ps.f)) {
   cat('\n',i) 
  
   # female
   fem<-cbind(ps.f[,i],fc.f[,i])
   colnames(fem)<-c('P-values','Fold Changes')
   setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/Female T v N/')
   file<-paste0(i,'.xlsx')
   write.xlsx2(fem,file=file)
   cat('\n',i,'female done') 
  
   # male
   male<-cbind(ps.m[,i],fc.m[,i])
   colnames(male)<-c('P-values','Fold Changes')
   setwd('~/googledrive/ExpressionSexProject/Results/Pathway Analysis/Male T v N/')
   file<-paste0(i,'.xlsx')
   write.xlsx2(male,file)
   cat('\n',i,'male done') 
}
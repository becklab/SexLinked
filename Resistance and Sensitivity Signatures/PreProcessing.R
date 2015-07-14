# Sort and save significant genes to prepare for connectivity map analysis
# Sadhika Malladi

# Choose the appropriate working directory
setwd('~/googledrive/ExpressionSexProject/Connectivity_Map_v2/')

maleAdj<-readRDS('adj.male.wilcP.RDS')
femaleAdj<-readRDS('adj.female.wilcP.RDS')

foldChangePVal<-function(exp,FC,fileName) {
  exp<-exp[exp<0.05]
  exp<-exp[names(exp)!='?']
  exp<-sort(exp,decreasing=F)
  FC.up<-FC[FC>1]
  FC.down<-FC[FC<1]
  exp.up<-exp[names(exp)%in%names(FC.up)]
  exp.up<-exp.up[1:250]
  exp.down<-exp[names(exp)%in%names(FC.down)]
  exp.down<-exp.down[1:250]
  fn<-gsub(' ','',paste(fileName,'_up.txt'))
  write(names(exp.up),fn)
  fn<-gsub(' ','',paste(fileName,'_down.txt'))
  write(names(exp.down),fn)
}

maleAdj.hugos<-maleAdj
rownames(maleAdj.hugos)<-gsub("([?A-Z0-9]+)\\|[0-9]+","\\1",rownames(maleAdj.hugos))
maleAdj.entrez<-maleAdj
rownames(maleAdj.entrez)<-gsub(".+\\|([0-9]+)","\\1",rownames(maleAdj.entrez))

femaleAdj.hugos<-femaleAdj
rownames(femaleAdj.hugos)<-gsub("([?A-Z0-9]+)\\|[0-9]+","\\1",rownames(femaleAdj.hugos))
femaleAdj.entrez<-femaleAdj
rownames(femaleAdj.entrez)<-gsub(".+\\|([0-9]+)","\\1",rownames(femaleAdj.entrez))

maleFC<-readRDS('male.fc.matrix.RDS')
maleFC<-maleFC[complete.cases(maleFC),]
maleFC.hugos<-maleFC
rownames(maleFC.hugos)<-gsub("([?A-Z0-9]+)\\|[0-9]+","\\1",rownames(maleFC.hugos))
maleFC.entrez<-maleFC
rownames(maleFC.entrez)<-gsub(".+\\|([0-9]+)","\\1",rownames(maleFC.entrez))

femaleFC<-readRDS('female.fc.matrix.RDS')
femaleFC<-femaleFC[complete.cases(femaleFC),]
femaleFC.hugos<-femaleFC
rownames(femaleFC.hugos)<-gsub("([?A-Z0-9]+)\\|[0-9]+","\\1",rownames(femaleFC.hugos))
femaleFC.entrez<-femaleFC
rownames(femaleFC.entrez)<-gsub(".+\\|([0-9]+)","\\1",rownames(femaleFC.entrez))

# sort by p-value and save txt files
for (i in colnames(maleAdj)) {
  foldChangePVal(maleAdj.hugos[,i],maleFC.hugos[,i],gsub(' ','',paste(i,'male_sighugos')))
  foldChangePVal(femaleAdj.hugos[,i],femaleFC.hugos[,i],gsub(' ','',paste(i,'female_sighugos')))
  foldChangePVal(maleAdj.entrez[,i],maleFC.entrez[,i],gsub(' ','',paste(i,'male_sigentrez')))
  foldChangePVal(femaleAdj.entrez[,i],femaleFC.entrez[,i],gsub(' ','',paste(i,'female_sigentrez')))
}

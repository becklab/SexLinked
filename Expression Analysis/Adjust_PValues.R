Tumor.wilcP=readRDS('WilcoxonRawP_Tumors_6.23.RDS')
adj.Tumor.wilcP=Tumor.wilcP
for(i in 1:ncol(Tumor.wilcP))
{
  pvals=Tumor.wilcP[,i]
  adj.pvals=p.adjust(pvals,method='fdr')
  adj.Tumor.wilcP[,i]=adj.pvals
}
saveRDS(adj.Tumor.wilcP,file='WilcoxonAdjustedP_Tumors_6.24.RDS')

Normal.wilcP=readRDS('WilcoxonRawP_Normals_6.23.RDS')
adj.Normal.wilcP=Normal.wilcP
for(i in 1:ncol(Normal.wilcP))
{
  pvals=Normal.wilcP[,i]
  adj.pvals=p.adjust(pvals,method='fdr')
  adj.Normal.wilcP[,i]=adj.pvals
}
saveRDS(adj.Normal.wilcP,file='WilcoxonAdjustedP_Normals_6.24.RDS')

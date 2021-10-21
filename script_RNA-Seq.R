#full analytical pipeline
#WYLai_2021_1021

rm(list=ls())
library(edgeR)
library(pheatmap)
library(biomaRt)
library(pathview)
library(RColorBrewer)
library(data.table)
library(tidyverse)
library(gplots)
library(magrittr)
library(dplyr)
library(ExpressionNormalizationWorkflow)
library(car)
library(scales)
library(colorspace)
setwd("/Volumes/cluster/Wei-Yun/")
####import basic function####
cont_table=function(query,background,classifyer){
  p1=length(Reduce(intersect,list(query,background,classifyer)))
  q1=length(intersect(query,background))-p1
  p0=length(setdiff(intersect(background,classifyer),intersect(query,classifyer)))
  q0=length(setdiff(background,query))-p0
  return(matrix(c(p1,p0,q1,q0),2,2))
}
genomic_enrichment=function(query,background,ann_summary,windowsize,steps,chr){
  X0=0
  X=windowsize
  S=steps
  pos=c()
  counts=c()
  fold_enrichment=c()
  p.val=c()
  x=ann_summary[ann_summary[,1]%in%chr,]
  while (X0<max(x[,c(2,3)])){
    gene_within=x[x[,2]<X0+X/2&x[,2]>X0-X/2,5]
    p=c()
    fe=c()
    for (j in 1:length(query)){
      p=cbind(p,fisher.test(cont_table(query[[j]],background,gene_within),alternative = "greater")$p.value) 
      fe=cbind(fe,(length(intersect(gene_within,query[[j]]))/length(gene_within))/(length(query[[j]])/length(background)))
    }
    p.val=rbind(p.val,p)
    fold_enrichment=rbind(fold_enrichment,fe)
    colnames(p.val)=names(query)
    colnames(fold_enrichment)=names(query)
    counts=c(counts,length(x[x[,2]<X0+X/2&x[,2]>X0-X/2,5]))
    pos=c(pos,X0)
    X0=X0+S
  }
  out_p=data.frame("mid_pos"=pos,"total_gene_counts"=counts,p.val)
  out_fe=data.frame("mid_pos"=pos,"total_gene_counts"=counts,fold_enrichment)
  out=list("p"=out_p,"fe"=out_fe)
  return(out)
}
ID_converter=function(ID,db,attributes,filters){
  getBM(attributes=attributes,filters=filters,mart=db,values=ID)
}


#############################
####RNA-Seq data analysis####
#############################
####loading count data####
dat.use_lot_B52_53_H49_58=read.table(file = "/Volumes/cluster/Wei-Yun/single_fly_B27_B28_H4_H9_filtered.csv",sep = ",",header = T,row.names = 1)

sample_ID=colnames(dat.use_lot_B52_53_H49_58);sample_ID=as.data.frame(sample_ID)
sample_ID[,2]=strsplit2(sample_ID[,1],"_")[,2];sample_ID[,3]=strsplit2(sample_ID[,1],"_")[,4]
colnames(sample_ID)=c("sample","evo","pop")
sample_ID[,4]=paste0(sample_ID$evo,"_",sample_ID$pop)
colnames(sample_ID)=c("sample","evo","pop","evo_pop")
table(sample_ID$evo_pop)

####normalization####
y=DGEList(counts=dat.use_lot_B52_53_H49_58,group = sample_ID$evo_pop)
y=calcNormFactors(y)

count.use_old=dat.use_lot_B52_53_H49_58[,1:43]
lib.size_old=colSums(count.use_old)
evo=strsplit2(colnames(count.use_old),"_")[,2]
y1=DGEList(counts=count.use_old,group = evo)
y1=calcNormFactors(y1)
count.use_new=dat.use_lot_B52_53_H49_58[,44:81]
lib.size_new=colSums(count.use_new)
evo=strsplit2(colnames(count.use_new),"_")[,2]
y2=DGEList(counts=count.use_new,group = evo)
y2=calcNormFactors(y2)

dat_B_H_all = log(cpm(y))
dat_B_H_all_old=log(cpm(y1))
dat_B_H_all_new=log(cpm(y2))
dat_B_WY_old = dat_B_H_all_old[,which(substr(colnames(dat_B_H_all_old),3,3)=="B")]
dat_H_WY_old = dat_B_H_all_old[,which(substr(colnames(dat_B_H_all_old),3,3)=="H")]
dat_B_WY_new = dat_B_H_all_new[,which(substr(colnames(dat_B_H_all_new),3,3)=="B")]
dat_H_WY_new = dat_B_H_all_new[,which(substr(colnames(dat_B_H_all_new),3,3)=="H")]



####F-test####
#old
ftest_result_old=matrix(data = NA,nrow = nrow(dat_B_H_all),ncol = 4)
ftest_result_old[,1]=row.names(dat_B_H_all)
colnames(ftest_result_old)=c("gene_name","F","p.val","p.adj")
for (i in 1:nrow(ftest_result_old)) {
  ftest_result_old[i,2]=as.numeric(var.test(x = dat_H_WY_old[i,],y = dat_B_WY_old[i,])$statistic)
  ftest_result_old[i,3]=as.numeric(var.test(x = dat_H_WY_old[i,],y = dat_B_WY_old[i,])$p.value)
}
ftest_result_old[,4]=as.numeric(p.adjust(ftest_result_old[,3],method = "BH"))
ftest_result_old_sig=ftest_result_old[which(as.numeric(paste(ftest_result_old[,4]))<0.05),]
ftest_result_old_sig=as.data.frame(ftest_result_old_sig)
#new
ftest_result_new=matrix(data = NA,nrow = nrow(dat_B_H_all),ncol = 4)
ftest_result_new[,1]=row.names(dat_B_H_all)
colnames(ftest_result_new)=c("gene_name","F","p.val","p.adj")
for (i in 1:nrow(ftest_result_new)) {
  ftest_result_new[i,2]=as.numeric(var.test(x = dat_H_WY_new[i,],y = dat_B_WY_new[i,])$statistic)
  ftest_result_new[i,3]=as.numeric(var.test(x = dat_H_WY_new[i,],y = dat_B_WY_new[i,])$p.value)
}
ftest_result_new[,4]=as.numeric(p.adjust(ftest_result_new[,3],method = "BH"))
ftest_result_new_sig=ftest_result_new[which(as.numeric(paste(ftest_result_new[,4]))<0.05),]
ftest_result_new_sig=as.data.frame(ftest_result_new_sig)

#ds
JI=c()
JI2=c()
set.seed(100)
j=1
while(j<=100){
  print(j)
  idx_H1=sample(1:21,15,replace = F)
  idx_B1=sample(1:22,15,replace = F)
  idx_H2=sample(1:21,15,replace = F)
  idx_B2=sample(1:22,15,replace = F)
  ftest_result_old_1=data.frame(matrix(data = NA,nrow = nrow(dat_B_H_all),ncol = 4))
  ftest_result_old_1[,1]=row.names(dat_B_H_all)
  colnames(ftest_result_old_1)=c("gene_name","F","p.val","p.adj")
  rownames(ftest_result_old_1)=ftest_result_old_1$gene_name
  for (i in 1:nrow(ftest_result_old_1)) {
    ftest_result_old_1[i,c(2,3)]=as.numeric(var.test(x = dat_H_WY_old[i,idx_H1],y = dat_B_WY_old[i,idx_B1])[c(1,3)])
  }
  ftest_result_old_1[,4]=as.numeric(p.adjust(ftest_result_old_1[,3],method = "BH"))
  
  ftest_result_old_2=data.frame(matrix(data = NA,nrow = nrow(dat_B_H_all),ncol = 4))
  ftest_result_old_2[,1]=row.names(dat_B_H_all)
  colnames(ftest_result_old_2)=c("gene_name","F","p.val","p.adj")
  rownames(ftest_result_old_2)=ftest_result_old_2$gene_name
  for (i in 1:nrow(ftest_result_old_2)) {
    ftest_result_old_2[i,c(2,3)]=as.numeric(var.test(x = dat_H_WY_old[i,idx_H2],y = dat_B_WY_old[i,idx_B2])[c(1,3)])
  }
  ftest_result_old_2[,4]=as.numeric(p.adjust(ftest_result_old_2[,3],method = "BH"))
  
  inter=sum((ftest_result_old_1$F<1&ftest_result_old_1$p.adj<0.05)&(ftest_result_old_2$F<1&ftest_result_old_2$p.adj<0.05))
  uni=sum((ftest_result_old_1$F<1&ftest_result_old_1$p.adj<0.05)|(ftest_result_old_2$F<1&ftest_result_old_2$p.adj<0.05))
  JI=c(JI,inter/uni)
  inter2=sum((ftest_result_old_1$F<1&ftest_result_old_1$p.adj<0.1)&(ftest_result_old_2$F<1&ftest_result_old_2$p.adj<0.1))
  uni2=sum((ftest_result_old_1$F<1&ftest_result_old_1$p.adj<0.1)|(ftest_result_old_2$F<1&ftest_result_old_2$p.adj<0.1))
  JI2=c(JI2,inter2/uni2)
  
  j=j+1
}

JI_e=c()
JI_e2=c()
set.seed(100)
j=1
while(j<=100){
  print(j)
  idx_H1=sample(1:21,15,replace = F)
  idx_B1=sample(1:22,15,replace = F)
  idx_H2=sample(1:19,15,replace = F)
  idx_B2=sample(1:19,15,replace = F)
  ftest_result_old_1=data.frame(matrix(data = NA,nrow = nrow(dat_B_H_all),ncol = 4))
  ftest_result_old_1[,1]=row.names(dat_B_H_all)
  colnames(ftest_result_old_1)=c("gene_name","F","p.val","p.adj")
  rownames(ftest_result_old_1)=ftest_result_old_1$gene_name
  for (i in 1:nrow(ftest_result_old_1)) {
    ftest_result_old_1[i,c(2,3)]=as.numeric(var.test(x = dat_H_WY_old[i,idx_H1],y = dat_B_WY_old[i,idx_B1])[c(1,3)])
  }
  ftest_result_old_1[,4]=as.numeric(p.adjust(ftest_result_old_1[,3],method = "BH"))
  
  ftest_result_old_2=data.frame(matrix(data = NA,nrow = nrow(dat_B_H_all),ncol = 4))
  ftest_result_old_2[,1]=row.names(dat_B_H_all)
  colnames(ftest_result_old_2)=c("gene_name","F","p.val","p.adj")
  rownames(ftest_result_old_2)=ftest_result_old_2$gene_name
  for (i in 1:nrow(ftest_result_old_2)) {
    ftest_result_old_2[i,c(2,3)]=as.numeric(var.test(x = dat_H_WY_new[i,idx_H2],y = dat_B_WY_new[i,idx_B2])[c(1,3)])
  }
  ftest_result_old_2[,4]=as.numeric(p.adjust(ftest_result_old_2[,3],method = "BH"))
  
  inter=sum((ftest_result_old_1$F<1&ftest_result_old_1$p.adj<0.05)&(ftest_result_old_2$F<1&ftest_result_old_2$p.adj<0.05))
  uni=sum((ftest_result_old_1$F<1&ftest_result_old_1$p.adj<0.05)|(ftest_result_old_2$F<1&ftest_result_old_2$p.adj<0.05))
  JI_e=c(JI_e,inter/uni)
  inter2=sum((ftest_result_old_1$F<1&ftest_result_old_1$p.adj<0.1)&(ftest_result_old_2$F<1&ftest_result_old_2$p.adj<0.1))
  uni2=sum((ftest_result_old_1$F<1&ftest_result_old_1$p.adj<0.1)|(ftest_result_old_2$F<1&ftest_result_old_2$p.adj<0.1))
  JI_e2=c(JI_e2,inter2/uni2)
  
  j=j+1
}




JI.emp=length(intersect(VD1.FDR.1,VD2.FDR.1))/length(union(VD1.FDR.1,VD2.FDR.1))

sum(JI<JI.emp)
#png("Fig.SX.png")

boxplot(JI_e,JI,col=c("salmon","grey"),names=c("between-population",
                                               "within-population"),
        ylab="Consistency")

boxplot(JI_e2,JI2,col=c("salmon","grey"),names=c("between-population",
                                                 "within-population"),
        ylab="Consistency")

####DE-test-limma####
evo_old=strsplit2(colnames(count.use_old),"_")[,2]
ModelDesign_old=model.matrix(~0+evo_old)
GLM_old=lmFit(dat_B_H_all_old,design = ModelDesign_old)
my_contrast_old=makeContrasts(evo_oldH-evo_oldB,levels = colnames(coef(GLM_old)))
myfit_old=contrasts.fit(GLM_old,my_contrast_old)
eBayes_old=eBayes(myfit_old)
res_old=topTable(eBayes_old, sort.by = "none", n = Inf)
sum(res_old$adj.P.Val<0.05)
res_old=cbind(row.names(res_old),res_old)
colnames(res_old)[1]="gene_name"
res_old_sig=subset(res_old,res_old$adj.P.Val<0.05)
write.csv(res_old,file = "/Volumes/cluster/Wei-Yun/review_ME/res_old.csv",sep = ",",quote = F,
          row.names = F)

#new
evo_new=strsplit2(colnames(count.use_new),"_")[,2]
ModelDesign_new=model.matrix(~0+evo_new)
GLM_new=lmFit(dat_B_H_all_new,design = ModelDesign_new)
my_contrast_new=makeContrasts(evo_newH-evo_newB,levels = colnames(coef(GLM_new)))
myfit_new=contrasts.fit(GLM_new,my_contrast_new)
eBayes_new=eBayes(myfit_new)
res_new=topTable(eBayes_new, sort.by = "none", n = Inf)
sum(res_new$adj.P.Val<0.05)
res_new=cbind(row.names(res_new),res_new)
colnames(res_new)[1]="gene_name"
res_new_sig=subset(res_new,res_new$adj.P.Val<0.05)
write.csv(res_new,file = "/Volumes/cluster/Wei-Yun/review_ME/res_new.csv",sep = ",",quote = F,
          row.names = F)

#comparison to Ana's candidates
candid=read.table("Ana_Expression_changes.txt",sep = "\t",header = T,stringsAsFactors = F)
candid_ind=candid$Gene[which(candid$FDR_BH<0.05)]
length(intersect(candid_ind,res_new_sig$gene_name))
fisher.test(cont_table(candid_ind,res_old$gene_name,res_old_sig$gene_name),alternative = "greater")
fisher.test(cont_table(candid_ind,res_new$gene_name,res_new_sig$gene_name),alternative = "greater")
fisher.test(cont_table(candid_ind,res_old$gene_name,res_old_edgeR_sig$gene_name),alternative = "greater")
fisher.test(cont_table(candid_ind,res_new$gene_name,res_new_edgeR_sig$gene_name),alternative = "greater")

####PCA analysis####
pca=prcomp(t(log(cpm(y))))#pca analysis
group_col=c()
for (i in 1:length(sample_ID$evo)) {group_col[i]=gsub("H","red",sample_ID$evo[i])}
for (i in 1:length(group_col)) {group_col[i]=gsub("B","blue",group_col[i])}
group_col_all=c(group_col,rep("blue",5),rep("red",6))
group_col
group_pch=c()
for (i in 1:length(sample_ID$pop)) {group_pch[i]=gsub("27",16,sample_ID$pop[i])}
for (i in 1:length(group_pch)) {group_pch[i]=gsub("4",16,group_pch[i])}
for (i in 1:length(group_pch)) {group_pch[i]=gsub("28",17,group_pch[i])}
for (i in 1:length(group_pch)) {group_pch[i]=gsub("9",17,group_pch[i])}
group_pch=as.numeric(group_pch)
group_col_all=c(group_col,rep("blue",5),rep("red",6))
group_pch=c(group_pch,rep(1,5),rep(1,3),rep(2,3))

#png("/Volumes/cluster/Wei-Yun/figure/plot_v1/plot_PCA_all.png",width = 8.7,height = 8.7,units = "cm",res=600,pointsize = 6)
par(mfrow=c(1,1))
par(mar=c(4.1,4.1,3,2))
plot(pca$x[,c(1,2)],pch=group_pch,col=group_col,xlab="PC1(11.9%)",ylab="PC2(9.97%)",main="PCA plot",cex.lab=1.2)
legend("topleft",legend = c("Anc.27","Anc.28","Evo.4","Evo.9"),col = c("blue","blue","red","red"),pch = c(16,17,16,17),bty = "n",horiz = T,cex = 1.2)
#dev.off()


####evolution of variance####
png("/Volumes/cluster/Wei-Yun/figure/plot_v1/F.png",width = 8.7,height = 10,units = "cm",res=600,pointsize = 6)
layout(matrix(c(1,2,3,4),2,2,byrow = T),widths = c(1,1),heights = c(0.5,1,2))
layout.show(4)
par(mar=c(0,4,4,0.5))
a=boxplot(log(as.numeric(paste(ftest_result_old[which(as.numeric(paste(ftest_result_old[,2]))<1&as.numeric(paste(ftest_result_old[,4]))<0.05),2]))),
          log(as.numeric(paste(ftest_result_old[which(as.numeric(paste(ftest_result_old[,4]))>0.05),2]))),
          log(as.numeric(paste(ftest_result_old[which(as.numeric(paste(ftest_result_old[,2]))>1&as.numeric(paste(ftest_result_old[,4]))<0.05),2]))),
          col = c("#fc8d62","gainsboro","#66c2a5"),main="Replicate1",horizontal=T,frame=F,axes=F)
text(a$stats[nrow(a$stats),]+1.32,c(1:3),c("n=125","n=10417","n=41"))
par(mar=c(0,0.5,4,4))
b=boxplot(log(as.numeric(paste(ftest_result_new[which(as.numeric(paste(ftest_result_new[,2]))<1&as.numeric(paste(ftest_result_new[,4]))<0.05),2]))),
          log(as.numeric(paste(ftest_result_new[which(as.numeric(paste(ftest_result_new[,4]))>0.05),2]))),
          log(as.numeric(paste(ftest_result_new[which(as.numeric(paste(ftest_result_new[,2]))>1&as.numeric(paste(ftest_result_new[,4]))<0.05),2]))),
          col = c("#fc8d62","gainsboro","#66c2a5"),main="Replicate2",horizontal=T,frame=F,axes=F)
text(b$stats[nrow(b$stats),]+1.32,c(1:3),c("n=97","n=10435","n=51"))
par(mar=c(5,4,0,0.5))
hist(log(as.numeric(paste(ftest_result_old[,2]))),main="",frame=F,ylab="",yaxt="n",xlab="ln(F)",cex.lab=1.5,breaks = 50,xlim = c(-4,6))
hist(log(as.numeric(paste(ftest_result_old[,2])))[which(as.numeric(paste(ftest_result_old[,2]))<1&as.numeric(paste(ftest_result_old[,4]))<0.05)],col="#fc8d62",add=T)
hist(log(as.numeric(paste(ftest_result_old[,2])))[which(as.numeric(paste(ftest_result_old[,2]))>1&as.numeric(paste(ftest_result_old[,4]))<0.05)],col="#66c2a5",add=T,breaks = 20)
abline(v = 0, col="red",lty=2,lwd=2)
par(mar=c(5,0.5,0,4))
hist(log(as.numeric(paste(ftest_result_new[,2]))),main="",frame=F,yaxt='n',ylab="",xlab="ln(F)",cex.lab=1.5,breaks=50,xlim = c(-4,6))
hist(log(as.numeric(paste(ftest_result_new[,2])))[which(as.numeric(paste(ftest_result_new[,2]))<1&as.numeric(paste(ftest_result_new[,4]))<0.05)],col="#fc8d62",add=T,breaks = 10)
hist(log(as.numeric(paste(ftest_result_new[,2])))[which(as.numeric(paste(ftest_result_new[,2]))>1&as.numeric(paste(ftest_result_new[,4]))<0.05)],col="#66c2a5",add=T,breaks = 20)
abline(v=0,col="red",lty=2,lwd=2)
par(mar=c(5,5,3,4))
diff_var_4=log(as.numeric(paste(ftest_result_old[,2])))
diff_var_9=log(as.numeric(paste(ftest_result_new[,2])))
plot(diff_var_4,diff_var_9,
     xlab='lnF (Replicate1)',
     ylab="lnF (Replicate2)",
     pch=19,asp = 1,col=alpha("grey",0.2),cex.lab=1.5,type="n")
abline(h=0,col="grey",lty=4)
abline(v=0,col="grey",lty=4)
points(diff_var_4[which(as.numeric(ftest_result_new[,4])>=0.05&as.numeric(ftest_result_old[,4])>=0.05)],
       diff_var_9[which(as.numeric(ftest_result_new[,4])>=0.05&as.numeric(ftest_result_old[,4])>=0.05)],
       col=alpha("gainsboro",1),pch=19)
points(diff_var_4[which(as.numeric(ftest_result_old[,4])<=0.05&as.numeric(ftest_result_old[,2])<=1)],
       diff_var_9[which(as.numeric(ftest_result_old[,4])<=0.05&as.numeric(ftest_result_old[,2])<=1)],
       col=alpha("#fc8d62",1),pch=19)
points(diff_var_4[which(as.numeric(ftest_result_new[,4])<=0.05&as.numeric(ftest_result_new[,2])<=1)],
       diff_var_9[which(as.numeric(ftest_result_new[,4])<=0.05&as.numeric(ftest_result_new[,2])<=1)],
       col=alpha("#fc8d62",1),pch=19)
points(diff_var_4[which(as.numeric(ftest_result_old[,4])<=0.05&as.numeric(ftest_result_old[,2])>1)],
       diff_var_9[which(as.numeric(ftest_result_old[,4])<=0.05&as.numeric(ftest_result_old[,2])>1)],
       col=alpha("#83c266",1),pch=19)
points(diff_var_4[which(as.numeric(ftest_result_new[,4])<=0.05&as.numeric(ftest_result_new[,2])>1)],
       diff_var_9[which(as.numeric(ftest_result_new[,4])<=0.05&as.numeric(ftest_result_new[,2])>1)],
       col=alpha("#83c266",1),pch=19)
points(diff_var_4[which(as.numeric(ftest_result_new[,4])<=0.05&as.numeric(ftest_result_old[,4])<=0.05&as.numeric(ftest_result_new[,2])<1&as.numeric(ftest_result_old[,2])<1)],
       diff_var_9[which(as.numeric(ftest_result_new[,4])<=0.05&as.numeric(ftest_result_old[,4])<=0.05&as.numeric(ftest_result_new[,2])<1&as.numeric(ftest_result_old[,2])<1)],
       col=alpha("#a83832",1),pch=19)
points(diff_var_4[which(as.numeric(ftest_result_new[,4])<=0.05&as.numeric(ftest_result_old[,4])<=0.05&as.numeric(ftest_result_new[,2])>1&as.numeric(ftest_result_old[,2])>1)],
       diff_var_9[which(as.numeric(ftest_result_new[,4])<=0.05&as.numeric(ftest_result_old[,4])<=0.05&as.numeric(ftest_result_new[,2])>1&as.numeric(ftest_result_old[,2])>1)],
       col=alpha("#286351",1),pch=19)
legend("topleft",legend = c("decrease in one replicate","decrease in both replicates","increase in one replicate","increase in both replicates","n.s."),
       col = c(alpha("#fc8d62",1),"#a83832",alpha("#83c266",1),"#286351",alpha("gainsboro",1)),
       pch = 19,bty = "n")
dev.off()

####tissue enrichment analysis####
#fig.4
#flyatlas2
#old
atlas_enrichment=read.table(file = "/Volumes/cluster/Wei-Yun/20190121/flyaltlas2_log2fc.txt",sep = " ",header = T)
flyatlas2=atlas_enrichment[,1:15]
Dat_VAR_mt_S_old=as.character(ftest_result_old_sig$gene_name[which(as.numeric(paste(ftest_result_old_sig$F))<1)])
background=row.names(res_old)
tet_result_VARS=matrix(NA,nrow = 13,ncol = 5)
for (i in 2:14) {
  tissue_specific_ID1=flyatlas2$FBgnID[!is.na(flyatlas2[,i])&flyatlas2[,i]>2]
  tet_result_VARS[(i-1),1]=strsplit2(colnames(flyatlas2)[i],"_")[3]
  tet_result_VARS[(i-1),2]=length(tissue_specific_ID1)
  tet_result_VARS[(i-1),3]=fisher.test(cont_table(Dat_VAR_mt_S_old,background,tissue_specific_ID1),alternative = "greater")$estimate
  tet_result_VARS[(i-1),4]=fisher.test(cont_table(Dat_VAR_mt_S_old,background,tissue_specific_ID1),alternative = "greater")$p.value
}
tet_result_VARS[,5]=p.adjust(as.numeric(tet_result_VARS[,4]))

Dat_VAR_mt_L_old=as.character(ftest_result_old_sig$gene_name[which(as.numeric(paste(ftest_result_old_sig$F))>1)])
tet_result_VARL=matrix(NA,nrow = 13,ncol = 5)
for (i in 2:14) {
  tissue_specific_ID1=flyatlas2$FBgnID[!is.na(flyatlas2[,i])&flyatlas2[,i]>2]
  tet_result_VARL[(i-1),1]=strsplit2(colnames(flyatlas2)[i],"_")[3]
  tet_result_VARL[(i-1),2]=length(tissue_specific_ID1)
  tet_result_VARL[(i-1),3]=fisher.test(cont_table(Dat_VAR_mt_L_old,background,tissue_specific_ID1),alternative = "greater")$estimate
  tet_result_VARL[(i-1),4]=fisher.test(cont_table(Dat_VAR_mt_L_old,background,tissue_specific_ID1),alternative = "greater")$p.value
}
tet_result_VARL[,5]=p.adjust(as.numeric(tet_result_VARL[,4]))
tet.ori=cbind(as.numeric(paste(tet_result_VARS[,5])),as.numeric(paste(tet_result_VARL[,5])))
tet.use=abs(log(tet.ori))

#new
atlas_enrichment=read.table(file = "/Volumes/cluster/Wei-Yun/20190121/flyaltlas2_log2fc.txt",sep = " ",header = T)
flyatlas2=atlas_enrichment[,1:15]
Dat_VAR_mt_S_new=as.character(ftest_result_new_sig$gene_name[which(as.numeric(paste(ftest_result_new_sig$F))<1)])
background=row.names(res_new)
tet_result_VARS=matrix(NA,nrow = 13,ncol = 5)
for (i in 2:14) {
  tissue_specific_ID1=flyatlas2$FBgnID[!is.na(flyatlas2[,i])&flyatlas2[,i]>2]
  tet_result_VARS[(i-1),1]=strsplit2(colnames(flyatlas2)[i],"_")[3]
  tet_result_VARS[(i-1),2]=length(tissue_specific_ID1)
  tet_result_VARS[(i-1),3]=fisher.test(cont_table(Dat_VAR_mt_S_new,background,tissue_specific_ID1),alternative = "greater")$estimate
  tet_result_VARS[(i-1),4]=fisher.test(cont_table(Dat_VAR_mt_S_new,background,tissue_specific_ID1),alternative = "greater")$p.value
}
tet_result_VARS[,5]=p.adjust(as.numeric(tet_result_VARS[,4]))

Dat_VAR_mt_L_new=as.character(ftest_result_new_sig$gene_name[which(as.numeric(paste(ftest_result_new_sig$F))>1)])
tet_result_VARL=matrix(NA,nrow = 13,ncol = 5)
for (i in 2:14) {
  tissue_specific_ID1=flyatlas2$FBgnID[!is.na(flyatlas2[,i])&flyatlas2[,i]>2]
  tet_result_VARL[(i-1),1]=strsplit2(colnames(flyatlas2)[i],"_")[3]
  tet_result_VARL[(i-1),2]=length(tissue_specific_ID1)
  tet_result_VARL[(i-1),3]=fisher.test(cont_table(Dat_VAR_mt_L_new,background,tissue_specific_ID1),alternative = "greater")$estimate
  tet_result_VARL[(i-1),4]=fisher.test(cont_table(Dat_VAR_mt_L_new,background,tissue_specific_ID1),alternative = "greater")$p.value
}
tet_result_VARL[,5]=p.adjust(as.numeric(tet_result_VARL[,4]))
tet.ori=cbind(as.numeric(paste(tet_result_VARS[,5])),as.numeric(paste(tet_result_VARL[,5])))
tet.use=cbind(tet.use,abs(log(tet.ori)))
colnames(tet.use)=c("Rep1.Dec","Rep1.Inc","Rep2.Dec","Rep2.Inc")
row.names(tet.use)=tet_result_VARL[,1]

####enrichment test for VC and non-VC####
D1=which(res_old$gene_name%in%ftest_result_old[ftest_result_old[,2]<1,1])
D2=which(res_new$gene_name%in%ftest_result_new[ftest_result_new[,2]<1,1])
I1=which(res_old$gene_name%in%ftest_result_old[ftest_result_old[,2]>1,1])
I2=which(res_new$gene_name%in%ftest_result_new[ftest_result_new[,2]>1,1])
DE1=which(res_old$gene_name%in%res_old_sig$gene_name)
DE2=which(res_new$gene_name%in%res_new_sig$gene_name)
VC1=which(res_old$gene_name%in%ftest_result_old_sig$gene_name)
VC2=which(res_new$gene_name%in%ftest_result_new_sig$gene_name)
VD1=which(res_old$gene_name%in%ftest_result_old_sig$gene_name[as.numeric(paste(ftest_result_old_sig[,2]))<1])
VD2=which(res_new$gene_name%in%ftest_result_new_sig$gene_name[as.numeric(paste(ftest_result_new_sig[,2]))<1])
VI1=which(res_old$gene_name%in%ftest_result_old_sig$gene_name[as.numeric(paste(ftest_result_old_sig[,2]))>1])
VI2=which(res_new$gene_name%in%ftest_result_new_sig$gene_name[as.numeric(paste(ftest_result_new_sig[,2]))>1])
VDp=intersect(VD1,VD2)
VIp=intersect(VI1,VI2)
VCp=c(VDp,VIp)
US1=union(DE1,VC1)
US2=union(DE2,VC2)


c1=cont_table(DE1,1:10583,VD1)
prop.test(t(c1))
c2=cont_table(DE2,1:10583,VD2)
prop.test(t(c2))
p.adjust(c(0.046,0.164),method = 'BH')

#correlation
par(mfrow=c(2,2))
plot(abs(log(as.numeric(paste(ftest_result_old[D1,2])))),abs(res_old$logFC[D1]),
     pch=19,xlim=c(0,3.5),ylim=c(0,3.5),col="gainsboro",xlab="",ylab="")
points(abs(log(as.numeric(paste(ftest_result_old[intersect(D1,DE1),2])))),
       abs(res_old$logFC)[intersect(D1,DE1)],pch=19,col="grey35")
points(abs(log(as.numeric(paste(ftest_result_old[intersect(D1,VD1),2])))),
       abs(res_old$logFC)[intersect(D1,VD1)],pch=19,col="#66c2a5")
points(abs(log(as.numeric(paste(ftest_result_old[Reduce(intersect,list(D1,VD1,DE1)),2])))),
       abs(res_old$logFC)[Reduce(intersect,list(D1,VD1,DE1))],pch=19,col="darkblue")

plot(abs(log(as.numeric(paste(ftest_result_new[D2,2])))),abs(res_new$logFC[D2]),
     pch=19,xlim=c(0,3.5),ylim=c(0,3.5),col="gainsboro",xlab="",ylab="")
points(abs(log(as.numeric(paste(ftest_result_new[intersect(D2,DE2),2])))),
       abs(res_new$logFC)[intersect(D2,DE2)],pch=19,col="grey35")
points(abs(log(as.numeric(paste(ftest_result_new[intersect(D2,VD2),2])))),
       abs(res_new$logFC)[intersect(D2,VD2)],pch=19,col="#66c2a5")
points(abs(log(as.numeric(paste(ftest_result_new[Reduce(intersect,list(D2,VD2,DE2)),2])))),
       abs(res_new$logFC)[Reduce(intersect,list(D2,VD2,DE2))],pch=19,col="darkblue")

plot(abs(log(as.numeric(paste(ftest_result_old[I1,2])))),abs(res_old$logFC[I1]),
     pch=19,xlim=c(0,3.5),ylim=c(0,3.5),col="gainsboro",xlab="",ylab="")
points(abs(log(as.numeric(paste(ftest_result_old[intersect(I1,DE1),2])))),
       abs(res_old$logFC)[intersect(I1,DE1)],pch=19,col="grey35")
points(abs(log(as.numeric(paste(ftest_result_old[intersect(I1,VI1),2])))),
       abs(res_old$logFC)[intersect(I1,VI1)],pch=19,col="#fc8d62")
points(abs(log(as.numeric(paste(ftest_result_old[Reduce(intersect,list(I1,VI1,DE1)),2])))),
       abs(res_old$logFC)[Reduce(intersect,list(I1,VI1,DE1))],pch=19,col="firebrick")

plot(abs(log(as.numeric(paste(ftest_result_new[I2,2])))),abs(res_new$logFC[I2]),
     pch=19,xlim=c(0,3.5),ylim=c(0,3.5),col="gainsboro",xlab="",ylab="")
points(abs(log(as.numeric(paste(ftest_result_new[intersect(I2,DE2),2])))),
       abs(res_new$logFC)[intersect(I2,DE2)],pch=19,col="grey35")
points(abs(log(as.numeric(paste(ftest_result_new[intersect(I2,VI2),2])))),
       abs(res_new$logFC)[intersect(I2,VI2)],pch=19,col="#fc8d62")
points(abs(log(as.numeric(paste(ftest_result_new[Reduce(intersect,list(I2,VI2,DE2)),2])))),
       abs(res_new$logFC)[Reduce(intersect,list(I2,VI2,DE2))],pch=19,col="firebrick")

cor.test(abs(log(as.numeric(paste(ftest_result_old[D1,2])))),abs(res_old$logFC[D1]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_new[D2,2])))),abs(res_new$logFC[D2]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_old[I1,2])))),abs(res_old$logFC[I1]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_new[I2,2])))),abs(res_new$logFC[I2]),method = "spearman")

cor.test(abs(log(as.numeric(paste(ftest_result_old[Reduce(union,list(D1,VD1,DE1)),2])))),
         abs(res_old$logFC[Reduce(union,list(D1,VD1,DE1))]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_new[Reduce(union,list(D2,VD2,DE2)),2])))),
         abs(res_new$logFC[Reduce(union,list(D2,VD2,DE2))]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_old[Reduce(union,list(I1,VI1,DE1)),2])))),
         abs(res_old$logFC[Reduce(union,list(I1,VI1,DE1))]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_new[Reduce(union,list(I2,VI2,DE2)),2])))),
         abs(res_new$logFC[Reduce(union,list(I2,VI2,DE2))]),method = "spearman")


plot(log(as.numeric(paste(ftest_result_old[DE1,2]))),abs(res_old$logFC[DE1]),pch=19)
plot(abs(log(as.numeric(paste(ftest_result_old[DE1,2])))),abs(res_old$logFC[DE1]),pch=19)
cor.test(log(as.numeric(paste(ftest_result_old[DE1,2]))),abs(res_old$logFC[DE1]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_old[DE1,2])))),abs(res_old$logFC[DE1]),method = "spearman")

plot(log(as.numeric(paste(ftest_result_new[DE2,2]))),abs(res_new$logFC[DE2]),pch=19)
plot(abs(log(as.numeric(paste(ftest_result_new[DE2,2])))),abs(res_new$logFC[DE2]),pch=19)
cor.test(log(as.numeric(paste(ftest_result_new[DE2,2]))),abs(res_new$logFC[DE2]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_new[DE2,2])))),abs(res_new$logFC[DE2]),method = "spearman")

plot(log(as.numeric(paste(ftest_result_old[VC1,2]))),abs(res_old$logFC[VC1]),pch=19)
plot(abs(log(as.numeric(paste(ftest_result_old[VC1,2])))),abs(res_old$logFC[VC1]),pch=19)
cor.test(log(as.numeric(paste(ftest_result_old[VC1,2]))),abs(res_old$logFC[VC1]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_old[VC1,2])))),abs(res_old$logFC[VC1]),method = "spearman")

plot(log(as.numeric(paste(ftest_result_new[VC2,2]))),abs(res_new$logFC[VC2]),pch=19)
plot(abs(log(as.numeric(paste(ftest_result_new[VC2,2])))),abs(res_new$logFC[VC2]),pch=19)
cor.test(log(as.numeric(paste(ftest_result_new[VC2,2]))),abs(res_new$logFC[VC2]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_new[VC2,2])))),abs(res_new$logFC[VC2]),method = "spearman")

plot(abs(log(as.numeric(paste(ftest_result_old[VD1,2])))),abs(res_old$logFC[VD1]),pch=19)
cor.test(abs(log(as.numeric(paste(ftest_result_old[VD1,2])))),abs(res_old$logFC[VD1]),method = "spearman")
plot(abs(log(as.numeric(paste(ftest_result_new[VD2,2])))),abs(res_new$logFC[VD2]),pch=19)
cor.test(abs(log(as.numeric(paste(ftest_result_new[VD2,2])))),abs(res_new$logFC[VD2]),method = "spearman")

plot(abs(log(as.numeric(paste(ftest_result_old[VI1,2])))),abs(res_old$logFC[VI1]),pch=19)
cor.test(abs(log(as.numeric(paste(ftest_result_old[VI1,2])))),abs(res_old$logFC[VI1]),method = "spearman")
plot(abs(log(as.numeric(paste(ftest_result_new[VI2,2])))),abs(res_new$logFC[VI2]),pch=19)
cor.test(abs(log(as.numeric(paste(ftest_result_new[VI2,2])))),abs(res_new$logFC[VI2]),method = "spearman")

plot(abs(log(as.numeric(paste(ftest_result_old[VCp,2])))),abs(res_old$logFC[VCp]),pch=19)
cor.test(abs(log(as.numeric(paste(ftest_result_old[VCp,2])))),abs(res_old$logFC[VCp]),method = "spearman")
plot(abs(log(as.numeric(paste(ftest_result_new[VCp,2])))),abs(res_new$logFC[VCp]),pch=19)
cor.test(abs(log(as.numeric(paste(ftest_result_new[VCp,2])))),abs(res_new$logFC[VCp]),method = "spearman")

plot(abs(log(as.numeric(paste(ftest_result_old[VDp,2])))),abs(res_old$logFC[VDp]),pch=19)
cor.test(abs(log(as.numeric(paste(ftest_result_old[VDp,2])))),abs(res_old$logFC[VDp]),method = "spearman")
plot(abs(log(as.numeric(paste(ftest_result_new[VDp,2])))),abs(res_new$logFC[VDp]),pch=19)
cor.test(abs(log(as.numeric(paste(ftest_result_new[VDp,2])))),abs(res_new$logFC[VDp]),method = "spearman")

plot(log(as.numeric(paste(ftest_result_old[US1,2]))),abs(res_old$logFC[US1]),pch=19)
plot(abs(log(as.numeric(paste(ftest_result_old[US1,2])))),abs(res_old$logFC[US1]),pch=19)
cor.test(log(as.numeric(paste(ftest_result_old[US1,2]))),abs(res_old$logFC[US1]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_old[US1,2])))),abs(res_old$logFC[US1]),method = "spearman")

plot(log(as.numeric(paste(ftest_result_new[US2,2]))),abs(res_new$logFC[US2]),pch=19)
plot(abs(log(as.numeric(paste(ftest_result_new[US2,2])))),abs(res_new$logFC[US2]),pch=19)
cor.test(log(as.numeric(paste(ftest_result_new[US2,2]))),abs(res_new$logFC[US2]),method = "spearman")
cor.test(abs(log(as.numeric(paste(ftest_result_new[US2,2])))),abs(res_new$logFC[US2]),method = "spearman")


boxplot(abs(res_old[which(as.numeric(paste(ftest_result_old[,2]))<1&as.numeric(paste(ftest_result_old[,4]))<0.05),2]),
        abs(res_old[-which(as.numeric(paste(ftest_result_old[,2]))<1&as.numeric(paste(ftest_result_old[,4]))<0.05),2]))

t.test(abs(res_old[which(as.numeric(paste(ftest_result_old[,2]))<1&as.numeric(paste(ftest_result_old[,4]))<0.05),2]),
       abs(res_old[-which(as.numeric(paste(ftest_result_old[,2]))<1&as.numeric(paste(ftest_result_old[,4]))<0.05),2]))


####TF enrichment analysis####
####RcisTarget####
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
# To support paralell execution:
BiocManager::install(c("doMC", "doRNG"))
# For the examples in the follow-up section of the tutorial:
BiocManager::install(c("DT", "visNetwork"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RcisTarget")
library(RcisTarget)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library(biomaRt)

install.packages("doParallel")
library(doParallel)

ID_VD=list(old=as.character(ftest_result_old_sig[which(as.numeric(paste(ftest_result_old_sig[,2]))<1),1]),
           new=as.character(ftest_result_new_sig[which(as.numeric(paste(ftest_result_new_sig[,2]))<1),1]))

ID_converter=function(ID,db,attributes,filters){
  getBM(attributes=attributes,filters=filters,mart=db,values=ID)
}

ensembl=useMart("ensembl",dataset="dmelanogaster_gene_ensembl",host = "http://useast.ensembl.org/" )
background_conv=ftest_result_old[,1]
latest_FB_ID=read.delim("/Volumes/cluster/Wei-Yun/TFBS/Rcistarget_DB/fbgn_annotation_ID_fb_2018_01_R.txt",header = T,stringsAsFactors = F,sep="\t")
for(i in 1:dim(latest_FB_ID)[1]){
  background_conv[which(background_conv%in%strsplit2(latest_FB_ID[i,4],","))]=latest_FB_ID[i,3]
}
query_ID_new=lapply(ID_VD,function(x) {
  for(i in 1:dim(latest_FB_ID)[1]){
    x[which(x%in%strsplit2(latest_FB_ID[i,4],","))]=latest_FB_ID[i,3]
  }
  return(x)
})
conv_query=lapply(query_ID_new,function(x) ID_converter(ID = x,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1])
conv_background=ID_converter(ID=background_conv,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1]
download.file("https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc8nr/gene_based/dm6-5kb-upstream-full-tx-11species.mc8nr.feather",
              destfile = "/Volumes/cluster/Wei-Yun/TFBS/Rcistarget_DB/dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
motifranking=importRankings("/Volumes/cluster/Wei-Yun/TFBS/Rcistarget_DB/dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
motifEnrichmentTable_VD=lapply(conv_query[1:2],function(x) {
  cisTarget(as.character(x),motifRankings = motifranking,motifAnnot = motifAnnotations_dmel_v8,nesThreshold = 3,nCores = 8,
            geneErnMethod = "iCisTarget",aucMaxRank = 0.01*ncol(motifranking))
})
all_TFs=intersect(unique(motifAnnotations_dmel_v8$TF),conv_background)
TFs_VD=lapply(motifEnrichmentTable_VD,function(x) {
  tmp=c(x$TF_highConf)
  genes <- gsub(" \\(.*\\). ", "; ", tmp, fixed=FALSE)
  genesSplit <- unique(unlist(strsplit(genes, "; ")))
  genesSplit <- unique(strsplit2(genesSplit, " ")[,1])
  return(genesSplit)
})

gut_tissue_gene=gut_tissue_gene=intersect(unique(as.character(atlas_enrichment_foldchange[which(atlas_enrichment_foldchange[,7]>1),1])),rownames(cpm(y))[apply(cpm(y),1,function(x) mean(x)>1)])
Dat_VAR_mt_S_old=ftest_result_old_sig$gene_name[which(as.numeric(paste(ftest_result_old_sig$F))<1)]
Dat_VAR_mt_S_old=as.data.frame(Dat_VAR_mt_S_old)
gut_gene_old=intersect(gut_tissue_gene,Dat_VAR_mt_S_old$Dat_VAR_mt_S_old)
Dat_VAR_mt_S_new=ftest_result_new_sig$gene_name[which(as.numeric(paste(ftest_result_new_sig$F))<1)]
Dat_VAR_mt_S_new=as.data.frame(Dat_VAR_mt_S_new)
gut_gene_new=intersect(gut_tissue_gene,Dat_VAR_mt_S_new$Dat_VAR_mt_S_new)

ID_VD=list(old=gut_gene_old,new=gut_gene_new)

ID_converter=function(ID,db,attributes,filters){
  getBM(attributes=attributes,filters=filters,mart=db,values=ID)
}

ensembl=useMart("ensembl",dataset="dmelanogaster_gene_ensembl",host = "http://useast.ensembl.org/" )
background=ftest_result_old[,1]
background_conv=background
latest_FB_ID=read.delim("/Volumes/cluster/Wei-Yun/TFBS/Rcistarget_DB/fbgn_annotation_ID_fb_2018_01_R.txt",header = T,stringsAsFactors = F,sep="\t")
for(i in 1:dim(latest_FB_ID)[1]){
  background_conv[which(background_conv%in%strsplit2(latest_FB_ID[i,4],","))]=latest_FB_ID[i,3]
}
query_ID_new=lapply(ID_VD,function(x) {
  for(i in 1:dim(latest_FB_ID)[1]){
    x[which(x%in%strsplit2(latest_FB_ID[i,4],","))]=latest_FB_ID[i,3]
  }
  return(x)
})
conv_query=lapply(query_ID_new,function(x) ID_converter(ID = x,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1])
conv_background=ID_converter(ID=background_conv,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1]
download.file("https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc8nr/gene_based/dm6-5kb-upstream-full-tx-11species.mc8nr.feather",
              destfile = "/Volumes/cluster/Wei-Yun/TFBS/Rcistarget_DB/dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
motifranking=importRankings("/Volumes/cluster/Wei-Yun/TFBS/Rcistarget_DB/dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
motifEnrichmentTable_VD=lapply(conv_query[1:2],function(x) {
  cisTarget(as.character(x),motifRankings = motifranking,motifAnnot = motifAnnotations_dmel_v8,nesThreshold = 5,nCores = 8,
            geneErnMethod = "iCisTarget",aucMaxRank = 0.01*ncol(motifranking))
})
all_TFs=intersect(unique(motifAnnotations_dmel_v8$TF),conv_background)
TFs_VD=lapply(motifEnrichmentTable_VD,function(x) {
  tmp=c(x$TF_highConf)
  genes <- gsub(" \\(.*\\). ", "; ", tmp, fixed=FALSE)
  genesSplit <- unique(unlist(strsplit(genes, "; ")))
  genesSplit <- unique(strsplit2(genesSplit, " ")[,1])
  return(genesSplit)
})
lapply(TFs_VD,function(x) sum(x%in%conv_query$new))

ID_DE=list(old=as.character(res_old_sig$gene_name),new=as.character(res_new_sig$gene_name))
DE_ID_new=lapply(ID_DE,function(x) {
  for(i in 1:dim(latest_FB_ID)[1]){
    x[which(x%in%strsplit2(latest_FB_ID[i,4],","))]=latest_FB_ID[i,3]
  }
  return(x)
})
conv_query=lapply(DE_ID_new,function(x) ID_converter(ID = x,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1])
TF_conv=lapply(TFs_VD,function(x) ID_converter(ID = x,db = ensembl,attributes = c("external_gene_name","flybase_gene_id"),filters = "external_gene_name"))
TF_conv
gtf_transform=read.delim(file = "/Volumes/cluster/Wei-Yun/reference_Dsim/Genome_annotation/gtf_resolved_table.txt",header = F,stringsAsFactors = F)
gene_chr=tapply(gtf_transform[,1],INDEX = gtf_transform$V5,unique)
gene_start=tapply(gtf_transform[,2],INDEX = gtf_transform$V5,min)
gene_end=tapply(gtf_transform[,3],INDEX = gtf_transform$V5,max)
gene_strand=tapply(gtf_transform[,4],INDEX = gtf_transform$V5,unique)
gene_summary=data.frame(gene_chr,gene_start,gene_end,gene_strand,sort(unique(gtf_transform$V5)))
colnames(gene_summary)=c("Chr","start","end","strand","FB_ID")
gene_summary_filtered=gene_summary[background,]
gene_summary_filtered$FB_ID_new=c(background_conv)
gene_summary_filtered=gene_summary_filtered[gene_summary_filtered[,1]%in%c("2L","2R","3L","3R","4","X"),]
rownames(gene_summary_filtered)=gene_summary_filtered$FB_ID

lapply(TF_conv,function(x) gene_summary_filtered[x[,2],])


####Enrichment for digestive genes####
setwd("/Users/weiyun/Dropbox (PopGen)/application_for_popgen_vienna/updated/")
query_ID_1=ftest_result_old_sig[as.numeric(paste0(ftest_result_old_sig[,2]))<1,1]
query_ID_2=ftest_result_new_sig[as.numeric(paste0(ftest_result_new_sig[,2]))<1,1]
background=ftest_result_old[,1]
DS_genes=read.table("./digestive_system.txt",stringsAsFactors = F,header = F)[,1]
DG_genes=read.table("./digestive_enzyme.txt",stringsAsFactors = F,header = F)[,1]
fisher.test(cont_table(query_ID_1,background,DG_genes),alternative = "greater")
fisher.test(cont_table(query_ID_2,background,DG_genes),alternative = "greater")


####overlay with genomic selection signature####
gtf_transform=read.table("/Volumes/cluster/Wei-Yun/reference_Dsim/Genome_annotation/gtf_resolved_table.txt",stringsAsFactors = F,header = F)
gene_chr=tapply(gtf_transform[,1],INDEX = gtf_transform$V5,unique)
gene_start=tapply(gtf_transform[,2],INDEX = gtf_transform$V5,min)
gene_end=tapply(gtf_transform[,3],INDEX = gtf_transform$V5,max)
gene_strand=tapply(gtf_transform[,4],INDEX = gtf_transform$V5,unique)
gene_summary=data.frame(gene_chr,gene_start,gene_end,gene_strand,sort(unique(gtf_transform$V5)))
colnames(gene_summary)=c("Chr","start","end","strand","FB_ID")
gene_summary_filtered=gene_summary[background,]

gene_summary_filtered_tmp=gene_summary_filtered[!duplicated(gene_summary_filtered$FB_ID),]
rownames(gene_summary_filtered_tmp)=gene_summary_filtered_tmp$FB_ID
# all_TFs_out=cbind(all_TFs_tmp[,1:2],gene_summary_filtered_tmp[all_TFs_tmp$flybase_gene_id,1:4],
#                   sex_bias=apply(sapply(sex_bias_ID_new[c(1,2)],function(x) ifelse(all_TFs_tmp$flybase_gene_id%in%x,1,0)),1,function(x) ifelse(sum(x)==0,NA,names(which(x==1)))),
#                   evolution=apply(sapply(query_ID_new[-c(5,8,9,12)],function(x) ifelse(all_TFs_tmp$flybase_gene_id%in%x,1,0)),1,function(x) ifelse(sum(x)==0,NA,names(which(x==1)))),
#                   description=all_TFs_tmp[,3])

AF_filter=readRDS("/Volumes/cluster/Wei-Yun/genomic_data/fl_Dsim_hot_cage_AF.rds")
#SNP_blockID=read.delim("/Volumes/Temp1/shengkai/common_info/Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync",header=F,stringsAsFactors = F)[,c(1:3,85,86)]
#FET_res=read.delim("/Volumes/Temp1/shengkai/common_info/Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync",header=F,stringsAsFactors = F)[,c(1:3,75:84)]
#saveRDS(FET_res,"/Volumes/Temp1/shengkai/common_info/fl_Dsim_hot_cage_FET.rds")
#saveRDS(SNP_blockID,"/Volumes/Temp1/shengkai/common_info/fl_Dsim_hot_cage_blockID.rds")
SNP_blockID=readRDS("/Volumes/cluster/Wei-Yun/genomic_data/fl_Dsim_hot_cage_blockID.rds")
FET_res=readRDS("/Volumes/cluster/Wei-Yun/genomic_data/fl_Dsim_hot_cage_FET.rds")

CMH_res=cbind(AF_filter[,c(1:3,24)],SNP_blockID[,-1:-3])
colnames(CMH_res)=c("CHR","BP","Allele","CMH","blockID_0.75","blockID_0.35")
CMH_res=pseudo_chr(CMH_res)
CMH_res=pseudo_pos(CMH_res)
CMH_res=SNPID_gen(CMH_res)
CMH_res$P=10^-CMH_res$CMH

colnames(FET_res)=c("CHR","BP","Allele",paste0("H",sprintf("%02d",1:10)))
FET_res=pseudo_chr(FET_res)
FET_res=pseudo_pos(FET_res)
FET_res=SNPID_gen(FET_res)

VD1_tab=data.frame(flybase_gene_id=as.character(query_ID_var1),gene_summary_filtered_tmp[as.character(query_ID_var1),1:4])
VD1_query=data.frame(gene_ID=VD1_tab$flybase_gene_id,Chr=VD1_tab$Chr,start=VD1_tab$start-5000,end=VD1_tab$end+5000)

SNP_VD1=c()
for (i in 1:dim(VD1_query)[1]){
  print(i)
  SNP_VD1=rbind(SNP_VD1,snp_identifier(VD1_query[i,],CMH_res,1e-32,1e-19))
}
sum(!SNP_VD1[,5]==0)/length(SNP_VD1[,5])

VD2_tab=data.frame(flybase_gene_id=as.character(query_ID_var2),gene_summary_filtered_tmp[as.character(query_ID_var2),1:4])
VD2_query=data.frame(gene_ID=VD2_tab$flybase_gene_id,Chr=VD2_tab$Chr,start=VD2_tab$start-5000,end=VD2_tab$end+5000)

SNP_VD2=c()
for (i in 1:dim(VD2_query)[1]){
  print(i)
  SNP_VD2=rbind(SNP_VD2,snp_identifier(VD2_query[i,],CMH_res,1e-32,1e-19))
}
sum(!SNP_VD2[,5]==0)/length(SNP_VD2[,5])

TFVD1=c("FBgn0003277","FBgn0010355","FBgn0003687","FBgn0020912","FBgn0032223","FBgn0038391","FBgn0003507","FBgn0003896","FBgn0037698","FBgn0002914","FBgn0037540","FBgn0014018","FBgn0051388","FBgn0040078","FBgn0002023","FBgn0028996","FBgn0003964","FBgn0011648")
TFVD1_tab=data.frame(flybase_gene_id=as.character(TFVD1),gene_summary_filtered_tmp[as.character(TFVD1),1:4])
TFVD1_query=data.frame(gene_ID=TFVD1_tab$flybase_gene_id,Chr=TFVD1_tab$Chr,start=TFVD1_tab$start-5000,end=TFVD1_tab$end+5000)

SNP_TFVD1=c()
for (i in 1:dim(TFVD1_query)[1]){
  print(i)
  SNP_TFVD1=rbind(SNP_TFVD1,snp_identifier(TFVD1_query[i,],CMH_res,1e-32,1e-19))
}
sum(!SNP_TFVD1[,5]==0)/length(SNP_TFVD1[,5])

TFVD2=c("FBgn0005612","FBgn0000577","FBgn0039938","FBgn0024288","FBgn0263112","FBgn0001269","FBgn0005613","FBgn0029711")
TFVD2_tab=data.frame(flybase_gene_id=as.character(TFVD2),gene_summary_filtered_tmp[as.character(TFVD2),1:4])
TFVD2_query=data.frame(gene_ID=TFVD2_tab$flybase_gene_id,Chr=TFVD2_tab$Chr,start=TFVD2_tab$start-5000,end=TFVD2_tab$end+5000)

SNP_TFVD2=c()
for (i in 1:dim(TFVD2_query)[1]){
  print(i)
  SNP_TFVD2=rbind(SNP_TFVD2,snp_identifier(TFVD2_query[i,],CMH_res,1e-32,1e-19))
}
sum(!SNP_TFVD2[,5]==0)/length(SNP_TFVD2[,5])


# tmp=apply(SNP_VD1,1,function(x) CMH_res[CMH_res$SNPID%in%unlist(strsplit(x[6],",")),])
# names(tmp)=SNP_VD1$gene_ID

r.set=list()
for (i in 1:100){
  r.set[[i]]=as.character(gene_summary_filtered$FB_ID[sample(1:dim(gene_summary_filtered)[1],100,replace = F)])
}
r.query=lapply(r.set,function(x) data.frame(gene_ID=x,Chr=gene_summary_filtered[x,]$Chr,start=gene_summary_filtered[x,]$start-5000,end=gene_summary_filtered[x,]$end+5000,stringsAsFactors = F))

random_res=lapply(r.query,function(x) {
  tmp.res=c()
  for (i in 1:dim(x)[1]){
    print(i)
    tmp.res=rbind(tmp.res,snp_identifier(x[i,],CMH_res,1e-32,1e-19))
  }
  return(tmp.res)
})
quantile(sapply(random_res,function(x) sum(!x[,5]==0)/length(x[,5])),probs=c(0.05,0.5,0.9,0.95,0.99))

hist(sapply(random_res,function(x) sum(!x[,5]==0)/length(x[,5])),breaks = 10,main = "",
     xlab = "Prop. of overlap")
abline(v=sum(!SNP_VD1[,5]==0)/length(SNP_VD1[,5]),col="red")
abline(v=sum(!SNP_VD2[,5]==0)/length(SNP_VD2[,5]),col="firebrick")


abline(v=sum(!SNP_TFVD1[,5]==0)/length(SNP_TFVD1[,5]),col="blue")

abline(v=sum(!SNP_TFVD2[,5]==0)/length(SNP_TFVD2[,5]),col="darkblue")


# SNP_TOI_all=c()
# for (i in 1:dim(TOI_query)[1]){
#   print(i)
#   tmp=cbind(CMH_res[CMH_res$CHR%in%TOI_query$Chr[i]&CMH_res$BP>=TOI_query$start[i]&CMH_res$BP<=TOI_query$end[i],],gene_ID=TOI_query$gene_ID[i])
#   SNP_TOI_all=rbind(SNP_TOI_all,tmp)
# }
# 
# table(SNP_TOI_all$gene_ID)
# 
# median_AF=apply(AF_filter[rownames(SNP_TOI_all),4:23],2,function(x) tapply(x,SNP_TOI_all$gene_ID,function(y) median(y,na.rm=T)))
# num_gene=c()
# for (i in 1:10){
#   tmp=median_AF[,-1:-10]-median_AF[,1:10]
#   num_gene=c(num_gene,sum(tmp[,i]>apply(test,2,function(x) quantile(x,0.99))[i]))
# }
# num_repl=c()
# for (i in 1:19){
#   tmp=median_AF[,-1:-10]-median_AF[,1:10]
#   num_repl=c(num_repl,sum(tmp[i,]>apply(test,2,function(x) quantile(x,0.99))))
# }
# 
# test=c()
# i=1
# while (i<=100){
#   random_index=sample(1:dim(AF_filter)[1],1000,replace = F)
#   test=rbind(test,apply(as.matrix(AF_filter[random_index,14:23]),2,median)-apply(as.matrix(AF_filter[random_index,4:13]),2,median))
#   i=i+1
# }
# apply(test,2,function(x) quantile(x,0.99))
# 
# num_sig_SNP=c()
# i=1
# while (i<=100){
#   random_SNP=CMH_res[sample(1:dim(CMH_res)[1],21454,replace = F),]
#   num_sig_SNP=c(num_sig_SNP,sum(random_SNP$CHR%in%"X"&random_SNP$P<1e-32|!random_SNP$CHR%in%"X"&random_SNP$P<1e-19))
#   i=i+1
# }
# quantile(num_sig_SNP,0.01)
# 
# e_p_FET_SNP=c()
# for(g in 1:19){
#   FET_SNP_TOI=FET_res[FET_res$SNPID%in%SNP_TOI_all[SNP_TOI_all$gene_ID%in%unique(SNP_TOI_all$gene_ID)[g],]$SNPID,]
#   #apply(FET_SNP_TOI[,4:13],2,function(x) sum(x>19,na.rm = T))
#   num_sig_SNP_FET=c()
#   for (i in 1:10){
#     num_sig_SNP_FET=c(num_sig_SNP_FET,sum(!FET_SNP_TOI$CHR%in%"X"&FET_SNP_TOI[,i+3]>c(28,32,27,34,38,35,38,27,44,36)[i]|FET_SNP_TOI$CHR%in%"X"&FET_SNP_TOI[,i+3]>c(42,Inf,33,Inf,36,Inf,36,28,Inf,37)[i],na.rm = T))
#   }
#   
#   num_sig_SNP_repl=c()
#   j=1
#   while (j<=1000){
#     random_seed=sample(1:dim(FET_res)[1],1,replace = F)
#     random_SNP=FET_res[random_seed:(random_seed+dim(FET_SNP_TOI)[1]),]
#     num_sig_SNP_FET_random=c()
#     for (i in 1:10){
#       num_sig_SNP_FET_random=c(num_sig_SNP_FET_random,sum(!random_SNP$CHR%in%"X"&random_SNP[,i+3]>c(28,32,27,34,38,35,20,27,44,36)[i]|random_SNP$CHR%in%"X"&random_SNP[,i+3]>c(42,Inf,33,Inf,36,Inf,36,28,Inf,37)[i],na.rm = T))
#     }
#     num_sig_SNP_repl=rbind(num_sig_SNP_repl,num_sig_SNP_FET_random)
#     j=j+1
#   }
#   num_sig_SNP_FET>apply(num_sig_SNP_repl,2,function(x) quantile(x,0.99))
#   tmp_e_p_FET_SNP=c()
#   for (i in 1:10){
#     tmp_e_p_FET_SNP=c(tmp_e_p_FET_SNP,sum(num_sig_SNP_repl[,i]>num_sig_SNP_FET[i],na.rm = T)/dim(num_sig_SNP_repl)[1])
#   }
#   e_p_FET_SNP=rbind(e_p_FET_SNP,tmp_e_p_FET_SNP)
# }
# colnames(e_p_FET_SNP)=paste0("H",sprintf("%02d",1:10))
# rownames(e_p_FET_SNP)=unique(SNP_TOI_all$gene_ID)
# apply(e_p_FET_SNP,2,function(x) sum(x<0.05))
# apply(e_p_FET_SNP,1,function(x) sum(x<0.05))
# 
# expressed_TFs_query=lapply(expressed_TFs_out[c(5,8,9,12)],function(x) data.frame(gene_ID=x$flybase_gene_id,Chr=x$Chr,start=x$start-5000,end=x$end+5000))
# expressed_TFs_s_query=lapply(expressed_TFs_s_out[c(5,8,9,12)],function(x) data.frame(gene_ID=x$flybase_gene_id,Chr=x$Chr,start=x$start-5000,end=x$end+5000))
# expressed_TFs_res=lapply(expressed_TFs_query,function(x) {
#   tmp.res=c()
#   for (i in 1:dim(x)[1]){
#     print(i)
#     tmp.res=rbind(tmp.res,snp_identifier(x[i,],CMH_res,1e-32,1e-19))
#   }
#   return(tmp.res)
# })
# 
# sapply(expressed_TFs_res,function(x) sum(!x[,5]==0)/length(x[,5]))
# lapply(expressed_TFs_res,function(x) sum(x[,5]>0))
# lapply(expressed_TFs_res,function(x) ID_converter(x[x[,5]>0,1],db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1])
# lapply(lapply(expressed_TFs_res,function(x) ID_converter(x[x[,5]>0,1],db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1]),function(x) intersect(x,unlist(TFs_sex_biased)))
# 
# qn=c()
# for (j in 1:4){
#   r.set=list()
#   for (i in 1:1000){
#     n=sapply(expressed_TFs_res,dim)[1,j]
#     r.set[[i]]=as.character(gene_summary_filtered$FB_ID[sample(1:dim(gene_summary_filtered),n,replace = F)])
#   }
#   r.query=lapply(r.set,function(x) data.frame(gene_ID=x,Chr=gene_summary_filtered[x,]$Chr,start=gene_summary_filtered[x,]$start-5000,end=gene_summary_filtered[x,]$end+5000,stringsAsFactors = F))
#   
#   random_res=lapply(r.query,function(x) {
#     tmp.res=c()
#     for (i in 1:dim(x)[1]){
#       print(i)
#       tmp.res=rbind(tmp.res,snp_identifier(x[i,],CMH_res,1e-32,1e-19))
#     }
#     return(tmp.res)
#   })
#   qn=rbind(qn,quantile(sapply(random_res,function(x) sum(!x[,5]==0)/length(x[,5])),probs=c(0.05,0.5,0.9,0.95,0.99)))
# }
# 
# 
# hist(sapply(random_res,function(x) sum(!x[,5]==0)/length(x[,5])))
# abline(v=sapply(expressed_TFs_res,function(x) sum(!x[,5]==0)/length(x[,5])))
# 
# hl_expressed_TFs_res=unique(unlist(sapply(expressed_TFs_res,function(x) strsplit2(x[!x[,5]==0,6],split = ","))))
# png("./florida_result_fullset_v7/genomic_signatures/manhattan_plot_expressed_TFs_res.png",width = 24,height = 8,units = "cm",pointsize = 8,res = 600)
# manhattan(CMH_res,snp = "SNPID",chr = "pseudo_CHR",chrlabs = c("X","2L","2R","3L","3R","4"),highlight = hl_expressed_TFs_res,
#           suggestiveline = F,genomewideline = 28,cex=0.75,logp = T)
# dev.off()
# 
# hl_SNP=sapply(expressed_TFs_res,function(x) {
#   tmp=strsplit2(x[!x[,5]==0,6],split = ",")
#   tmp_ID=ID_converter(x[!x[,5]==0,1],db = ensembl,attributes = c("external_gene_name","flybase_gene_id"),filters = "flybase_gene_id")
#   rownames(tmp_ID)=tmp_ID[,2]
#   rownames(tmp)=tmp_ID[as.character(x[!x[,5]==0,1]),1]
#   return(tmp)
# })
# png("./florida_result_fullset_v7/genomic_signatures/manhattan_plot.png",width = 24,height = 8,units = "cm",pointsize = 8,res = 600)
# plot(CMH_res$pseudo_POS,CMH_res$CMH,col=ifelse(CMH_res$CHR%in%c("X","2R","3R"),"grey10","grey60"),pch=19,axes=F,cex=0.75
#      ,ylab=expression(-log[10](italic(p))),xlab="Chromosome",ylim=c(0,max(CMH_res$CMH)*1.6))
# axis(1,at=tapply(CMH_res$pseudo_POS,CMH_res$pseudo_CHR,quantile,probs=0.5),labels=c("X","2L","2R","3L","3R","4"))
# axis(2)
# abline(h=c(19,32),col=c("blue","red"),lty=2)
# dev.off()
# 
# png("./florida_result_fullset_v7/genomic_signatures/manhattan_plot_labeled.png",width = 24,height = 8,units = "cm",pointsize = 8,res = 600)
# plot(CMH_res$pseudo_POS,CMH_res$CMH,col=ifelse(CMH_res$CHR%in%c("X","2R","3R"),"grey10","grey60"),pch=19,axes=F,cex=0.75
#      ,ylab=expression(-log[10](italic(p))),xlab="Chromosome",ylim=c(0,max(CMH_res$CMH)*1.6))
# axis(1,at=tapply(CMH_res$pseudo_POS,CMH_res$pseudo_CHR,quantile,probs=0.5),labels=c("X","2L","2R","3L","3R","4"))
# axis(2)
# for (i in 1:length(hl_SNP)){
#   with(CMH_res[CMH_res$SNPID%in%unlist(hl_SNP[[i]]),],points(pseudo_POS,CMH,col=alpha(c("orange","salmon","skyblue","royalblue")[i],0.6),pch=19,cex=0.75))
#   for (j in 1:dim(hl_SNP[[i]])[1]){
#     set.seed(50)
#     with(CMH_res[CMH_res$SNPID%in%hl_SNP[[i]][j,1],],
#          text(pseudo_POS,80+which(unique(unlist(lapply(hl_SNP,rownames)))==rownames(hl_SNP[[i]])[j])*100/length(unique(unlist(lapply(hl_SNP,rownames)))),
#               labels = rownames(hl_SNP[[i]])[j],col=alpha(c("orange","salmon","skyblue","royalblue")[i],0.6),
#               cex=ifelse(rownames(hl_SNP[[i]])[j]%in%unlist(TFs_sex_biased),1,0.5),
#               font = ifelse(rownames(hl_SNP[[i]])[j]%in%unlist(TFs_sex_biased),4,3)))
#   }
# }
# legend("topleft",legend = rev(c("F.down (14/27)","F.up (38/60)","M.down (25/38)","M.up (34/62)")),pch=19,col = rev(alpha(c("orange","salmon","skyblue","royalblue"),0.6)))
# #abline(h=c(19,32),col=c("blue","red"),lty=2)
# dev.off()
# #[sample(1:length(unique(unlist(lapply(hl_SNP,rownames)))),length(unique(unlist(lapply(hl_SNP,rownames)))),replace = F)]
# 
# png("/Volumes/Temp1/shengkai/remap_run/florida_result_fullset_v7/genomic_signatures/Figure_S4.png",height = 25,width = 20,units = "cm",pointsize = 10,res = 600)
# par(mfrow=c(5,2))
# for(i in which(TOI_query$gene_ID%in%SNP_TOI$gene_ID[SNP_TOI$SNP_counts>0])){
#   plot(NA,xlim=c(0,TOI_query[i,4]-TOI_query[i,3]),ylim=c(-5,55),xlab=TOI_query[i,2],ylab=expression(-log[10](p)),main=TOI_query[i,1],axes=F)
#   axis(1,at=round(summary(0:(TOI_query[i,4]-TOI_query[i,3]))),labels = round(summary(TOI_query[i,3]:TOI_query[i,4])))
#   axis(2)
#   rect(xleft = gtf_transform[gtf_transform$V5%in%TOI_query[i,1],2]-TOI_query[i,3],
#        xright = gtf_transform[gtf_transform$V5%in%TOI_query[i,1],3]-TOI_query[i,3],
#        ybottom = -6,ytop = -2,col = "black")
#   Arrows(x0 = 0,x1 = TOI_query[i,4]-TOI_query[i,3],y0 = -4,y1 = -4,
#          code = ifelse(gene_summary_filtered[as.character(TOI_query[i,1]),]$strand=="+",2,1),lwd=3,arr.type = "triangle",arr.length = 0.1,arr.width = 0.125)
#   points(SNP_TOI_all[SNP_TOI_all$gene_ID%in%TOI_query[i,1],2]-TOI_query[i,3],SNP_TOI_all[SNP_TOI_all$gene_ID%in%TOI_query[i,1],4],
#          pch=19,col=ifelse(SNP_TOI_all[SNP_TOI_all$gene_ID%in%TOI_query[i,1],4]>19,"green","grey"),cex=0.75)
#   abline(h=19,col="red",lty=2)
# }
# dev.off()
# 123



####figure2.update####
png(filename = "/Volumes/cluster/Wei-Yun/review_ME/figure2_v1.png",width = 8.7*1.5,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(1,2))
VD_rep1=ftest_result_old_sig$gene_name[which(ftest_result_old_sig$F<1)]
non_VD_rep1=setdiff(row.names(count.use_old),VD_rep1)
plot.use_rep1=matrix(data=c(length(intersect(VD_rep1,res_old_sig$gene_name))/length(VD_rep1),
                            length(setdiff(VD_rep1,res_old_sig$gene_name))/length(VD_rep1),
                            length(intersect(non_VD_rep1,res_old_sig$gene_name))/length(non_VD_rep1),
                            length(setdiff(non_VD_rep1,res_old_sig$gene_name))/length(non_VD_rep1)),
                     2,2,byrow = F)
barplot(plot.use_rep1, main="population 1",
        ylab="Proportion of DE/non-DE genes", col=c("grey30","grey"),
        names = c("variance decrease\n(n=125)","background\n(n=10458)"))
mean_VD_1=length(intersect(VD_rep1,res_old_sig$gene_name))/length(VD_rep1)
se_VD_1=sqrt(mean_VD_1*(1-mean_VD_1)/97)
arrows(0.7, mean_VD_1 - se_VD_1 *2, 0.7,
       mean_VD_1 + se_VD_1 *2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
mean_nVD_1=length(intersect(non_VD_rep1,res_old_sig$gene_name))/length(non_VD_rep1)
se_nVD_1=sqrt(mean_nVD_1*(1-mean_nVD_1)/10486)
arrows(1.9, mean_nVD_1 - se_nVD_1 *2, 1.9,
       mean_nVD_1 + se_nVD_1 *2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)



VD_rep2=ftest_result_new_sig$gene_name[which(ftest_result_new_sig$F<1)]
non_VD_rep2=setdiff(row.names(count.use_old),VD_rep2)
plot.use_rep2=matrix(data=c(length(intersect(VD_rep2,res_new_sig$gene_name))/length(VD_rep2),
                            length(setdiff(VD_rep2,res_new_sig$gene_name))/length(VD_rep2),
                            length(intersect(non_VD_rep2,res_new_sig$gene_name))/length(non_VD_rep2),
                            length(setdiff(non_VD_rep2,res_new_sig$gene_name))/length(non_VD_rep2)),
                     2,2,byrow = F)
barplot(plot.use_rep2, main="population 2",
        ylab="Proportion of DE/non-DE genes", col=c("grey30","grey"),
        names = c("variance decrease\n(n=97)","background\n(n=10486)"))
mean_VD_2=length(intersect(VD_rep2,res_new_sig$gene_name))/length(VD_rep2)
se_VD_2=sqrt(mean_VD_2*(1-mean_VD_2)/97)
arrows(0.7, mean_VD_2 - se_VD_2 *2, 0.7,
       mean_VD_2 + se_VD_2 *2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
mean_nVD_2=length(intersect(non_VD_rep2,res_new_sig$gene_name))/length(non_VD_rep2)
se_nVD_2=sqrt(mean_nVD_2*(1-mean_nVD_2)/10486)
arrows(1.9, mean_nVD_2 - se_nVD_2 *2, 1.9,
       mean_nVD_2 + se_nVD_2 *2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off()






####supplementary information |logFC| and logF
png(filename = "/Volumes/cluster/Wei-Yun/review_ME/SF2_v1.png",width = 15.4,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(1,2))
plot(log(as.numeric(paste(ftest_result_old[,2]))),abs(res_old$logFC),pch=19,col="gainsboro",
     xlab="variance change",ylab="absolute mean change",xlim=c(-4,4),ylim=c(0,2),main="population 1")
points(log(as.numeric(paste(ftest_result_old[,2])))[res_old$adj.P.Val<0.05],
       abs(res_old$logFC)[res_old$adj.P.Val<0.05],pch=19,col="grey35")
points(log(as.numeric(paste(ftest_result_old[,2])))[as.numeric(paste(ftest_result_old[,4]))<0.05&as.numeric(paste(ftest_result_old[,2]))>1],
       abs(res_old$logFC)[as.numeric(paste(ftest_result_old[,4]))<0.05&as.numeric(paste(ftest_result_old[,2]))>1],pch=19,col="#66c2a5")
points(log(as.numeric(paste(ftest_result_old[,2])))[as.numeric(paste(ftest_result_old[,4]))<0.05&as.numeric(paste(ftest_result_old[,2]))<1],
       abs(res_old$logFC)[as.numeric(paste(ftest_result_old[,4]))<0.05&as.numeric(paste(ftest_result_old[,2]))<1],pch=19,col="#fc8d62")
points(log(as.numeric(paste(ftest_result_old[,2])))[as.numeric(paste(ftest_result_old[,4]))<0.05&as.numeric(paste(ftest_result_old[,2]))>1&res_old$adj.P.Val<0.05],
       abs(res_old$logFC)[as.numeric(paste(ftest_result_old[,4]))<0.05&as.numeric(paste(ftest_result_old[,2]))>1&res_old$adj.P.Val<0.05],pch=19,col="darkblue")
points(log(as.numeric(paste(ftest_result_old[,2])))[as.numeric(paste(ftest_result_old[,4]))<0.05&as.numeric(paste(ftest_result_old[,2]))<1&res_old$adj.P.Val<0.05],
       abs(res_old$logFC)[as.numeric(paste(ftest_result_old[,4]))<0.05&as.numeric(paste(ftest_result_old[,2]))<1&res_old$adj.P.Val<0.05],pch=19,col="firebrick")


plot(log(as.numeric(paste(ftest_result_new[,2]))),abs(res_new$logFC),pch=19,col="gainsboro",
     xlab="variance change",ylab="absolute mean change",xlim=c(-4,4),ylim=c(0,2),main="population 2")
points(log(as.numeric(paste(ftest_result_new[,2])))[res_new$adj.P.Val<0.05],
       abs(res_new$logFC)[res_new$adj.P.Val<0.05],pch=19,col="grey35")
points(log(as.numeric(paste(ftest_result_new[,2])))[as.numeric(paste(ftest_result_new[,4]))<0.05&as.numeric(paste(ftest_result_new[,2]))>1],
       abs(res_new$logFC)[as.numeric(paste(ftest_result_new[,4]))<0.05&as.numeric(paste(ftest_result_new[,2]))>1],pch=19,col="#66c2a5")
points(log(as.numeric(paste(ftest_result_new[,2])))[as.numeric(paste(ftest_result_new[,4]))<0.05&as.numeric(paste(ftest_result_new[,2]))<1],
       abs(res_new$logFC)[as.numeric(paste(ftest_result_new[,4]))<0.05&as.numeric(paste(ftest_result_new[,2]))<1],pch=19,col="#fc8d62")
points(log(as.numeric(paste(ftest_result_new[,2])))[as.numeric(paste(ftest_result_new[,4]))<0.05&as.numeric(paste(ftest_result_new[,2]))>1&res_new$adj.P.Val<0.05],
       abs(res_new$logFC)[as.numeric(paste(ftest_result_new[,4]))<0.05&as.numeric(paste(ftest_result_new[,2]))>1&res_new$adj.P.Val<0.05],pch=19,col="darkblue")
points(log(as.numeric(paste(ftest_result_new[,2])))[as.numeric(paste(ftest_result_new[,4]))<0.05&as.numeric(paste(ftest_result_new[,2]))<1&res_new$adj.P.Val<0.05],
       abs(res_new$logFC)[as.numeric(paste(ftest_result_new[,4]))<0.05&as.numeric(paste(ftest_result_new[,2]))<1&res_new$adj.P.Val<0.05],pch=19,col="firebrick")
dev.off()

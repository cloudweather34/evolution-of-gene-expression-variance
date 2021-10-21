png("/Users/weiyun/Dropbox (PopGen)/application_for_popgen_vienna/manuscript/fig2_ms2.png",width = 10,height = 6,units = "cm",res=600,pointsize = 6)
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

dev.off()


png("/Users/weiyun/Dropbox (PopGen)/application_for_popgen_vienna/manuscript/fig3_ms2.png",width = 8.7,height = 8.7,units = "cm",res=600,pointsize = 6)
par(mfrow=c(2,1))
barplot(t(tet.use[,c(1,2)]),main="Replicate 1", xlab = "Tissue",
        ylab="-log(FDR)", col=c("#fc8d62","#66c2a5"),beside=TRUE)
legend("topright",legend = c("Decreased","Increased"),fill = c("#fc8d62","#66c2a5"), bty="n")

barplot(t(tet.use[,c(3,4)]),main="Replicate 2", xlab = "Tissue",
        ylab="-log(FDR)", col=c("#fc8d62","#66c2a5"), beside=TRUE)
dev.off()


setwd("/Volumes/cluster/Wei-Yun/microbiome_seq/pair_end/")
dir.list=list.files(pattern = "brackG")
temp=read.table(file = dir.list[1],header = T)[,c(1,6)]
colnames(temp)=c("name","B1")
for (i in 2:length(dir.list)) {
  aa=read.table(file = dir.list[i],sep = "\t",header = T,stringsAsFactors = F)[,c(1,6)]
  colnames(aa)=c("name",strsplit2(dir.list[i],"_")[1,1])
  temp=merge(temp,aa,by="name",all=T)
}

microb_noWol=temp[-49,-c(5,6)]

for (i in 2:9) {
  microb_noWol[,i]=as.numeric(paste(microb_noWol[,i]))
}

richness=c()
for (i in 2:9) {
  richness[i]=length(which(microb_noWol[,i]>5))
}
microb_noWol[is.na(microb_noWol)]=0

a=sample(5:9,3)
microb_noWol_base=microb_noWol[apply(microb_noWol[,2:4],1,function(x) any(x>5)),]

microb_noWol_hot=microb_noWol[apply(microb_noWol[,a],1,function(x) any(x>5)),]


gamma_diversity=c(56,46)
alpha_diversity=c()
alpha_diversity[1]=mean(richness[2:4])
alpha_diversity[2]=mean(richness[a])
beta_diversity=gamma_diversity/alpha_diversity


data=microb_noWol[,-1]
row.names(data)=microb_noWol[,1]
data[microb_noWol[,2:9]<10]=0
colnames(data)=c("B1","B2","B3","H1","H2","H3","H4","H5")

#Transform this data in %
data_percentage=apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
data_percentage=data_percentage[rowSums(data_percentage)!=0,]


#create color palette:
library(RColorBrewer)
n <- 62
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


# Make a stacked barplot--> it will be in %!
png(filename = "/Users/weiyun/Dropbox (PopGen)/application_for_popgen_vienna/manuscript/fig4_ms2.png",height = 12,width = 22,units = "cm",res = 600,pointsize = 9)
set.seed(444)
mycol=sample(col_vector,n)
layout(matrix(1:2,1,2),widths = c(2,1))
layout.show(2)
par(mar=c(5,4,4,2))
barplot(data_percentage, col=mycol , border=NA, xlab="population",ylab = "relative abundance (%)")
par(mar=c(0,0,5,0))
plot(NA,type="n",axes = F,xlim = c(0,1),ylim=c(0,1),xlab = "",ylab="",main="genus level")
legend("topleft",legend = rownames(data_percentage)[1:22],fill=mycol[1:22],bty="n")
legend("topright",legend = rownames(data_percentage)[23:44],fill=mycol[23:44],bty="n")
dev.off()


setwd("/Volumes/cluster/Wei-Yun/simulation/") 
ff_0.5=lapply(list.files(path = "formal_0opt_100ef_0.5st_20/",pattern = ".gpf",full.names = T)[-86],function(x) read.table(x))
#ff_3.6=lapply(list.files(path = "test_opt0_100ef_3.6st_20/",pattern = ".gpf",full.names = T)[-6],function(x) read.table(x))
#ff_5.4=lapply(list.files(path = "test_0opt_100ef_5.4st_20/",pattern = ".gpf",full.names = T),function(x) read.table(x))
ff_neutral=lapply(list.files(path = "test_neutral_100ef_20/",pattern = ".gpf",full.names = T),function(x) read.table(x))

list_ff=lapply(list(ff_0.5,ff_neutral),function(x) x[-length(x)])
names(list_ff)=c(0.5,"neutral")
list_loci_no_g0=lapply(list_ff,function(x) lapply(x,function(y) y[!y$V2==0,]))
#phenotypic variance
F_list=lapply(list_ff,function(x) sapply(x,function(y) {
  var_dat=tapply(y$V5,y$V2,var)
  var_dat1=var_dat[-1]
  return(var_dat1/var_dat1[1])
}))

png("~/Dropbox (PopGen)/application_for_popgen_vienna/manuscript/fig5_ms2.png",
    width = 8.7,height = 11.6,units = "cm",res=600,pointsize = 8)
par(mfrow=c(2,1))
par(mar=c(5,4,2,4))
plot(density(pheno),xlim=c(-10,10),ylim=c(0,0.5),xlab=expression("trait value "(sigma[anc])),main="")
lines(density(ff1)$x,density(ff1)$y/max(density(ff1)$y)*max(density(pheno)$y),col="blue",lty=3)
#lines(density(ff2)$x,density(ff2)$y/max(density(ff2)$y)*max(density(pheno)$y),col="blue",lty=3)
axis(4,at=seq(0,max(density(pheno)$y),length.out = 3),labels = c(0.5,1,1.5))
mtext("Fitness",side = 4,line = 2.5)
legend("topleft",lty=c(1,3,3),col=c("black","blue"),bty="n",
       legend = c("Anc. Phenotypic distribution","F.F. (0.5 sd)"))

F_list_mean_ff=sapply(F_list[1:2], function(x) apply(x,1,mean))
F_list_sd_ff=sapply(F_list[1:2], function(x) apply(x,1,sd))

plot(names(F_list_mean_ff[,1]),F_list_mean_ff[,2],type="b",ylim=c(0.6,1.1),
     col="grey",pch=19,xlab="generation",ylab="F-statistics")
points(names(F_list_mean_ff[,1]),F_list_mean_ff[,1],type="b",col="blue",pch=19)


legend("topright",pch=19,col=c("grey","blue"),bty="n",
       legend = c("neutral","F.F. (0.5 sd)"))
dev.off()



png(filename = "~/Dropbox (PopGen)/application_for_popgen_vienna/manuscript/fig4_ms2_ff.png",
    width=12,height = 4.45,units = "cm",pointsize = 9,res=600)
par(mfrow=c(1,3))
a=density(runif(100000,1,2))
plot(a,axes=F,main="",xlab="",ylab="")
polygon(a,col="grey",border = "brown",lty=2,lwd=2)
axis(side=2, at=c(0,0.5,1))
mtext("Fitness", side = 2, line=2.5, at=0.5)
axis(side=1, at=c(1,1.5,2))
b=density(rnorm(100000,1.5,0.125))
plot(b,axes=F,main="Fitness Function",xlab="",ylab="")
polygon(b,col="grey",border = "brown",lty=2,lwd=2)
axis(side=1, at=c(1,1.5,2))
c=density(rnorm(100000,1.5,0.125))
plot(c,axes=F,main="",xlab="",ylab="")
polygon(c,col="grey",border = "brown",lty=2,lwd=2)
axis(side=1, at=c(1,1.5,2))
dev.off()

png(filename = "~/Dropbox (PopGen)/application_for_popgen_vienna/manuscript/fig4_ms2_pv.png",
    width=12,height = 4.45,units = "cm",pointsize = 9,res=600)
par(mfrow=c(1,3))
a=density(runif(100000,1,2))
plot(a,axes=F,main="",xlab="",ylab="")
polygon(a,col="sandybrown",border = "brown",lty=2,lwd=2)
axis(side=2, at=c(0,0.5,1))
mtext("Density", side = 2, line=2.5, at=0.5)
axis(side=1, at=c(1,1.5,2))
b=density(runif(100000,1,2))
plot(b,axes=F,main="Phenotypic distribution",xlab="phenotypic value",ylab="")
polygon(b,col="sandybrown",border = "brown",lty=2,lwd=2)
axis(side=1, at=c(1,1.5,2))
c=density(rnorm(100000,1.5,0.125))
plot(c,axes=F,main="",xlab="",ylab="")
polygon(c,col="sandybrown",border = "brown",lty=2,lwd=2)
axis(side=1, at=c(1,1.5,2))
dev.off()

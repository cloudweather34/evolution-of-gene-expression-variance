####diversity####
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
row.names(microb_noWol)=microb_noWol[,1]

richness=c()
for (i in 2:9) {
  richness[i]=length(which(microb_noWol[,i]>5))
}
microb_noWol[is.na(microb_noWol)]=0


a=sample(5:9,3)
microb_noWol_filtered=microb_noWol[apply(microb_noWol[,c(2,3,4,6,7,8)],1,function(x) any(x>5)),]
microb_noWol_base=microb_noWol[apply(microb_noWol[,2:4],1,function(x) any(x>5)),]
microb_noWol_hot=microb_noWol[apply(microb_noWol[,a],1,function(x) any(x>5)),]



richness=c()
for (i in 2:9) {
  richness[i]=length(which(microb_noWol[,i]>5))
}

set.seed(666)

boot_base=matrix(NA,100,3)
for (i in 1:100) {
  a=sample(2:4,3,replace = T)
  boot_base[i,1]=mean(richness[a])
  boot_base[i,2]=nrow(microb_noWol[apply(microb_noWol[,a],1,function(x) any(x>5)),])
  boot_base[i,3]=boot_base[i,2]/boot_base[i,1]
}
apply(boot_base,2,mean)
apply(boot_base,2,function(x) sd(x)/sqrt(length(x)))*2

boot_hot=matrix(NA,100,3)
for (i in 1:100) {
  a=sample(6:8,3,replace = T)
  boot_hot[i,1]=mean(richness[a])
  boot_hot[i,2]=nrow(microb_noWol[apply(microb_noWol[,a],1,function(x) any(x>5)),])
  boot_hot[i,3]=boot_hot[i,2]/boot_hot[i,1]
}
apply(boot_hot,2,mean)
apply(boot_hot,2,function(x) sd(x)/sqrt(length(x)))*2



####barplot####
data=microb_noWol_filtered[,c(2:4,6:8)]
data[microb_noWol_filtered[,c(2:4,6:8)]<3]=0
colnames(data)=c("B1","B2","B3","H1","H2","H3")

#create color palette:
library(RColorBrewer)
n <- 57
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Transform this data in %
data_percentage=apply(data, 2, function(x){x*100/sum(x,na.rm=T)})

# Make a stacked barplot--> it will be in %!
#png(filename = "/Volumes/cluster/Wei-Yun/figure/plot_v1/relative_abundance.png",height = 12,width = 22,units = "cm",res = 600,pointsize = 9)
set.seed(1234)
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


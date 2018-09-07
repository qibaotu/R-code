
### Sankey plot
library(googleVis)
library(mclust)
setwd('Z:/QUBIC2/Results/Juan/Yan/RST/Sankey')
Ref <-read.table('Yan_cell.csv',header=T,sep=',')
Predict <-read.table('Yan_KL_f0.85_c0.85_k13.csv',header=T,sep=',')  # this is the results from QUBIC2 f0.85k13c0.85

a <-strsplit(as.character(Ref$Cell_type),'_')
Ref$Cells <-sapply(a,'[',1)

b <-strsplit(as.character(Predict$cell),'_')
Predict$Cells <-sapply(a,'[',1)

SC3 <-read.table('Yan_sc3_cluster_result.csv',header=T,sep=',')   # the results provided by Yawei Niu
adjustedRandIndex(SC3$cluster,Ref$Cluster)

### the flow  between Ref and QUBIC2
N1 <-unique(Ref$Cluster)
N2 <-unique(Predict$cluster)

COUNT <-list()
for (i in 1:length(N1)){
	df <-data.frame(From=rep(N1[i],length(N2)),To=rep(-1,length(N2)),weight=rep(-1,length(N2)))
	for (j in 1:length(N2)){	
	cell_target <-Ref[which(Ref$Cluster==N1[i]),]$Cell_type
	cell_pred <-Predict[which(Predict$cluster==N2[j]),]$cell
	matched <-length(intersect(cell_target,cell_pred))
	df$To[j] <-N2[j]
	df$weight[j] <-matched
	}
	COUNT[[i]] <-df
}
RES <-as.data.frame(do.call(rbind,COUNT))
RES$From <-paste('Cluster',RES$From,sep='')
RES$To <-paste('C',RES$To,sep='')
RES <-subset(RES,weight!=0)


### the flow 
COUNT2 <-list()
N3 <-unique(SC3$cluster)
for (i in 1:length(N1)){
	df <-data.frame(From=rep(N1[i],length(N3)),To=rep(-1,length(N3)),weight=rep(-1,length(N3)))
	for (j in 1:length(N3)){	
	cell_target <-Ref[which(Ref$Cluster==N1[i]),]$Cell_type
	cell_sc3 <-SC3[which(SC3$cluste==N3[j]),]$Cell_type
	matched <-length(intersect(cell_target,cell_sc3))
	df$To[j] <-N3[j]
	df$weight[j] <-matched
	}
	COUNT2[[i]] <-df
}

res <-as.data.frame(do.call(rbind,COUNT2))
res <-subset(res,weight!=0)
res$From <-paste('Cluster',res$From,sep='')
res$To <-paste('SC3',res$To,sep='_')
names(res)<-c('To','From','weight')
res <-res[c('From','To','weight')]

RES <-rbind(RES,res)

sk1 <- gvisSankey(RES, from="From", to="To", weight="weight")
=======
### Sankey plot
library(googleVis)
library(mclust)
setwd('Z:/QUBIC2/Results/Juan/Yan/RST/Sankey')
Ref <-read.table('Yan_cell.csv',header=T,sep=',')
Predict <-read.table('Yan_KL_f0.85_c0.85_k13.csv',header=T,sep=',')  # this is the results from QUBIC2 f0.85k13c0.85

a <-strsplit(as.character(Ref$Cell_type),'_')
Ref$Cells <-sapply(a,'[',1)

b <-strsplit(as.character(Predict$cell),'_')
Predict$Cells <-sapply(a,'[',1)

SC3 <-read.table('Yan_sc3_cluster_result.csv',header=T,sep=',')   # the results provided by Yawei Niu
adjustedRandIndex(SC3$cluster,Ref$Cluster)

### the flow  between Ref and QUBIC2
N1 <-unique(Ref$Cluster)
N2 <-unique(Predict$cluster)

COUNT <-list()
for (i in 1:length(N1)){
	df <-data.frame(From=rep(N1[i],length(N2)),To=rep(-1,length(N2)),weight=rep(-1,length(N2)))
	for (j in 1:length(N2)){	
	cell_target <-Ref[which(Ref$Cluster==N1[i]),]$Cell_type
	cell_pred <-Predict[which(Predict$cluster==N2[j]),]$cell
	matched <-length(intersect(cell_target,cell_pred))
	df$To[j] <-N2[j]
	df$weight[j] <-matched
	}
	COUNT[[i]] <-df
}
RES <-as.data.frame(do.call(rbind,COUNT))
RES$From <-paste('Cluster',RES$From,sep='')
RES$To <-paste('C',RES$To,sep='')
RES <-subset(RES,weight!=0)


### the flow 
COUNT2 <-list()
N3 <-unique(SC3$cluster)
for (i in 1:length(N1)){
	df <-data.frame(From=rep(N1[i],length(N3)),To=rep(-1,length(N3)),weight=rep(-1,length(N3)))
	for (j in 1:length(N3)){	
	cell_target <-Ref[which(Ref$Cluster==N1[i]),]$Cell_type
	cell_sc3 <-SC3[which(SC3$cluste==N3[j]),]$Cell_type
	matched <-length(intersect(cell_target,cell_sc3))
	df$To[j] <-N3[j]
	df$weight[j] <-matched
	}
	COUNT2[[i]] <-df
}

res <-as.data.frame(do.call(rbind,COUNT2))
res <-subset(res,weight!=0)
res$From <-paste('Cluster',res$From,sep='')
res$To <-paste('SC3',res$To,sep='_')
names(res)<-c('To','From','weight')
res <-res[c('From','To','weight')]

RES <-rbind(RES,res)

sk1 <- gvisSankey(RES, from="From", to="To", weight="weight")

plot(sk1)

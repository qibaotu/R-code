## combine TFgene list

setwd('C:/Users/Juan.XIe/Desktop/geneid result - human+mouse/result/human')
files <-list.files(getwd(),pattern='GSM*')
MAP <-read.table('summary.csv',header=T,sep=',')

a <-basename(files)
b <-sapply(strsplit(a,'.',1),'[',1)

RES <-list()
for (i in 1:length(files)){
	temp <-read.table(files[i],header=F)
	NAME <-b[i]
	tf <-MAP[which(MAP$Accession %in% NAME),]$ChIP.antibody
	tf2 <-paste(NAME,tf,sep='_')
	Genes <-temp$V1
	TF <-rep(tf2,length(Genes))
	df <-data.frame(TF=TF,Genes=Genes)
	RES[[i]] <-df		
}

DF <-as.data.frame(do.call(rbind,RES))

write.csv(DF,'Human_modules.csv')

### generate regulation table ##
P <-unique(DF$TF)  ## unique pathway
G <-unique(DF$Genes)  ## unique genes

Regulat <-data.frame(matrix(0,nrow=length(P),ncol=length(G)))
row.names(Regulat)<-as.character(P)
colnames(Regulat) <-as.character(G)

a <-colnames(Regulat)

for (i in 1:length(P)){
		Regulated <-DF[c(which(DF$TF %in% P[i])),2]  ## extract genes regulated by ith TF
		colindex <-which(a %in% as.character(Regulated))  ## extract the column index
		Regulat[i,c(colindex)]=1  ## if TF regulate that gene, denote it as 1
	}
write.table(Regulat,"Human.csv",sep=",",row.names=T,col.names=T)

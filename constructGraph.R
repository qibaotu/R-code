## scRNA-seq DATA #
setwd("/home/xiej/Juan/Experiment/Yan/KL/7")
files <- list.files(path="/home/xiej/Juan/Experiment/Yan/KL/7/Conds/", pattern="Yan*", full.names=T, recursive=FALSE)

for (i in 1:length(files)){
	BC <- readLines(files[i])  # (i-1)th BC
	BCnum <-length(BC)
	newTxt <- unlist(strsplit(BC, split = " "))  # extract all the conditions
	conds <-unique(newTxt[newTxt!=""])   # unique conditions
	if (BCnum >2 & length(conds==90)){   # only consider results that with more than 2 BCs and cover all conditions
		CONDS <-as.character()   # store the conditions 
		label_C <-as.numeric()   # store the occurence of one condistions
		for (j in 1:BCnum){
		BCcond <-unlist(strsplit(BC[j], split = " "))
		BCcond <-BCcond[BCcond!=""]  # exclude the blank string
		CONDS <-c(BCcond,CONDS)
		label_C <-c(label_C,rep(j,length(BCcond)))
		}
	df_C <-data.frame(conds=CONDS,label=label_C)
	uniq_C <-df_C$conds[!duplicated(df_C$conds)]   # unique conditions
	Node <-t(combn(uniq_C,2))
	
	Wt <-rep(-1,dim(Node)[1])
	for (k in 1:dim(Node)[1]){
		member1 <-df_C[which(df_C$conds %in% Node[k,1]),]   # identify which BC the k th Node appear
		member2 <-df_C[which(df_C$conds %in% Node[k,2]),]
		Wt[k] <-length(intersect(member1[,2],member2[,2])) # the weight between two node
	}
	GRAPH <-data.frame(Node[,1],Node[,2],Wt)
	if (dim(GRAPH)[1]!=0)	{
		write.csv(subset(GRAPH,Wt!=0,select=Node...1.:Wt),paste(basename(files[i]),".csv",sep=","),row.names=F)
		}
	}
}

print("Done")

	


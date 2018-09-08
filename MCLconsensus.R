### 

setwd('C:/Users/Juan.XIe/Desktop/test1/Consensus')
library(igraph)
library(mclust)
library(MCL)

target <-read.table("Yan_cell.csv",header=T,sep=",")   # ground truth cell type
CellNum <-dim(target)[1]  # the number of cells 
files <-list.files(path="CSV/", pattern="*csv", full.names=T, recursive=FALSE)

ARI <-list()   # a list to store the max ARI for each file
sorted_LABEL <-list()  # a list to store the cell label that led to max ARI for each file
INFLA <-list()
FileID <-as.numeric()  # a vector to store the fileID
CLUSTERNum <-as.numeric()  # a vector to store the # of clusters

for (J in 1:length(files)){
	Graph <-read.csv(files[J],header=T,sep=",")
	names(Graph) <-c("Node1","Node2","Weight")
	G <-graph.data.frame(Graph,directed = FALSE)  # convert file into graph
	A <- as_adjacency_matrix(G,type="both",attr="Weight",names=TRUE)  # convert graph into adjacency matrix
	MCL_name <-rownames(A)   # the vertix
	
	if (length(rownames(A)) ==CellNum){   # only consider the case that cover all cells
		CLUST <-list()
		for (i in 1:100){
		CLUST[[i]] <-mcl(A,addLoops = F,inflation =i,max.iter=200)
		}
		K <- as.data.frame(do.call(rbind,lapply(CLUST,'[[',1)))  # extract the number of clusters
		CAN_I <-c(which(as.numeric(as.character(K$V1))>=2))  # results that has more than 5 clusters
		tt <-as.numeric(as.character(K$V1))
		tt <-sort(table(tt),decreasing=T)[1]
		Final_K <-as.numeric(names(tt))		
		if(length(CAN_I!=0)){   # some may have less than 5 clusters, just to make sure
			ARI_sc <-rep(-1,length(CAN_I))
			MATRIX <-rep(0,CellNum)%o%rep(0,CellNum)
			for (k in 1:length(CAN_I)){
				MCL_label <-CLUST[[CAN_I[k]]]$Cluster  # record the label
				ClusterNum <-unique(MCL_label)   # record the number of clusters
				TEMP <-rep(0,CellNum)%o%rep(0,CellNum)
				temp <-rep(0,CellNum) %o% rep(0,length(ClusterNum))
				  for (n in 1:length(ClusterNum)){
				    index <-which(MCL_label==ClusterNum[n])
				    temp[index,n] <-1
				   TEMP <-TEMP+temp[,n]%o%temp[,n] 
				  }
				MATRIX <-MATRIX+TEMP
			}
			MATRIX <-MATRIX/length(CAN_I)
		}
	}
	rownames(MATRIX) <-colnames(MATRIX) <-rownames(A)
	hc <-hclust(dist(MATRIX))
	memb <-cutree(hc,k=Final_K)		

 	sorted <-memb[match(target$Cell_type,names(memb))]
 	ARI[[J]] <-adjustedRandIndex(sorted,target$Cluster)
	FileID <-c(FileID,J)
	CLUSTERNum <-c(CLUSTERNum ,Final_K)		
}

ARI_all <-as.data.frame(do.call(rbind,ARI))
DF <- data.frame(ARI=ARI_all,FileID =basename(files[FileID]),K=CLUSTERNum)
write.csv(DF,"Yan_split_KL_ARI_I.csv")			
setwd("/home/xiej/Juan/Experiment/Yan/KL")
library(igraph)
library(mclust)
library(MCL)
target <-read.table("/home/xiej/Juan/Experiment/Yan/Yan_cell.csv",header=T,sep=",")   # ground truth cell type
CellNum <-dim(target)[1]  # the number of cells 
files <-list.files(path="/home/xiej/Juan/Experiment/Yan/KL/Graph/", pattern="*csv", full.names=T, recursive=FALSE)

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
		for (i in 1:200){
		CLUST[[i]] <-mcl(A,addLoops = F,inflation =i,max.iter=200)
		}
		K <- as.data.frame(do.call(rbind,lapply(CLUST,'[[',1)))  # extract the number of clusters
		CAN_I <-c(which(as.numeric(as.character(K$V1))>=3))  # results that has more than 3 clusters
		
		if(length(CAN_I!=0)){   # some may have less than 3 clusters, just to make sure
			ARI_sc <-rep(-1,length(CAN_I))
			for (k in 1:length(CAN_I)){
				MCL_label <-CLUST[[CAN_I[k]]]$Cluster  # record the label
				ClusterNum <-CLUST[[CAN_I[k]]]$K   # record the number of clusters
				df_cell_label <-data.frame(cell=MCL_name,cluster=MCL_label,K=ClusterNum)
				sorted_label <-df_cell_label[match(target$Cell_type,df_cell_label$cell),]$cluster   # sort the label based on the order of target sample ID
				ARI_sc[k] <-adjustedRandIndex(sorted_label,target$Cluster)  # calculate ARI
			}
			BESTI <- CAN_I[c(which(ARI_sc==max(ARI_sc)))][1]   # find which CLUST results highest ARI, if multiple, choose the first one
			LABEL <- CLUST[[BESTI]]$Cluster
			INFLA[[J]] <-BESTI   # the value of parameter inflaction for that case
			# df_cell_label <-data.frame(cell=MCL_name,cluster=LABEL)
			# sorted_LABEL[[J]] <-df_cell_label[match(target$Cell_type,df_cell_label$cell),]$cluster 
			ARI[[J]] <- max(ARI_sc)
			FileID <-c(FileID,J)
			K <- CLUST[[BESTI]]$K
			CLUSTERNum <-c(CLUSTERNum,K)	
		}	
	}
}

ARI_all <-as.data.frame(do.call(rbind,ARI))
INFLA_all <-as.data.frame(do.call(rbind,INFLA))
DF <- data.frame(ARI=ARI_all,I =INFLA_all,FileID =basename(files[FileID]),K=CLUSTERNum)
write.csv(DF,"Yan_split_KL_ARI_I.csv")






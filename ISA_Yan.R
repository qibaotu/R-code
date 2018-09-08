### ISA biclustering 

library(isa2)
library(biclust)

setwd("/pylon5/cc5fpcp/xiej/Yan_GeneEnrich/ISA")
Yan <-read.table("/pylon5/cc5fpcp/xiej/Yan_GeneEnrich/ISA/Yan_expression.csv",header=T,sep=',')

matrix.please<-function(x) {
    m<-as.matrix(x[,-1])
    rownames(m)<-x[,1]
    m
}

M <-matrix.please(Yan)

C <-seq(10,600,10)  # thr.col
for (i in C){
			set.seed(i)
			res <-isa(M)   # ISA clustering
			bc <- isa.biclust(res)  # the number of BCs
			BCnum <-bc@Number
			GENES <-list()
			for (k in 1:BCnum){
				BC <-bicluster(M,bc,k)[[1]] 
				GENES[[k]] <-rownames(BC)
			}	
		writeLines(unlist(lapply(GENES, paste, collapse=" ")),paste('ISA_Yan',i,sep='_'))
		}
	}
	

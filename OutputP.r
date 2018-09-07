## hyper geometric 
# https://stackoverflow.com/questions/8382806/r-hypergeometric-test-phyper
#########################  ENrich2 output every regulator and corresponding P-value for each BC   ###
PATHS <-list.files(path="/pylon5/cc5fpcp/xiej/modules/",pattern="*csv",full.names=T,recursive=FALSE)

# loop through all files in a directory  #

files <- list.files(path="/pylon5/cc5fpcp/xiej/Yan_GeneModule/KLDual/Genes/1", pattern="*Yan*", full.names=T, recursive=FALSE)


ENrich2 <-function(Regulation,blocks){
	BCList <-list()  # list to store the results for each BC
	for (i in 1:length(blocks)){	# loop through each bicluster
		BC <- blocks[i]  # (i-1)th BC
		newTxt <- unlist(strsplit(BC, split = " "))
		#newTxt <-newTxt[grep('G',newTxt)]
		tBC <- newTxt
		tBC[tBC==""] <- NA	# replace the filling value (blank, "") with "NA" 
		tBC <- tBC[!is.na(tBC)] # remove the NAs
		tBC <-sapply(strsplit(tBC,split='_'),'[',1)
		numGN <- length(tBC)	# num. of gene in each BC
	
		tmplnc <-as.numeric()
		for(j in 1:numGN){ # loop through each gene
			GN <- tBC[j]	# pick up one gene in the bicluster "BC"
			#GN <-substr(GN,1,5)  # 2.0biclustering results b0040_1, we just need b0040
			if (!length(which(colnames(Regulation)%in% GN)) == 0){    # only process Genes that has TF 
				tb_GN_col <- Regulation[,GN]	# extract the column titled as 'GN'
				reg_ind <- which(tb_GN_col==1) # extract indices of all matched regulator by gene 'GN'
				tmplnc <-c(tmplnc,reg_ind)
			}
		}
	
## For some BC we cannot find corresponding regulators, thus need the length==0 judgement 
			if (!length(tmplnc)==0){
			LNC <-sort(unique(tmplnc))  # unique lncRNA for each BC
			# tmpX record the hyperP values for regulators in a BC
			# tmpM record the overall # of regulated gene for each regulator
			# tmpxx record the # of regulated gene for each regulator in the BC
			# tmpk record the BC size
			# tmpn record the regulon size
			tmpX <-tmpM <-tmpxx <-tmpk <-tmpn <-BCid <-BCNum <-as.numeric() 
			 
			for(k in 1:length(LNC)){ # loop through regulators to find number of genes of on regulator in the bicluster "BC"
					reg_genes <- Regulation[LNC[k],] 		
					reg_genes <- colnames(Regulation)[c(which(reg_genes==1))]  ## genes that regulated by this regulator
					matchedGenes <- reg_genes[reg_genes %in% tBC] 
					numMatched <- length(matchedGenes)	# number of genes of on regulator in the bicluster "BC"			
			
					m <-length(reg_genes)
					x <- numMatched
					K <-numGN
					n <-ncol(Regulation)-1-m
					tmp_P <-1-phyper(x-1, m,n,K)
					tmpX <-c(tmpX,tmp_P)
					tmpM <-c(tmpM,m)
					tmpxx <-c(tmpxx,x)
					tmpk <-c(tmpk,K)
					tmpn <-c(tmpn,n)
					BCid <-c(BCid,i)
					BCNum <-c(BCNum,length(blocks))
			}
			
			LNCRNA <-Regulation[c(LNC),1]  # the name of LncRNA 
			rst <-data.frame(Regulator=LNCRNA,Pvalue=tmpX,x=tmpxx,M=tmpM,n=tmpn,K=tmpk,BCID=BCid,BCNUM=BCNum)
			BCList[[i]] <- rst	
		}	
	}
	BCList	
}

### read functional pathways ##

LISTsplit <- function(LIST){
		P_by_BC <-split(LIST,LIST$BCID)  # split rst based on BCID
		bb <- lapply(P_by_BC,function(x) x[which(x$Pvalue==min(x$Pvalue)),])  # extract the minimum P-values
		rst_bb <- as.data.frame(do.call(rbind, bb))
		rst <- rst_bb[!duplicated(rst_bb$BCID),]   # several regulator may have the same smallest P value,just keep one
		rst		
}

lapply(files,function(x) {
		GeneBC <-readLines(x)  ## BC genes
		numBC <-length(GeneBC) 
		
		GETRSTlist <- function(PATH){
		Paths <- read.delim(PATH,sep=",",header=T)
		RST <- ENrich2(Paths,GeneBC)
		LS <-as.data.frame(do.call(rbind,RST))  # combine the lists into a dataframe
		LS
		}

		for (i in 1:length(PATHS)){
			LS <-GETRSTlist(PATHS[i])  # combine the lists into a dataframe
			rst <- LISTsplit(LS)
			write.csv(LS,paste(sep="_",x,basename(PATHS[i]),"all.csv"),row.names=F)
			#write.csv(rst,paste(sep="_",x,basename(PATHS[i]),".csv") ,row.names=F)
		}
})
print("Done")

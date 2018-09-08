## various scores  ##
fIn <-"/home/xiej/Juan/Experiment/FABIA/Yan/Genes/tmpCSV/001/"
files <- list.files(path=fIn, pattern="*csv", full.names=T, recursive=FALSE)

Percent <-function(seed){
sdstat <- matrix(0,nrow=nrow(seed),ncol=5)
for (i in 1:nrow(seed)){
	if (seed$Pvalue[i]<=0.05){
		sdstat[i,1]=1}
	if (seed$Pvalue[i]<=0.01){
		sdstat[i,2]=1}
	if (seed$Pvalue[i]<=0.001){
		sdstat[i,3]=1}
	if (seed$Pvalue[i]<=0.0001){
		sdstat[i,4]=1}
	if (seed$Pvalue[i]<=0.00001){
		sdstat[i,5]=1}
	}
return(sdstat)
}

LISTsplit <- function(LIST){
		P_by_BC <-split(LIST,LIST$BCID)  # split rst based on BCID
		bb <- lapply(P_by_BC,function(x) x[which(x$Pvalue==min(x$Pvalue)),])  # extract the minimum P-values
		rst_bb <- as.data.frame(do.call(rbind, bb))
		rst <- rst_bb[!duplicated(rst_bb$BCID),]   # several regulator may have the same smallest P value,just keep one
		rst		
}

PATHpercent <- function(P){
		DFsub <-as.character((Full_sub[which(Full_sub$Pvalue<=P),"Regulator"])) # subset those with significant Pvalues
		PERCENT <-length(unique(DFsub))/635
		PERCENT
		}

mm <-c(20,30,40)
CUTS <-c(0.05,0.01,0.001,0.0001,0.00001)

Pcutoffs <- Recall <-list()
for (i in 1:length(files)){
	FullList <-read.delim(files[i],sep=",",header=T)  
	BCnum <-FullList$BCNUM[1]   
	INList <-LISTsplit(FullList)  # each BC has one Pvalue
	
	## recall, count for the first 100 BCs or all BCs,whichever is smaller
	if (BCnum <100){
		Full_sub <-subset(FullList,BCID <=BCnum,select=Regulator:BCNUM)		
	}else{
		Full_sub <-subset(FullList,BCID <=100,select=Regulator:BCNUM)
	}
	Recall[[i]] <-mapply(PATHpercent,CUTS)
	
	## Precision, count for the first m BCs
	A <- list()
	k=0		
	for (m in mm){
		k <- k+1
		SUB <-subset(INList,BCID<=m,select=Regulator:BCNUM)  # select the first m th BC
		Subrow <-nrow(SUB)  # just in case that some biclusters don't have m rows
 		A[[k]] <- colSums(Percent(SUB)/Subrow)
	}				
	Pcutoffs[[i]] <- A
}

fout <- fIn

b <-as.data.frame(do.call(rbind, Recall))
names(b) <-CUTS
rownames(b) <-basename(files)
write.csv(b,paste(fout,'Recall.csv'),sep='')

for (j in 1:length(mm)){
	a <-as.data.frame(do.call(rbind,lapply(Pcutoffs,'[[',j)))
	names(a) <-CUTS
	rownames(a) <-basename(files)
	write.csv(a,paste(fout,mm[j],"Precision.csv",sep=""))
 }
		

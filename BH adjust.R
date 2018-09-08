### BH adjustment  ##
setwd('C:/Users/Juan.XIe/Desktop/test1/')
files <-list.files(getwd(),pattern='*csv')

for (i in 1:length(files)){
	F <-read.table(files[i],header=T,sep=',')
	F_by_BC <-split(F,F$BCID)   # split it by BCID
	bb <-lapply(F_by_BC,function(x) x$p.adj <-p.adjust(x$Pvalue,method='BH'))
	F$p.adj <-unlist(bb)
	write.csv(F,paste(basename(files[i]),'adjusted.csv',sep=''))
}
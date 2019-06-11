readMACS = function(filename){
	peaks=read.table(filename,header=T,fill=T,sep='\t')
	peaks$summit=peaks$summit+peaks$start
	peaks$p_value=10^(-peaks$X.10.log10.pvalue./10)
	peaks=peaks[,c('chr_no','start','end','length','summit','FDR...','p_value', 'tags', 'fold_enrichment')]
	peaks$MACS_SICER='MACS'
	colnames(peaks)[5:6]=c('summit/centre','FDR') # Change column names from FDR(%) to FDR and summit to summit/centre
	peaks
}

readSICER = function(filename){
	peaks=read.table(filename,fill=T)
	colnames(peaks)=c("chrom", "start", "end", "ChIP_island_read_count", "CONTROL_island_read_count","p_value", "fold_change", "FDR","chr_no","centre") #10 columns
	peaks$length=peaks$end-peaks$start+1
	peaks=peaks[,c('chr_no','start','end','length','centre','FDR','p_value')]
	peaks$MACS_SICER='SICER'
	colnames(peaks)[5]='summit/centre'
	peaks=peaks[which(peaks$chr_no!=-1),]
	peaks
}
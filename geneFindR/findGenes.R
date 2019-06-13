#This function takes an output peak file from either MACS or SICER (specified using the programme argument) which has been modified using the relevant chromosome_modifier file, and finds genes from a tab-delimited genefile, which is in a certain format (path hardcoded above), that are within a certain distance, dist, from each peak. It only includes peaks with a pvalue and FDR beneath the specified thresholds and writes a peaktable and genetable to tab-delimited files with the prefix output that can be read into Excel.
findGenes = function(peaksfile,programme,genefile,dist,pvalue,FDR,output){
	#Read in gene file
  print(programme)
	genes=read.table(genefile,header=T,sep='\t',quote="", fill = FALSE)
	
	#Read in peak file
	if(programme=='MACS'){
		peaks=readMACS(peaksfile)
	} else if(programme=='SICER'){
		peaks=readSICER(peaksfile)
	} else{
		print('Please provide the programme that was used to generate the peak file. This should either be "MACS" or "SICER".')
		}
		
	#Remove insignificant peaks.
	peaks=findSignificantPeaks(peaks,FDR,pvalue)
	#browser()
	#Find gene indices
	GeneInds=findGeneInds(peaks,genes,dist)
	
	#Create peak table and gene table
	peaktable=generatePeakTable(peaks,genes,GeneInds)
	genetable=generateGeneTable(peaks,genes,GeneInds,peaktable,output)
	
	#Write both tables to tab-delimited files which can be read into Excel
	writetoxls(peaktable,genetable,output)
	
}
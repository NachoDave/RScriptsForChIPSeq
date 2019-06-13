#This function removes peaks that have an FDR or pvalue that is higher than the threshold.
findSignificantPeaks = function(peaks,FDR,pvalue){
	peaks=peaks[which((peaks$FDR<=FDR)&(peaks$p_value<=pvalue)),]
	peaks
}

#This function finds the same Associate Gene Name in genes2 as in genes1 and adds its description entry to genes1's description column, overwriting anything that may have already been there
changeGeneDesc = function(genes1,genes2){
	genes1$Description=as.character(lapply(as.character(genes1$Gene_name),function(x) paste(as.character(genes2$Description[which(as.character(genes2$Gene_name)==x)]),collapse='; ')))
	genes1$Description[which(genes1$Description=='character(0)')]=''
	genes1
}

#This function finds gene indices for the genes that each peak is within distance bp of, separating them into three groups where the peak is either 'in gene' or 'upstream'/'downstream' of the gene.
findGeneInds = function(peaks,genes,dist){
	#in gene: gene start < peak summit < gene end
  ingene=apply(peaks,1,function(x) which((genes$Chromosome==as.numeric(x[1]))&(genes$Gene_start_bp<as.numeric(x[2]))&(genes$Gene_end_bp>as.numeric(x[3]))))
	#ingene=apply(peaks,1,function(x) which((genes$Chromosome==as.numeric(x[1]))&(genes$Gene.Start..bp.<as.numeric(x[2]))&(genes$Gene.End..bp.>as.numeric(x[3]))))
	#downstream: 0 < (gene start - peak end) < dist, or, 0 (gene start - peak summit) < dist
  downstream=apply( peaks,1,function(x) which( (genes$Chromosome==as.numeric(x[1])) & ( ( ((genes$Gene_start_bp-as.numeric(x[3]))<dist) & ((genes$Gene_start_bp-as.numeric(x[3]))>0) ) | ( ((genes$Gene_start_bp-as.numeric(x[5]))<dist) & ((genes$Gene_start_bp-as.numeric(x[5]))>0) ) ) ) )
	#downstream=apply( peaks,1,function(x) which( (genes$Chromosome.Number==as.numeric(x[1])) & ( ( ((genes$Gene.Start..bp.-as.numeric(x[3]))<dist) & ((genes$Gene.Start..bp.-as.numeric(x[3]))>0) ) | ( ((genes$Gene.Start..bp.-as.numeric(x[5]))<dist) & ((genes$Gene.Start..bp.-as.numeric(x[5]))>0) ) ) ) )
	#upstream: -dist < (gene end - peak start) < 0, or, -dist < (gene end - peak summit) < 0
  upstream=apply(peaks,1,function(x) which( (genes$Chromosome==as.numeric(x[1])) & ( ( ((genes$Gene_end_bp-as.numeric(x[2]))<0) & ((genes$Gene_end_bp-as.numeric(x[2]))>-dist) ) | ( ((genes$Gene_end_bp-as.numeric(x[5]))<0) & ((genes$Gene_end_bp-as.numeric(x[5]))>-dist) ) ) ))
  #upstream=apply(peaks,1,function(x) which( (genes$Chromosome.Number==as.numeric(x[1])) & ( ( ((genes$Gene.End..bp.-as.numeric(x[2]))<0) & ((genes$Gene.End..bp.-as.numeric(x[2]))>-dist) ) | ( ((genes$Gene.End..bp.-as.numeric(x[5]))<0) & ((genes$Gene.End..bp.-as.numeric(x[5]))>-dist) ) ) ))
	GeneInds=list(ingene,downstream,upstream)
	#browser()
	GeneInds
}

#This function produces a table where each row contains information about a peak along with its associated genes
generatePeakTable = function(peaks,genes,GeneInds){
	ingene=GeneInds[[1]]
	downstream=GeneInds[[2]]
	upstream=GeneInds[[3]]
	
	peakTable=array(NA,dim=c(length(peaks[,1]),13))
	peakTable=data.frame(peakTable)
	colnames(peakTable)=c(colnames(peaks),"Genes","Peak.Pos","Dist.to.Start")
  #browser()
	for(peak in 1:length(peaks[,1])){
		peakTable[peak, 2:10]=peaks[peak,2:10]
		
		if (peaks[peak, 1] == 998){
		  peakTable[peak, 1] = as.character('X')
		} else if (peaks[peak, 1] == 999){
		  peakTable[peak, 1] = as.character('Y')
		} else {
		  peakTable[peak, 1] = peaks[peak,1]
		}
		
		ingene_genes=as.character(genes$Gene_name[ingene[[peak]]])
		downstream_genes=as.character(genes$Gene_name[downstream[[peak]]])
		upstream_genes=as.character(genes$Gene_name[upstream[[peak]]])
		genelist=paste(c(ingene_genes,downstream_genes,upstream_genes),collapse=", ")
		ingene_tags=rep("in gene",length(ingene[[peak]]))
		downstream_tags=rep("downstream",length(downstream[[peak]]))
		upstream_tags=rep("upstream",length(upstream[[peak]]))
		taglist=paste(c(ingene_tags,downstream_tags,upstream_tags),collapse=", ")
		dists=paste(genes$Gene_start_bp[c(ingene[[peak]],downstream[[peak]],upstream[[peak]])]-peaks$summit[peak],collapse=", ")
		peakTable[peak,11:13]=c(genelist,taglist,dists)
	}
	
	peakTable
}

generateGeneTable = function(peaks,genes,GeneInds,peakTable,fileName){
	allGeneInds=sort(unique(unlist(GeneInds)))
	genelist=unlist(strsplit(as.character(peakTable$Genes),split=", "))
	poslist=unlist(strsplit(as.character(peakTable$Peak.Pos),split=", "))
	distlist=unlist(strsplit(as.character(peakTable$Dist.to.Start),split=", "))
  fileName=tail(strsplit(fileName,split='/')[[1]],n=1)
	
	geneTable=array(NA,dim=c(length(allGeneInds),11))
	geneTable=data.frame(geneTable)
  geneColNms = colnames(genes)
	colnames(geneTable)=c(geneColNms[c(6,3,1,4,5,7,2,8)],"Peak.Pos","Dist.to.Gene.Start","Peaks")
	#browser()
	for(gene in 1:length(allGeneInds)){
	  # Chromosome
	  if(genes$Chromosome[allGeneInds[gene]] == 998){
	    geneTable[gene, 1] = as.character('X')
	  } else if(genes$Chromosome[allGeneInds[gene]] == 999){
	    geneTable[gene, 1] = as.character('Y')
	  } else {
	  geneTable[gene, 1] = as.character((genes$Chromosome[allGeneInds[gene]]))
	  }
	  # Gene name
	  geneTable[gene, 2] = as.character((genes$Gene_name[allGeneInds[gene]]))
	  # Ensemble ID 
	  geneTable[gene, 3] = as.character((genes$Gene_stable_ID[allGeneInds[gene]]))
	  # Start bp
	  geneTable[gene, 4] = as.character((genes$Gene_start_bp[allGeneInds[gene]]))
	  # End bp
	  geneTable[gene, 5] = as.character((genes$Gene_end_bp[allGeneInds[gene]]))
	  # Strand
	  geneTable[gene, 6] = as.character((genes$Strand[allGeneInds[gene]]))
	  # Description
	  geneTable[gene, 7] = as.character((genes$Gene_description[allGeneInds[gene]]))
	  # Type
	  geneTable[gene, 8] = as.character((genes$Gene_type[allGeneInds[gene]]))

	  goi=as.character(genes$Gene_name[allGeneInds[gene]])
		#geneTable[gene,2]=goi
		#geneTable[gene,c(1,3:5)]=as.character(genes[allGeneInds[gene],c(1,3:5)])
		#geneTable[gene,6]=as.character(genes[allGeneInds[gene],6])
		peaks=paste(paste(fileName,'_',grep(paste('(^|[ ,]+)',goi,'([ ,]|$)',sep=''),peakTable$Genes),sep=''),collapse=", ")
		pos=paste(poslist[which(genelist==goi)],collapse=", ")
		dists=paste(distlist[which(genelist==goi)],collapse=", ")
		geneTable[gene,9:11]=c(pos,dists, peaks)

	}
	
	geneTable
}

#This function writes both peakTable and geneTable to tab-delimited files which can be read into Excel
writetoxls = function(peakTable,geneTable,outputPath){
  fileName=tail(strsplit(outputPath,split='/')[[1]],n=1)
  peakColNames=colnames(peakTable)
  peakTable$Peak.ID = paste(fileName,'_',rownames(peakTable),sep='')
  peakTable = peakTable[,c('Peak.ID',peakColNames)]
	write.table(peakTable,file=paste(outputPath,'_peakTable.xls',sep=''),sep="\t",row.names=F)
	write.table(geneTable,file=paste(outputPath,'_geneTable.xls',sep=''),sep="\t",row.names=F)
}

#This function compares two lists of genes to find common genes and unique genes to each list, writes them to file, and produces a crude venn diagram. A better venn diagram can be produced using the venn.diagram function from the library "VennDiagram".














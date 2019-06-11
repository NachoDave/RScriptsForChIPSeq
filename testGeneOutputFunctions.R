
source('/home/rstudio/scripts/findGeneReadFunctions.R')
source('/home/rstudio/scripts/findGenes.R')
source('/home/rstudio/scripts/findGeneFunctions.r')


peaksFileMod = '/home/rstudio/data/CurryExample/results/SRR6730206_SRR6730208ctrl_PstdMACS14_peaksmodified.xls'
geneFile = '/home/rstudio/scripts/geneTable/GeneCoordinates_ENSEMBL_Biomart150319sorted.txt'
outputPath = '/home/rstudio/data/CurryExample/results/trial_out'

genes=read.table(geneFile,header=T,sep='\t')
peaks = peaks=readMACS(peaksFileMod)

#gnInds = findGeneInds(peaks,genes,10000)

#ingene=apply(peaks,1,function(x) which((genes$Chromosome==as.numeric(x[1]))&(genes$Gene_start_bp<as.numeric(x[2]))&(genes$Gene_end_bp>as.numeric(x[3])))) #

#pkTb = generatePeakTable(peaks,genes,gnInds)

#gnTb = generateGeneTable(peaks, genes, gnInds, pkTb, 'deletethisfile')

writetoxls(pkTb,gnTb,outputPath)
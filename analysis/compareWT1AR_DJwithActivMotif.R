library("readxl")

geneListDJ_PCOS <-  read.table('/home/rstudio/data/WT1/results/genes/1_PCOS_AR_i7_geneTable.xls',header=T,fill=T,sep='\t')
geneListDJ_Fertile <- read.table('/home/rstudio/data/WT1/results/genes/2_Fertile_WT1_geneTable.xls',header=T,fill=T,sep='\t')

geneListActifMotive <- read_excel('/home/rstudio/data/WT1/data/AM Swansea AR WT1 ChIP-Seq analysis with Input 16428/2093Swansea_AR-WT1_genes.xlsx')

genesDJ = union(geneListDJ_Fertile$Gene_name, geneListDJ_PCOS$Gene_name)

genesIntersectDJ_AM = intersect(genesDJ, geneListActifMotive$`Gene Name`)
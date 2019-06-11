

library("readxl")
library(VennDiagram)

rm(list=ls())

# Load PCOS Data from Me and From actif motif------------------------------------------------------------------#
outP = 'hg19'

geneListDJ_PCOS <-  read.table('/home/rstudio/data2/WT1_300519/results/genes/1_PCOS_AR_i7_filteredhg19_trim_bowtie2UPBlkLstRm_MACS_geneTable.xls',header=T,fill=T,sep='\t')
geneListDJ_Fer <- read.table('/home/rstudio/data2/WT1_300519/results/genes/2_fertile_WT1_i6_filteredhg19_trim_bowtie2UPBlkLstRm_MACS_geneTable.xls',header=T,fill=T,sep='\t')

geneListActifMotive <- read_excel('/home/rstudio/data2/WT1_ActiveMotif/AM Swansea AR WT1 ChIP-Seq analysis with Input 16428/2093Swansea_AR-WT1_genes.xlsx')

geneListAM_PCOS <- geneListActifMotive[which(geneListActifMotive$`1_PCOS_AR::1 Present` == 1),]
geneListAM_Fer <- geneListActifMotive[which(geneListActifMotive$`2_fertile_WT1::1 Present` == 1),]
#genesDJ = union(geneListDJ_Fertile$Gene_name, geneListDJ_PCOS$Gene_name)

#Compare Active motif genes vs DJ genes PCOS cells -------------------------------#
# Get Intersect
DJAMOverlapPCOS = intersect(geneListAM_PCOS$`Gene Name`, geneListDJ_PCOS$Gene_name)
AMOnlyPCOS = setdiff(geneListAM_PCOS$`Gene Name`, geneListDJ_PCOS$Gene_name)
DJOnlyPCOS = setdiff(geneListDJ_PCOS$Gene_name,geneListAM_PCOS$`Gene Name`)
unionDJ_AM_PCOS = union(geneListAM_PCOS$`Gene Name`, geneListDJ_PCOS$Gene_name)

# Draw Venn diagram
require("VennDiagram")
pdf(paste("/home/rstudio/scripts/analysis/PCOS_DJvsAM",outP,".pdf", sep = ""),height=6,width=6)
grid.newpage()
draw.pairwise.venn(area1 = length(geneListDJ_PCOS$Gene_name), area2 = length(geneListAM_PCOS$`Gene Name`), cross.area = length(DJAMOverlapPCOS), category = 
                     c('DJ_PCOS', 'AM_PCOSDJ'),lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2))
dev.off()

# Tables of overlapping/nonoverlapping genes

presentDJPCOS = integer(length(unionDJ_AM_PCOS))
presentAMPCOS = integer(length(unionDJ_AM_PCOS))

for(dx in 1:length(unionDJ_AM_PCOS)){
  if (unionDJ_AM_PCOS[dx] %in% AMOnlyPCOS){
    presentAMPCOS[dx] = 1
  }  
  if (unionDJ_AM_PCOS[dx] %in% DJOnlyPCOS){
    presentDJPCOS[dx] = 1
  }
  if (unionDJ_AM_PCOS[dx] %in% DJAMOverlapPCOS){
    presentDJPCOS[dx] = 1
    presentAMPCOS[dx] = 1
  }
}

presDJAM = data.frame(unionDJ_AM_PCOS, presentDJPCOS, presentAMPCOS)

write.table(presDJAM,paste('/home/rstudio/scripts/analysis/PCOS_genes_DJvsAM',outP,'.xls', sep = ''),sep="\t",row.names=F)

#Compare Active motif genes vs DJ genes Fertile cells -------------------------------#
# Get intersect
DJAMOverlapFer = intersect(geneListAM_Fer$`Gene Name`, geneListDJ_Fer$Gene_name)
AMOnlyFer = setdiff(geneListAM_Fer$`Gene Name`, geneListDJ_Fer$Gene_name)
DJOnlyFer = setdiff(geneListDJ_Fer$Gene_name,geneListAM_Fer$`Gene Name`)
unionDJ_AM_Fer = union(geneListAM_Fer$`Gene Name`, geneListDJ_Fer$Gene_name)
# Draw Venn diagram
require("VennDiagram")
pdf(paste("/home/rstudio/scripts/analysis/Fer_DJvsAM",outP,".pdf"),height=6,width=6)
grid.newpage()
draw.pairwise.venn(area1 = length(geneListDJ_Fer$Gene_name), area2 = length(geneListAM_Fer$`Gene Name`), cross.area = length(DJAMOverlapFer), category = 
                     c('DJ_Fer', 'AM_Fer'),lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2))
dev.off()

# Tables of overlapping/nonoverlapping genes
presentDJFer = integer(length(unionDJ_AM_Fer))
presentAMFer = integer(length(unionDJ_AM_Fer))

for(dx in 1:length(unionDJ_AM_Fer)){
  if (unionDJ_AM_Fer[dx] %in% AMOnlyFer){
    presentAMFer[dx] = 1
  }  
  if (unionDJ_AM_Fer[dx] %in% DJOnlyFer){
    presentDJFer[dx] = 1
  }
  if (unionDJ_AM_Fer[dx] %in% DJAMOverlapFer){
    presentDJFer[dx] = 1
    presentAMFer[dx] = 1
  }
}

presDJAMFer = data.frame(unionDJ_AM_Fer, presentDJFer, presentAMFer)

write.table(presDJAMFer,paste('/home/rstudio/scripts/analysis/Fer_genes_DJvsAM',outP,'.xls', sep = ''),sep="\t",row.names=F)
## Compare PCOS vs Fertile for overlapping genes

PCOS_FerOverlapAM <- intersect(geneListAM_PCOS$`Gene Name`, geneListAM_Fer$`Gene Name`)
PCOS_FerOverlapDJ <- intersect(geneListDJ_PCOS$Gene_name, geneListDJ_Fer$Gene_name)

# Draw Venn diagram
require("VennDiagram")
pdf("/home/rstudio/scripts/analysis/FerPCOSOverlapAM.pdf",height=6,width=6)
grid.newpage()
draw.pairwise.venn(area1 = length(geneListAM_PCOS$`Gene Name`), area2 = length(geneListAM_Fer$`Gene Name`), cross.area = length(PCOS_FerOverlapAM), category = 
                     c('AM_PCOS', 'AM_Fer'),lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2))
dev.off()

require("VennDiagram")
pdf(paste("/home/rstudio/scripts/analysis/FerPCOSOverlapAMDJ", outP, ".pdf"),height=6,width=6)
grid.newpage()
draw.pairwise.venn(area1 = length(geneListDJ_PCOS$Gene_name), area2 = length(geneListDJ_Fer$Gene_name), cross.area = length(PCOS_FerOverlapDJ), category = 
                     c('DJ_PCOS', 'AM_PCOSDJ'),lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2))
dev.off()

# 


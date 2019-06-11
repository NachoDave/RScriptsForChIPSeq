
library(VennDiagram)

dat1= read.table('/home/rstudio/data2/WT1/results/genes/1_PCOS_AR_i7_geneTable.xls',header=T,fill=T,sep='\t')
dat2= read.table('/home/rstudio/data2/WT1/results/genes/2_Fertile_WT1_geneTable.xls',header=T,fill=T,sep='\t')


gn1 = as.vector(dat1$Gene_name)
gn2 = as.vector(dat2$Gene_name)

unionGn1Gn2 = union(gn1, gn2)
intsectGn1Gn2 = intersect(gn1, gn2)
inG1Only = setdiff(gn1, gn2)
inG2Only = setdiff(gn2, gn1)



noGnsCommon = length(intsectGn1Gn2)
noPCOSOnly = length(inG1Only)
noFertileWT1Only = length(inG2Only)

#draw.single.venn(area = 22, category = "Dog People")

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
require("VennDiagram")

pdf("/home/rstudio/data/PCOS_FERTILEWT1_OVERLAP.pdf",height=6,width=6)
grid.newpage()
  draw.pairwise.venn(area1 = 5528, area2 = 750, cross.area = length(intsectGn1Gn2), category = 
                       c('PCOS_AR', 'Fertile_WT1'),lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2))
  dev.off()
  
# Output list of common genes
  
l = max(c(noGnsCommon, noPCOSOnly, noFertileWT1Only))
X = vector(mode="character", length=l)
XPCOS = X
XPCOS[1:noGnsCommon] = intsectGn1Gn2
XWT = X
XWT[1:noFertileWT1Only] = inG2Only

diffExGns = data.frame('Expressed in both PCOS and Fertile' = X, 'PCOS only' = X, 'Fertile WT1 only' = X)
diffExGns$Expressed.in.both.PCOS.and.Fertile = XPCOS
diffExGns$PCOS.only = inG1Only
diffExGns$Fertile.WT1.only = XWT

write.table(diffExGns, '/home/rstudio/data/commonGenesFertileWT1_PCOSAR.xls',sep="\t",row.names=F)

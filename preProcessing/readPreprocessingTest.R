library(ShortRead) # import short read library

source('/home/rstudio/scripts/preProcessing/readTrimming.R') # import trimming functions



args = commandArgs(trailingOnly = TRUE)
fnIn1 <- "/home/rstudio//data/FABP4/Example/FABP4-FABP4_S5_L001_R1_001_10000reads.fastq.gz" # <- args[1] # input filename -this is the only compulsory argument
fnIn2 <- "/home/rstudio//data/FABP4/Example/FABP4-FABP4_S5_L001_R2_001_10000reads.fastq.gz"


# Check for filename arguments
# if (length(args) == 1){
#   fnOut <- sprintf("%s_trim.fastq", gsub('.fastq','', fnIn))
# } else {fnOut <- args[2]}
# 
# # Check for number of undefined bases
# if (length(args) < 3){
#   uDefB = 1
# } else {uDefB = as.integer(args[3])}
# 
# # Check for number of nucleotides below a certain Phred score
# if (length(args) < 4){
#   nBs = 5
# } else {nBs = as.integer(args[4])}
# 
# # Threshold on Phred score
# if (length(args) < 5){
#   thres = 30
# } else {thres = as.integer(args[5])}
thres = 30
nBs = 5
uDefB = 1

fnOut1 = "/home/rstudio//data/FABP4/Example/FABP4-FABP4_S5_L001_R1_001_10000readsTrim.fastq.gz"
fnOut2 = "/home/rstudio//data/FABP4/Example/FABP4-FABP4_S5_L001_R2_001_10000readsTrim.fastq.gz"

#cnts = fastqFilt(fnIn, thres, nBs, uDefB, fnOut)

cnts = fastFilt_PE(fnIn1, fnIn2, thres, nBs, uDefB, fnOut1, fnOut2)
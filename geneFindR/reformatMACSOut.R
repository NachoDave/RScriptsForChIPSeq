# Change the output from MACs to add an additional column with chromosome number
# remove the wierd chromosome columns and set X and Y to 998 and 999 respectively 
# in the new chromosome column

# If no FDR column present (i.e no ctrl) add one 

reformatMACSOut <- function(peaksFile){
  
  xx = read.table(peaksFile,header=T,fill=T,sep='\t')
  
  sze <- dim(xx)
  
  if (exists('xxx')){
    rm(xxx)
  }
  
  xxx <- xx
  
  if (!"FDR..." %in% colnames(xxx)){ # add the FDR... column
    xxx$FDR... <- 0
    
  }
  
  xxx$chr_no <-0 # add column of chromosome ints to dataframe (or set to 0 if chr_no already exists)
  
  #xxx$chr_no = as.integer(gsub("chr", "", xxx$chr))
  
  xxx$chr_no = (gsub("chr", "", xxx$chr))
  xxx$chr_no[which(xxx$chr_no == 'X' )] = 998
  xxx$chr_no[which(xxx$chr_no == 'Y' )] = 999
  
  validChr <- c(1:22, 998, 999) # allowable chromosome numbers
  xxx <- xxx[c(sapply(xxx$chr_no, function(.vec){any(.vec == validChr)} )), ]
  
  write.table(xxx, file = gsub(".xls", "_modified.xls", peaksFile), sep='\t', row.names = FALSE)
}
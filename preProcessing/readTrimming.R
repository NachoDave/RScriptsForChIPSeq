library("ShortRead")

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

fastqFilt <- function(fid, thres = 30, nBs = 5, uDefB = 1, destination=sprintf("%s_rmNsAndLowQul.fastq", gsub('.fastq','', fid)))
{
  
  if (substrRight(fid, 2) == 'gz'){
    cmprs <- TRUE
  } else {
    cmprs <- FALSE
  }
  
  print(cmprs)
  
  dest2 = sprintf("%s", gsub('.fastq', 'X.fastq', destination) )
  # Need to remove file if it already exists ------------------------------------ #
  
  if(file.exists(destination)){
    print(paste('Deleting file:', destination))
    file.remove(destination)
  }
  
  if(file.exists(dest2)){
    print(paste('Deleting file:', dest2))
    file.remove(dest2)
  } 
  
  #-------------------------------------------------------------------------------#
  cnts = c(0, 0, 0)
  
  stream <- FastqStreamer(fid) # open stream object
  on.exit(close(stream))
  
  repeat{
    
    fq <-yield(stream) # load in the files
    if (length(fq) == 0) # if nothing left to load from file
    {break}
    
    cnts[1] <- cnts[1] + length(fq)
    
    filter <- nFilter(threshold = 1) # create a filter to remove reads (1 N in read)
    fq <- fq[filter(fq)] # use indexing to remove reads
    
    # remove reads with N bases
    qcount <- rowSums(as(quality(fq), "matrix") <= thres) 
    
    fqx <- fq[qcount >= nBs] # Number of reads where all Phred scores <= thres
    fq <- fq[qcount < nBs] # Number of reads where all Phred scores > thres
    
    cnts[2] <- cnts[2] + length(fq)
    cnts[3] <- cnts[3] + length(fqx)
    #browser()
    writeFastq(fqx, dest2, "a", compress = cmprs)
    writeFastq(fq, destination, "a", compress = cmprs)
    
  }
  cnts
}

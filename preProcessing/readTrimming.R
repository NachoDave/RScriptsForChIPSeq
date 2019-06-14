library("ShortRead")

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# Function to trim single ended reads ==================================================#
fastqFilt <- function(fid, thres = 30, nBs = 5, uDefB = 0, destination=gsub('.fastq','_Trim.fastq', fid))
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
    
    filter <- nFilter(threshold = uDefB) # create a filter to remove reads (1 N in read)
    fq <- fq[filter(fq)] # use indexing to remove reads
    
    # remove reads with N bases
    
    #qcount <- rowSums(as(quality(fq), "matrix") <= thres) 
    
    #fqx <- fq[qcount >= nBs] # Number of reads where all Phred scores <= thres
    #fq <- fq[qcount < nBs] # Number of reads where all Phred scores > thres
    
    q <- as(quality(fq), "matrix") <= thres # find base pairs with Phred score > thrs
    q[is.na(q)] = FALSE # set any NA to FALSE (NAs occur when the read length varies between reads)
    
    qcount <- rowSums(q) # sum the rows to find number of bases that don't meet quality threshold
    
    qok <- qcount < nBs # get indeices of ok reads
    
    fqx <- fq[!qok] # Number of reads where all Phred scores <= thres
    fq <- fq[qok] # Number of reads where all Phred scores > thres
    
    cnts[2] <- cnts[2] + length(fq)
    cnts[3] <- cnts[3] + length(fqx)
    #browser()
    writeFastq(fqx, dest2, "a", compress = cmprs)
    writeFastq(fq, destination, "a", compress = cmprs)
    
  }
  cnts
}

# Function to trim pair ended reads ==================================================#

fastFilt_PE <- function(fid, fid2, thres = 30, nBs = 5, uDefB = 0, dest1 = gsub('.fastq','_Trim.fastq', fid), 
dest2 =  gsub('.fastq','_Trim.fastq', fid2))
{
  
  
  if (substrRight(fid, 2) == 'gz'){
    cmprs <- TRUE
  } else {
    cmprs <- FALSE
  }
  
  print(cmprs)
  
  dest1x = sprintf("%s", gsub('.fastq', 'X.fastq', dest1) )
  dest2x = sprintf("%s", gsub('.fastq', 'X.fastq', dest2) )
  # Need to remove file if it already exists ------------------------------------ #
  
  if(file.exists(dest1)){
    print(paste('Deleting file:', dest1))
    file.remove(dest1)
  }
  
  if(file.exists(dest2)){
    print(paste('Deleting file:', dest2))
    file.remove(dest2)
  } 

  if(file.exists(dest1x)){
    print(paste('Deleting file:', dest1x))
    file.remove(dest1x)
  }
  
  if(file.exists(dest2x)){
    print(paste('Deleting file:', dest2x))
    file.remove(dest2x)
  } 
  #-------------------------------------------------------------------------------#
  cnts1 = c(0, 0, 0)
  cnts2 = c(0, 0, 0)
  stream1 <- FastqStreamer(fid) # open stream object
  on.exit(close(stream1))
  stream2 <- FastqStreamer(fid2) # open stream object
  on.exit(close(stream2))  
  
  repeat{
    
    # Read in file in sections
    fq1 <-yield(stream1) # load in the files
    if (length(fq1) == 0) # if nothing left to load from file
    {break}
    fq2 <-yield(stream2) # load in the files
    if (length(fq2) == 0) # if nothing left to load from file
    {break}   
    
    # Count number of reads
    cnts1[1] <- cnts1[1] + length(fq1)
    cnts2[2] <- cnts2[2] + length(fq2)
    
    # Remove reads with N bases
    filter <- nFilter(threshold = 0) # create a filter to remove reads (1 N in read)
    filt1 <- filter(fq1) # filter on file 1
    filt2 <- filter(fq2) # filter on file2
    filt <- filt1 & filt2 # combined filter
    #browser()
    fq1 <-fq1[filt]
    fq2 <- fq2[filt]
    
    #browser()
    
    # Remove reads with x bases less than y Phred score
    
    # remove reads with N bases
    q1 <- as(quality(fq1), "matrix") <= thres
    q2 <- as(quality(fq2), "matrix") <= thres
    q1[is.na(q1)] = FALSE
    q2[is.na(q2)] = FALSE
    
    qcount1 <- rowSums(q1) 
    qcount2 <- rowSums(q2) 
    
    q1ok <- qcount1 < nBs
    q2ok <- qcount2 < nBs
    qok <- q1ok & q2ok
    
    fq1x <- fq1[!qok] # Number of reads where all Phred scores <= thres
    fq1 <- fq1[qok] # Number of reads where all Phred scores > thres
    fq2x <- fq2[!qok] # Number of reads where all Phred scores <= thres
    fq2 <- fq2[qok] # Number of reads where all Phred scores > thres
    
    cnts1[2] <- cnts1[2] + length(fq1)
    cnts1[3] <- cnts1[3] + length(fq1x)#
    cnts2[2] <- cnts2[2] + length(fq2)
    cnts2[3] <- cnts2[3] + length(fq2x)# 
    
    writeFastq(fq1x, dest1, "a", compress = cmprs)
    writeFastq(fq1, dest1x, "a", compress = cmprs)
    writeFastq(fq2x, dest2, "a", compress = cmprs)
    writeFastq(fq2, dest2x, "a", compress = cmprs)    
    #browser()
  }
  cnts1
}
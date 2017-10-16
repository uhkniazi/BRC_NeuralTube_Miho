# File: 09.2_importMethylationExtractorFork.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: an external script called via methylation extractor to process a file
# Date: 07/09/2017


# Function: f_oGRReadBismarkMethylExtractor
# Desc: as input it takes the name of the file (tab separated) created by methyl_extractor script
#       in bismark and the strand i.e. + for OT files (Original Top) and - for OB. returns the data
#       in a GRagnes sorted object
# Args: file: name of bismark file; strand: strand + or -
# Rets: single stranded, sorted GRanges object with a Methylated mcols boolean variable
f_oGRReadBismarkMethylExtractor = function(file, strand){
  # read the data from the tab separated file
  # also skip the first line as it has comments
  require(GenomicRanges)
  # define string splitting function
  f1 = function(str) strsplit(str, '\t', fixed=T)[[1]]
  gr = GRanges()
  # open the file 
  infile = file(file, 'r')
  input = readLines(infile, n = 1)
  options(warn = -1)
  while(TRUE){
    input = readLines(infile, n = 2000000)
    if (length(input) == 0) break; 
    ## create a array/data frame from the string
    df = vapply(input, f1, FUN.VALUE = character(5), USE.NAMES = F)
    ## mark strings with non methylated cytosines
    f = df[2,] != '-'
    # create a GRanges object
    grtemp = GRanges(df[3,], IRanges(as.numeric(df[4,]), as.numeric(df[4,])), strand=strand)
    # add methylated flag
    grtemp$Methylated = f
    gr = append(gr, grtemp)
  }
  close(infile)
  options(warn = 0)
  gr = sort(gr)
  # garbage collector
  gc(verbose = FALSE)
  return(gr)  
} # function

args = commandArgs(trailingOnly = TRUE)

gr = f_oGRReadBismarkMethylExtractor(args[1], '*')
save(gr, file=paste0('~/Data/Temp/', args[1], '.rds'))
print('done')


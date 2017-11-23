# File: header.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: global variables
# Date: 08/11/2017


## variables
g_pid = 8
#g_did = 20
gcswd = getwd()
gcRemoteDir = "/run/user/1000/gvfs/sftp:host=10.202.64.29,user=k1625253/users/k1625253/brc_scratch/Data/ProjectsData/BRC_NeuralTube_Miho/"

p.old = par()

###### utility functions
## reads the report.txt file produced by bismark and bowtie2 and extracts some information
## about the alignments
dfParseBismarkReport = function(reportFile, title=''){
  # read the data in *report.txt file from bismark
  infile = file(reportFile, 'rt')
  input = readLines(infile, n = -1)
  # keep the lines with a : symbol in it as those containing information
  i = grep(':\t', input)
  # split the string based on this marker
  input = input[i]
  ## extract the relevant sections
  iIndices = c(1, 2, 3, 4, 12, 13, 14, 15, 17, 18, 19, 21, 22, 23)
  input = input[iIndices]
  lInput = strsplit(input, ':\t')
  close(infile)
  df = do.call(rbind, lInput)
  rownames(df) = df[,1]
  df = data.frame(df[,-1])
  colnames(df) = title
  iIndices = c(1, 2, 4:11 )
  df[,1] = as.character(df[,1])
  df[iIndices,1] = signif(as.numeric(df[iIndices,1])/1e6, 4)
  ## remove white space
  df[,1] = gsub(" ", "", df[,1], fixed = T)
  return(df)
}

f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}


# Function: f_bin_vector
# DESC: Takes the start and end values of the vector; and the number of
#       bins, and returns a data.frame with the start and ends of the 
#       bins
# ARGS: start coordinate, end coordinate, and number of bins
# RETS: a data.frame object with starts and ends of each bin
f_bin_vector = function(start, end, bins){
  s = floor(seq(start, end, length.out=bins+1))
  e = s-1
  e[length(e)] = s[length(s)]
  length(s) = length(s)-1
  e = e[2:length(e)]
  return(data.frame(start=s, end=e))
}# f_bin_vector


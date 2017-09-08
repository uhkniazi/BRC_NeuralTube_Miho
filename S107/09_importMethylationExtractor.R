# File: 09_importMethylationExtractor.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: serial import of methylation extractor output and convert to GRanges object
# Date: 07/09/2017


## set variables and source libraries
source('header.R')

library(GenomicRanges)

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

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
# get the query
g_did
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.group3, Sample.title, File.* from Sample, File
           where (Sample.idData = 14) AND (File.idSample = Sample.id AND File.type like "%duplicates removed%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)

## for each sample create the relevant files 
f1 = function(sid){
  i = which(dfSample$sid == sid)
  # get file name
  x = gsub('.bam$', '_sort2.txt.gz', dfSample$name[i])
  cpg.ot = c(idSample=sid, name=paste0('CpG_OT_', x), type='CpG_OT', group1='CpG Original Top Strand')
  cpg.ob = c(idSample=sid, name=paste0('CpG_OB_', x), type='CpG_OB', group1='CpG Original Bottom Strand')
  chh.ot = c(idSample=sid, name=paste0('CHH_OT_', x), type='CHH_OT', group1='CHH Original Top Strand')
  chh.ob = c(idSample=sid, name=paste0('CHH_OB_', x), type='CHH_OB', group1='CHH Original Bottom Strand')
  chg.ot = c(idSample=sid, name=paste0('CHG_OT_', x), type='CHG_OT', group1='CHG Original Top Strand')
  chg.ob = c(idSample=sid, name=paste0('CHG_OB_', x), type='CHG_OB', group1='CHG Original Bottom Strand')
  df = rbind(cpg.ot, cpg.ob, chh.ot, chh.ob, chg.ot, chg.ob)
  return(df)
}

dfNewData = lapply(dfSample$sid, f1)
dfNewData = do.call(rbind, dfNewData)
dfNewData = data.frame(dfNewData)

#### set working directory to appropriate location with methylation extractor files
setwd('~/Downloads/Temp/S107/')
csFiles = list.files('.', pattern = '*.gz$')

# sanity check
table(csFiles %in% dfNewData$name)
dim(dfNewData)

## create entry in database for these new files
## comment out as this has been done once
## dbWriteTable(db, name = 'File', value=dfNewData, append=T, row.names=F)

## process each file at a time after loading from database
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.group3, Sample.title, File.* from Sample, File
           where (Sample.idData = 14) AND (File.idSample = Sample.id AND File.type like "%CpG%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
dbDisconnect(db)
## CpGs first
lResults = lapply(dfSample$name, f_oGRReadBismarkMethylExtractor, '*')

## take the results in pairs OT then OB 
## assign strands and merge
## all the even ones are minus strands and odd ones are plus
lResults = lapply(seq_along(lResults), function(x){
  bOdd = ifelse((x %% 2) == 0, '-', '+')  
  strand(lResults[[x]]) = bOdd
  return(lResults[[x]])
})

## merge ranges from same sample
oGRLmerged = GRangesList()
options(warn=-1)
i = 1;
while(i <= nrow(dfSample)){
  gr = do.call(append, lResults[i:(i+1)])
  oGRLmerged = append(oGRLmerged, GRangesList(gr))
  i = i+2;
}
options(warn=0)

oGRLmerged = sort(oGRLmerged)
oGRLmerged.chg = oGRLmerged
names(oGRLmerged.chg) = unique(dfSample$sid)
rm(oGRLmerged)
rm(lResults)
gc(verbose = F)

## this should be one file for CpG or Chh or Chg methylation
setwd(gcswd)
n = make.names(paste('GRangesList object for CpG methylation from bs seq S107 miho rds'))
metadata(oGRLmerged.chg) = dfSample
n2 = paste0('~/Data/MetaData/', n)
save(oGRLmerged.chg, file=n2)
rm(oGRLmerged.chg)
gc(verbose = F)
# comment out as this has been done once
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
                comment='GRangesList object for CpG methylation from S107 run for bs seq data for miho ishida produced by the methylation extractor script from bismark')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

#################### repeat from Chh methylation
setwd('~/Downloads/Temp/S107/')

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
## process each file at a time after loading from database
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.group3, Sample.title, File.* from Sample, File
           where (Sample.idData = 14) AND (File.idSample = Sample.id AND File.type like "%CHH%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
dbDisconnect(db)
## CHH second
lResults = lapply(dfSample$name, f_oGRReadBismarkMethylExtractor, '*')

## take the results in pairs OT then OB 
## assign strands and merge
## all the even ones are minus strands and odd ones are plus
lResults = lapply(seq_along(lResults), function(x){
  bOdd = ifelse((x %% 2) == 0, '-', '+')  
  strand(lResults[[x]]) = bOdd
  return(lResults[[x]])
})

## merge ranges from same sample
oGRLmerged = GRangesList()
options(warn=-1)
i = 1;
while(i <= nrow(dfSample)){
  gr = do.call(append, lResults[i:(i+1)])
  oGRLmerged = append(oGRLmerged, GRangesList(gr))
  i = i+2;
}
options(warn=0)

oGRLmerged = sort(oGRLmerged)
oGRLmerged.chh = oGRLmerged
names(oGRLmerged.chh) = unique(dfSample$sid)
rm(oGRLmerged)
rm(lResults)
gc(verbose = F)

## this should be one file for CpG or Chh or Chg methylation
setwd(gcswd)
n = make.names(paste('GRangesList object for CHH methylation from bs seq S107 miho rds'))
metadata(oGRLmerged.chh) = dfSample
n2 = paste0('~/Data/MetaData/', n)
save(oGRLmerged.chh, file=n2)
rm(oGRLmerged.chh)
gc(verbose = F)

# comment out as this has been done once
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
                comment='GRangesList object for CHH methylation from S107 run for bs seq data for miho ishida produced by the methylation extractor script from bismark')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)


################ repeat for CHG methylation
setwd('~/Downloads/Temp/S107/')

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
## process each file at a time after loading from database
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.group3, Sample.title, File.* from Sample, File
           where (Sample.idData = 14) AND (File.idSample = Sample.id AND File.type like "%CHG%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
dbDisconnect(db)
## CHG third
lResults = lapply(dfSample$name, f_oGRReadBismarkMethylExtractor, '*')

## take the results in pairs OT then OB 
## assign strands and merge
## all the even ones are minus strands and odd ones are plus
lResults = lapply(seq_along(lResults), function(x){
  bOdd = ifelse((x %% 2) == 0, '-', '+')  
  strand(lResults[[x]]) = bOdd
  return(lResults[[x]])
})

## merge ranges from same sample
oGRLmerged = GRangesList()
options(warn=-1)
i = 1;
while(i <= nrow(dfSample)){
  gr = do.call(append, lResults[i:(i+1)])
  oGRLmerged = append(oGRLmerged, GRangesList(gr))
  i = i+2;
}
options(warn=0)

oGRLmerged = sort(oGRLmerged)
oGRLmerged.chg = oGRLmerged
names(oGRLmerged.chg) = unique(dfSample$sid)
rm(oGRLmerged)
rm(lResults)
gc(verbose = F)

## this should be one file for CpG or Chh or Chg methylation
setwd(gcswd)
n = make.names(paste('GRangesList object for CHG methylation from bs seq S107 miho rds'))
metadata(oGRLmerged.chg) = dfSample
n2 = paste0('~/Data/MetaData/', n)
save(oGRLmerged.chh, file=n2)
rm(oGRLmerged.chg)
gc(verbose = F)

# comment out as this has been done once
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
                comment='GRangesList object for CHG methylation from S107 run for bs seq data for miho ishida produced by the methylation extractor script from bismark')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)


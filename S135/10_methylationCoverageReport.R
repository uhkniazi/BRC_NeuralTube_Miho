# File: 10_methylationCoverageReport.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: import GRanges object created from methylation extractor and get coverage report
# Date: 7/11/2017


## set variables and source libraries
source('header.R')

library(GenomicRanges)

## choose the sequence names over which to calculate coverage
cvSeqnames = c(paste0('chr', 1:19), 'chrX', 'chrY')
bs.cut = 3 

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
# get the query
g_did
q = paste0('select * from MetaFile where (idData = 20) AND (comment like "%methylation%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)
n = paste0(dfSample$location, dfSample$name)
oGRLbis = f_LoadObject(n)
names(oGRLbis)

## one sample at a time load from the list 
oGRbis.Munique = oGRLbis[[6]]
iSampleID = 810;
#rm(oGRLbis); gc(reset = T)

oGRbis.Munique = oGRbis.Munique[seqnames(oGRbis.Munique) %in% cvSeqnames]
gc(reset = T)

### these steps have been performed in previous analysis to reduce size of data
# ## get the regions that are methylated
# oGRbis.M = oGRbis[oGRbis$Methylated]
# oGRbis.Munique = unique(oGRbis.M)
# 
# ## count the methylated regions i.e. oGRbis.Munique seen methylated and unmethylated
# ivMethylated = countOverlaps(oGRbis.Munique, oGRbis.M)
# ivUnmethylated = countOverlaps(oGRbis.Munique, oGRbis[!oGRbis$Methylated])
# ## sanity check with total coverage
# ivTotal = countOverlaps(oGRbis.Munique, oGRbis)
# 
# ## add this information to the object
# oGRbis.Munique$ivMethylated = ivMethylated
# oGRbis.Munique$ivUnmethylated = ivUnmethylated
# oGRbis.Munique$ivTotal = ivTotal
# rm(list = c('oGRbis', 'oGRbis.M', 'ivTotal', 'ivMethylated', 'ivUnmethylated'))

## if we want to remove low count potentially noisy data
f = oGRbis.Munique$ivTotal <= bs.cut
table(f)
oGRbis.Munique = oGRbis.Munique[!f]
gc(reset = T)

### save the information for this sample in a matrix
mCoverage = cbind(methylated=oGRbis.Munique$ivMethylated, unmethylated=oGRbis.Munique$ivUnmethylated, total=oGRbis.Munique$ivTotal)
mCoverage = cbind(mCoverage, prop=round(oGRbis.Munique$ivMethylated/oGRbis.Munique$ivTotal, 2))

## create a list to store this information
#lCoverage = list()
lCoverage$"810" = mCoverage
## repeat until all 8 samples done


# comment out as this has been done once
n = make.names(paste('List for CpG methylation coverage from bs seq S135 miho rds'))
n2 = paste0('~/Data/MetaData/', n)
#save(lCoverage, file=n2)
rm(list=c('oGRLbis', 'f', 'mCoverage'))
gc(verbose = F, reset = T)
## save to db
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='List object for CpG methylation coverage from S135 run for bs seq data for miho ishida produced by the methylation extractor script from bismark')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

### make plots 
s = seq(0.1, to = 1, by = 0.1)

## load the sample annotations
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'Sample')
# get the query
g_did
q = paste0('select * from Sample where (idData = 20)')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

pdf('results/ChgCoverageReport.pdf')

## plot at 4X coverage
par(mfrow=c(2,2))
temp = lapply(names(lCoverage), function(xn){
  x = lCoverage[[xn]][,'prop']
  x = x[x >= 0.1]
  sid = dfSample$title[dfSample$id == as.numeric(xn)]
  hist(x, breaks=s, ylab='Frequency', xlab='Proportion', main=paste('CpG Meth prop of 5mCs in', sid, 'at 4X'), xaxt='n')
  axis(1, at = s)
})

par(mfrow=c(2,2))
## at 10X
temp = lapply(names(lCoverage), function(xn){
  x = lCoverage[[xn]]
  x = x[x[,'total'] >= 10,]
  x = x[,'prop']
  x = x[x >= 0.1]
  sid = dfSample$title[dfSample$id == as.numeric(xn)]
  hist(x, breaks=s, ylab='Frequency', xlab='Proportion', main=paste('CpG Meth prop of 5mCs in', sid, 'at 10X'), xaxt='n')
  axis(1, at = s)
})

dev.off(dev.cur())







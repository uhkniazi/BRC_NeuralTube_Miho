# 02_combineReplicates.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: create a ranged object covering all peaks and add proportions data
# Date: 09/11/2017


source('header.R')

library(GenomicRanges)

## choose the sequence names over which to calculate coverage
cvSeqnames = c(paste0('chr', 1:19), 'chrX', 'chrY')

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
dbListFields(db, 'Data')
g_pid
df = dbGetQuery(db, 'select * from Data where Data.idProject = 8')
# get the query
g_did = df$id

q = paste0('select * from MetaFile where (idData = 14 OR idData = 15 OR idData = 20) AND (comment like "%Merged GRangesList%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

## process one sample at a time
n = paste0(dfSample$location[3], dfSample$name[3])
oGRLbis = f_LoadObject(n)
gc(reset = T)

## reduce to create one granges object
oGRbis = reduce(unlist(oGRLbis))
mMethylated = matrix(0, nrow = length(oGRbis), ncol = length(oGRLbis), dimnames = list(NULL, names(oGRLbis)))
mTotal = matrix(0, nrow = length(oGRbis), ncol = length(oGRLbis), dimnames = list(NULL, names(oGRLbis)))
mProp = matrix(0, nrow = length(oGRbis), ncol = length(oGRLbis), dimnames = list(NULL, names(oGRLbis)))
# find overlapping ranges in each element/sample of granges list
for (i in 1:length(oGRLbis)){
  df = findOverlaps(oGRbis, oGRLbis[[i]])
  mMethylated[queryHits(df), i] = oGRLbis[[i]][subjectHits(df)]$ivMethylated 
  mTotal[queryHits(df), i] = oGRLbis[[i]][subjectHits(df)]$ivTotal 
  mProp[queryHits(df), i] = round(oGRLbis[[i]][subjectHits(df)]$ivMethylated / oGRLbis[[i]][subjectHits(df)]$ivTotal, 3)
}

lMeta = list(methylated=mMethylated, total=mTotal, proportion=mProp, samples=metadata(oGRLbis))
metadata(oGRbis) = lMeta

save(oGRbis, file='results/GRangesMethylatedPositionsWithMetaData.rds')
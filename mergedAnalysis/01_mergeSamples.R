# File: 01_mergeSamples.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: merge the ranges and proportions data from all the samples
# Date: 08/11/2017


source('header.R')

library(GenomicRanges)

## choose the sequence names over which to calculate coverage
cvSeqnames = c(paste0('chr', 1:19), 'chrX', 'chrY')
bs.cut = 10 

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

q = paste0('select * from MetaFile where (idData = 14 OR idData = 15 OR idData = 20) AND (comment like "%GRangesList%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

## process one sample at a time
n = paste0(dfSample$location[3], dfSample$name[3])
oGRLbis = f_LoadObject(n)
names(oGRLbis)

## one sample at a time load from the list 
cSampleID = names(oGRLbis)[2]
#oGRbis = oGRLbis[[cSampleID]] # for S107 only
oGRbis.Munique = oGRLbis[[cSampleID]]

oGRbis.Munique = oGRbis.Munique[seqnames(oGRbis.Munique) %in% cvSeqnames]
gc(reset = T)

### these steps have been performed in previous analysis to reduce size of data for S126 and S135
# ## get the regions that are methylated
# oGRbis = oGRbis[seqnames(oGRbis) %in% cvSeqnames]
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
f = oGRbis.Munique$ivTotal < bs.cut
oGRbis.Munique = oGRbis.Munique[!f]
gc(reset = T)

#oGR.609 = oGRbis.Munique
#oGR.610 = oGRbis.Munique
#oGR.611 = oGRbis.Munique
#oGR.612 = oGRbis.Munique
#oGR.613 = oGRbis.Munique
#oGR.614 = oGRbis.Munique
#oGR.615 = oGRbis.Munique
#oGR.616 = oGRbis.Munique
#oGR.617 = oGRbis.Munique
#oGR.618 = oGRbis.Munique
#oGR.619 = oGRbis.Munique
#oGR.620 = oGRbis.Munique
#oGR.621 = oGRbis.Munique
#oGR.622 = oGRbis.Munique
#oGR.623 = oGRbis.Munique
#oGR.624 = oGRbis.Munique
#oGR.805 = oGRbis.Munique
#oGR.806 = oGRbis.Munique
#oGR.807 = oGRbis.Munique
#oGR.808 = oGRbis.Munique
#oGR.809 = oGRbis.Munique
#oGR.810 = oGRbis.Munique

## merge all these ranges into one list
oGRLbis = GRangesList(oGR.609, oGR.610, oGR.611, oGR.612, oGR.613, oGR.614, oGR.615, oGR.616, 
                      oGR.617, oGR.618, oGR.619, oGR.620, oGR.621, oGR.622, oGR.623, oGR.624,
                      oGR.805, oGR.806, oGR.807, oGR.808, oGR.809, oGR.810)

# set names equal to sample names
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dfSample = dbGetQuery(db, 'select * from Sample where Sample.idProject = 8')
cvSamples = as.character(dfSample$id)

names(oGRLbis) = cvSamples
metadata(oGRLbis) = dfSample

setwd(gcswd)
n = make.names(paste('GRangesList object for CpG methylation Merged from 3 bs seq runs miho rds'))
n2 = paste0('~/Data/MetaData/', n)
save(oGRLbis, file=n2)

# comment out as this has been done once
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
                comment='Merged GRangesList object for CpG methylation at 10X coverage from S107, S126 and S135 runs for bs seq data for miho ishida produced by the methylation extractor script from bismark')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)


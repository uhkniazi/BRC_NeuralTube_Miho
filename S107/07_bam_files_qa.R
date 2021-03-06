# File: 07_bam_files_qa.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the bam files
# Date: 15/03/2017


## set variables and source libraries
source('header.R')
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CBamQuality/master/CBamQuality.R'
download(url, 'CBamQuality.R')

# load the required packages
source('CBamQuality.R')
# delete the file after source
unlink('CBamQuality.R')

## choose the sequence names over which to calculate coverage
cvSeqnames = c(paste0('chr', 1:19), 'chrX', 'chrY')

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
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)

#### set working directory to appropriate location with bam files
setwd('DataExternal/bamFiles/')
csFiles = list.files('.', pattern = '*.bam$')
# check if these files match the file names in database
table(dfSample$name %in% csFiles)

csFiles = dfSample$name

## perform the analysis one sample at a time
## function to write the qa files
write.qa = function(bfn, sn){
  # open bam file
  bf = BamFile(bfn)
  lBam = lapply(sn, function(x){
    return(CBamScaffold(bf, x))
  })
  cat(paste('done', bfn, '\n'))
  return(lBam)
}

lAllBams = lapply(csFiles, function(x){
  write.qa(x, cvSeqnames)
})
names(lAllBams) = dfSample$id

## load additional metadata from bismark reports
csFiles = list.files('.', '*report.txt$')
## set the files in the same order as the samples
n = gsub('PE_report.txt', 'pe.bam_q10_sort_rd.bam', csFiles)
# sanity check
table(n %in% dfSample$name)
names(csFiles) = n
i = match(dfSample$name, n)
csFiles = csFiles[i]
identical(names(csFiles), dfSample$name)
names(csFiles) = dfSample$id

df = lapply(seq_along(csFiles), function(x) dfParseBismarkReport(csFiles[x], names(csFiles)[x]))
dfBismark = do.call(cbind, df)

setwd(gcswd)
n = make.names(paste('CBamScaffold rd picard bs seq S107 miho rds'))
lAllBams$metaBismark = dfBismark
lAllBams$meta = dfSample
lAllBams$desc = paste('CBamScaffold object from S107 run for bs seq data for miho ishida with quality 10 duplicates removed using picard tools', date())
n2 = paste0('~/Data/MetaData/', n)
save(lAllBams, file=n2)

# comment out as this has been done once
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='CBamScaffold object from S107 run for bs seq data for miho ishida with quality 10 duplicates removed and standard settings on trimmomatic
#                 with bismark alignment report as metadata. Processed using Picard tools instead of samtools.')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

### create the plots of interest
getwd()
lAllBams$desc = NULL
lAllBams$meta = NULL
dfBismark = lAllBams$metaBismark
lAllBams$metaBismark = NULL
names(lAllBams) = dfSample$title
colnames(dfBismark) = dfSample$title

pdf(file='Results/bam.q10.rd.qa.picard.pdf')
par(mfrow=c(2,2))
## returns the coverage for each scaffold in the bam file in matrix
f1 = function(ob){
  mat = sapply(ob, CBamScaffold.getCoverage)
  n = sapply(ob, CBamScaffold.getSeqname)
  colnames(mat) = n
  return(log(mat+1))
}

## call the function f1
lMat = lapply(lAllBams, f1)
i = seq_along(cvSeqnames)

## reorder the matrix, so that same chromosome/scaffold from all bam files are in one matrix
lMat.ordered = lapply(i, function(x){
  m = sapply(lMat, function(y) y[,x])
})

names(lMat.ordered) = cvSeqnames
iCol = rainbow(length(lAllBams))
## plot the lowess profile 
temp = sapply(cvSeqnames, function(x){
  m = apply(lMat.ordered[[x]], 2, function(y) lowess(y, f=2/10)$y )
  matplot(m, type='l', main=paste('Lowess fit to coverage', x), col=iCol, xlab='Bins', sub='Binned Coverage', ylab='Coverage', 
          lty=1:length(iCol))
})

par(mfrow=c(1,1))
plot.new()
legend('center', legend = names(lAllBams), ncol=3, col = iCol, lty=1:length(iCol), cex=0.6, lwd=2)


## average coverage, read width
## modify function with replace = T 
#CBamScaffold.getReadWidthSample = function(obj, size=1000) sample(obj@ivReadWidth, size = size, replace = T)

f1 = function(ob){
  mat = sapply(ob, CBamScaffold.getReadWidthSample)
  n = sapply(ob, CBamScaffold.getSeqname)
  colnames(mat) = n
  return(mean(colMeans(mat)))
}

## call the function f1
ivMat = sapply(lAllBams, f1)
barplot(ivMat, las=2, main='Average Read Width', ylab = 'Average Read Width', cex.names =0.6)

# number of reads aligned
f1 = function(ob){
  n = sapply(ob, function(x) length(CBamScaffold.getReadWidth(x)))
  n = sum(n)/1e+6
  return(n)
}

iReadCount = sapply(lAllBams, f1)

barplot(iReadCount, las=2, main='No. of reads aligned', ylab = 'No. of Reads in Millions', cex.names =0.8)
mReadCount = rbind(iReadCount, as.numeric(dfBismark[2,])*2)
rownames(mReadCount) = c('Post', 'Pre')

barplot(mReadCount, beside=T, las=2, main='No. of reads aligned', ylab = 'No. of Reads in Millions', cex.names =0.8,
        ylim=c(0, 250))
legend('topright', legend = c('Dup Rem', 'Orig'), fill=c('black', 'grey'))

# average binned coverage distribution
f1 = function(ob){
  mat = sapply(ob, getCoverageGammaParam)['shape',]
  n = sapply(ob, CBamScaffold.getSeqname)
  names(mat) = n
  return(mat)
}

## call the function f1
lMat = lapply(lAllBams, function(x) tryCatch(expr = f1(x), error=function(e)NULL))

boxplot(lMat, las=2, main='Average Binned Coverage', ylab = 'Average', outline=F, xaxt='n')
axis(1, at=1:length(lMat), labels = names(lMat), cex.axis=0.7, las=2)

#### average lowess profile for each group
f1 = function(ob){
  mat = sapply(ob, CBamScaffold.getCoverage)
  n = sapply(ob, CBamScaffold.getSeqname)
  colnames(mat) = n
  return(log(mat+1))
}

lMat = lapply(lAllBams, f1)
i = seq_along(cvSeqnames)

## reorder the matrix, so that same chromosome/scaffold from all bam files are in one matrix
lMat.ordered = lapply(i, function(x){
  m = sapply(lMat, function(y) y[,x])
})

names(lMat.ordered) = cvSeqnames
# choose the grouping to colour with
head(dfSample)
fGroups = factor(dfSample$group1)
levels(fGroups)
names(lAllBams)
# make sure the names and factors are in same order

iCol = rainbow(nlevels(fGroups))
## plot the lowess profile 
temp = sapply(cvSeqnames, function(x){
  m = lMat.ordered[[x]]
  m.av = sapply(levels(fGroups), function(l) {
    i = which(fGroups == l)
    rowMeans(m[,i])
  })
  m = apply(m.av, 2, function(y) lowess(y, f=2/10)$y)
  matplot(m, type='l', main=paste('Lowess fit to coverage', x), lwd=2, col=iCol, xlab='Bins', sub='Binned Coverage', ylab='Coverage', 
          lty=1)
  legend('bottomright', legend = levels(fGroups), ncol=1, col = iCol, lty=1, cex=1, lwd=2)
})
# 
# par(mfrow=c(1,1))
# plot.new()
# legend('center', legend = levels(fGroups), ncol=3, col = iCol, lty=1:length(iCol), cex=0.6, lwd=2)

dev.off(dev.cur())



# bf = BamFile('/run/user/1000/gvfs/sftp:host=10.202.64.28,user=k1625253/users/k1625253/Data/ProjectsData/BRC_Keloid/Aligned/S014/K2-2nd_S014_R1_.fastq.gz_q10_sort_rd.bam')
# si = seqinfo(bf)
# si = as.data.frame(si)
# si = si[order(si$seqlengths, decreasing = T),]
# sn = rownames(si[1:24,])
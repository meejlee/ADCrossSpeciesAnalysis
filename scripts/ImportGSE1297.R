## ImportGSE1297.R
# Download raw microarray files for GSE1297
# Finish with exporting normalized + annotated data

## Load libraries (Install if necessary)
library("GEOquery")
library("affy")
library("limma")
library(hgu133a.db)
library(hgu133acdf)

## Set Working Directory
mainwd <-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mainwd)

## Download and unzip dataset
setwd("../data/external")
getGEOSuppFiles("GSE1297")
setwd("../external/GSE1297")
untar("GSE1297_RAW.tar", exdir = "data")
## Import phenotype information
setwd("../..")
cels = read.csv("GSMList_GSE1297human.csv") #from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1297

setwd("../data/external/GSE1297")
sapply(paste("data", paste(cels$GSMID,".cel.gz",sep = ""), sep = "/"), gunzip)

## Normalize and annotate
setwd("../GSE1297/data")
data_raw = ReadAffy(verbose = FALSE, filenames = paste(cels$GSMID,".cel",sep = ""), cdfname = "hgu133acdf")
data_rmanorm = rma(data_raw)
data_df <- data.frame(exprs(data_rmanorm))

## Format hgu133a.db
hgu133Annot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), 
                          SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=", "), 
                          UNIGENE=sapply(contents(hgu133aUNIGENE), paste, collapse=", "), 
                          ENTREZID=sapply(contents(hgu133aENTREZID), paste, collapse=", "), 
                          ENSEMBL=sapply(contents(hgu133aENSEMBL), paste, collapse=", "), 
                          DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "))

## Annotate probes
comb_probes = merge(hgu133Annot, data_df, by.x = 0, by.y = 0, all = T)
comb_probes[1:4, 1:4]
write.table(comb_probes, file = "GSE1297_rmaprobes.txt", quote = FALSE, sep = "\t")

## OVERVIEW: Download + process GSE64398 data
# Download raw files from GEO
# Finish with exporting normalized + annotated data

## Load libraries (Install if necessary)
library("GEOquery")
library("affy")
library("beadarray")
library("illuminaMousev2.db")
library(mouse4302.db)
library(mouse4302cdf)

## Set Working Directory
mainwd <-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mainwd)

## Download dataset files
## ------- NOTE: issue with URL ------------------------------------------------------
# Manually download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64398
# Place in ../data/external/GSE64398
# getGEOSuppFiles("GSE64398") 
## -----------------------------------------------------------------------------------
setwd("../data/external/GSE64398") # Unzip downloaded files
untar("GSE64398_RAW.tar", exdir = "data")

## Import phenotype information and file name list
setwd("../..")
cels = read.csv("GSMList_GSE64398mouse.csv")
setwd("../data/external/GSE64398")
sapply(paste("data",paste(cels$geo_accession,cels$description,cels$sample.section.ch1,"Grn.idat.gz",sep = "_"), sep = "/"), gunzip)

## Import bead manifest file and annotate
## ------- NOTE: Need to download bead manifest file ---------------------------------
# Manually download bead manifest file from Illumina
# Place in ../data/external/GSE64398/data
## -----------------------------------------------------------------------------------
bgxfile = dir(path = "../GSE64398/data", pattern = ".bgx") #bead manifest file
setwd("../GSE64398/data")
data_raw = read.idat(paste(cels$geo_accession,cels$description,cels$sample.section.ch1,"Grn.idat",sep = "_"),bgxfile,annotation = "Symbol")
data_raw_vals = data_raw$E #pull matrix of raw intensities
data_raw_norm = nec(data_raw) #background correction
write.table(data_raw_norm, file = "GSE64398_raw_norm.txt", quote = FALSE, sep = "\t")

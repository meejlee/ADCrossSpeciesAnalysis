## ImportGSE48350.R
# Download raw files from GEO for GSE48350
# Finish with exporting normalized + annotated data

## Load libraries (Install if necessary)
library("GEOquery")
library("affy")
library("limma")
library(hgu133plus2.db)
library(hgu133plus2cdf)

## Set Working Directory
mainwd <-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mainwd)

## Download and unzip dataset
setwd("../data/external")
getGEOSuppFiles("GSE48350")
setwd("../external/GSE48350")
untar("GSE48350_RAW.tar", exdir = "data")

## Import phenotype
setwd("../..")
cels = read.csv("GSMList_GSE48350human.csv") # Sample info from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48350
setwd("../data/external/GSE48350")
sapply(paste(paste("data", cels$GSMID, sep = "/"),".CEL.gz",sep=""), gunzip)

## Normalize
setwd("../GSE48350/data")
data_raw = ReadAffy(verbose = FALSE, filenames = paste(cels$GSMID,".CEL",sep = ""), cdfname = "hgu133plus2cdf")
data_rmanorm = rma(data_raw)
data_df <- data.frame(exprs(data_rmanorm))

## Format hgu133plus2.db
hgu133plus2Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), 
                               SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), 
                               UNIGENE=sapply(contents(hgu133plus2UNIGENE), paste, collapse=", "), 
                               ENTREZID=sapply(contents(hgu133plus2ENTREZID), paste, collapse=", "), 
                               ENSEMBL=sapply(contents(hgu133plus2ENSEMBL), paste, collapse=", "), 
                               DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))

## Annotate probes
comb_probes = merge(hgu133plus2Annot, data_df, by.x = 0, by.y = 0, all = T)
comb_probes[1:4, 1:4]
write.table(comb_probes, file = "GSE48350_rmaprobes.txt", quote = FALSE, sep = "\t")

## DEG analysis
# Set up experimental design + fit model
dis <- factor(cels$Disease, levels = c("Control","AD"))
sex <- factor(cels$Gender, levels = c("Female","Male"))
age <- as.numeric(cels$Age)
design <- model.matrix(~dis + sex + age)
colnames(design) <- c("Int","Disease","Sex","Age")
rownames(design) <- cels$GSMID
fit <- lmFit(comb_probes[,8:69], design) #62 samples from hippocampus
# Generate category-specific contrasts
cont.matrix <- makeContrasts(Disease,Sex,Age, levels = design)
# COEF List:
# 1: Disease
# 2: Sex
# 3: Age
# Generate contrasts and save DEG list (DEG = adjusted p-value < 0.20)
fitdis <- contrasts.fit(fit,cont.matrix)
fitdis <- eBayes(fitdis)
fitdis2 <- data.frame(fitdis)
dtdis <- decideTests(fitdis, p.value = 0.20, adjust.method = "BH")
summary(dtdis)

res <- topTable(fitdis,coef=1,adjust.method="BH",n=15000,p.value=0.20,genelist=comb_probes$SYMBOL)
write.table(res, file="GSE48350degs_coef1_pvalonly_Dis.txt", quote=FALSE, sep="\t", col.names = NA) 
res <- topTable(fitdis,coef=2,adjust.method="BH",n=15000,p.value=0.20,genelist=comb_probes$SYMBOL)
write.table(res, file="GSE48350degs_coef2_pvalonly_Gen.txt", quote=FALSE, sep="\t", col.names = NA) 
res <- topTable(fitdis,coef=3,adjust.method="BH",n=15000,p.value=0.20,genelist=comb_probes$SYMBOL)
write.table(res, file="GSE48350degs_coef3_pvalonly_Age.txt", quote=FALSE, sep="\t", col.names = NA)

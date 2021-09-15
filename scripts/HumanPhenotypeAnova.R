## Load Dataset
library("readxl")

## Set Working Directory
mainwd <-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mainwd)
setwd("../modeling outputs")
file <- 'PCScoresForAnova_h1297m64398.xls'

## Before reading data, we will return the names of the sheets for later use:
sheets.name <- excel_sheets(file.path(file))

# GSE1297
hu1297_scores <- read_excel(file.path(file), sheet = sheets.name[1], col_names = FALSE)
colnames(hu1297_scores) <- paste0("PC_", 1:dim(hu1297_scores)[2])
hu1297mm64398_scores <- read_excel(file.path(file), sheet = sheets.name[2], col_names = FALSE)
colnames(hu1297mm64398_scores) <- paste0("PC_", 1:dim(hu1297mm64398_scores)[2])
hu1297_pheno <- read_excel(file.path(file), sheet = sheets.name[3], col_names = TRUE)
mm64398_pheno <- read_excel(file.path(file), sheet = sheets.name[4], col_names = TRUE)

## Anova
# GSE1297
PctExp = array(NA, dim = c(dim(hu1297_scores)[2], 4))
Pvalue = array(NA, dim = c(dim(hu1297_scores)[2], 4))
Pvalue.lm = array(NA, dim = c(dim(hu1297_scores)[2], 6))
tvalue.lm = array(NA, dim = c(dim(hu1297_scores)[2], 6))
for(index in 1 : dim(hu1297_scores)[2]){
  Dataframe <- data.frame(hu1297_scores[, index], hu1297_pheno)
  colnames(Dataframe) <- c('Reponse', 'disease', 'sex', 'age')
  Dataframe$disease <- factor(Dataframe$disease)
  Dataframe$sex <- factor(Dataframe$sex)
  fit = lm(Reponse ~ disease + sex + age, data = Dataframe)
  summary.fit <- summary(fit)
  Pvalue.lm[index, ] = summary.fit$coefficients[, 4]
  tvalue.lm[index, ] = summary.fit$coefficients[, 3]
  
  af <- anova(fit)
  afss <- af$"Sum Sq"
  res_2 <- cbind(af,PctExp=afss/sum(afss)*100)
  PctExp[index, ] = afss/sum(afss)*100
  Pvalue[index, ] = af$`Pr(>F)`
}
rownames(Pvalue.lm) <- paste0('PC.', 1:dim(hu1297_scores)[2])
colnames(Pvalue.lm) <- rownames(summary.fit$coefficients)
rownames(tvalue.lm) <- paste0('PC.', 1:dim(hu1297_scores)[2])
colnames(tvalue.lm) <- rownames(summary.fit$coefficients)
rownames(PctExp) <- paste0('PC.', 1:dim(hu1297_scores)[2])
colnames(PctExp) <- c("disease", "sex", "age", "residual")
rownames(Pvalue) <- paste0('PC.', 1:dim(hu1297_scores)[2])
colnames(Pvalue) <- c("disease", "sex", "age", "residual")

#Create a folder to save the results.
saved.dir <- './HumanPhenotypes'
if (!dir.exists(saved.dir)){
  dir.create((saved.dir))
}
write.table(tvalue.lm, file = "GSE1297_tvalueLM.txt", quote = FALSE, sep = "\t")
write.table(Pvalue.lm, file = "GSE1297_pvalueLM.txt", quote = FALSE, sep = "\t")
write.table(PctExp, file = "GSE1297_PctExp.txt", quote = FALSE, sep = "\t")
write.table(Pvalue, file = "GSE1297_Pvalue.txt", quote = FALSE, sep = "\t")

#GSE48350
file <- 'PCScoresForAnova_h48350m64398.xls'
## Before reading data, we will return the names of the sheets for later use:
sheets.name <- excel_sheets(file.path(file))

hu48350_scores <- read_excel(file.path(file), sheet = sheets.name[1], col_names = FALSE)
colnames(hu48350_scores) <- paste0("PC_", 1:dim(hu48350_scores)[2])
hu48350mm64398_scores <- read_excel(file.path(file), sheet = sheets.name[2], col_names = FALSE)
colnames(hu48350mm64398_scores) <- paste0("PC_", 1:dim(hu48350mm64398_scores)[2])
hu48350_pheno <- read_excel(file.path(file), sheet = sheets.name[3], col_names = TRUE)
hu48350_pheno$disease <- factor(hu48350_pheno$disease, levels = c("Control", "AD"))
hu48350_pheno$sex <- as.factor(hu48350_pheno$sex)

# GSE48350
PctExp = array(NA, dim = c(dim(hu48350_scores)[2], 4))
Pvalue = array(NA, dim = c(dim(hu48350_scores)[2], 4))
Pvalue.lm = array(NA, dim = c(dim(hu48350_scores)[2], 4))
tvalue.lm = array(NA, dim = c(dim(hu48350_scores)[2], 4))
for(index in 1 : dim(hu48350_scores)[2]){
  Dataframe <- data.frame(hu48350_scores[, index], hu48350_pheno)
  colnames(Dataframe) <- c('Reponse', 'disease', 'sex', 'age')
  Dataframe$disease <- factor(Dataframe$disease)
  Dataframe$sex <- factor(Dataframe$sex)
  fit = lm(Reponse ~ disease + sex + age, data = Dataframe)
  summary.fit <- summary(fit)
  Pvalue.lm[index, ] = summary.fit$coefficients[, 4]
  tvalue.lm[index, ] = summary.fit$coefficients[, 3]
  
  af <- anova(fit)
  afss <- af$"Sum Sq"
  res_2 <- cbind(af,PctExp=afss/sum(afss)*100)
  PctExp[index, ] = afss/sum(afss)*100
  Pvalue[index, ] = af$`Pr(>F)`
}
rownames(Pvalue.lm) <- paste0('PC.', 1:dim(hu48350_scores)[2])
colnames(Pvalue.lm) <- rownames(summary.fit$coefficients)
rownames(tvalue.lm) <- paste0('PC.', 1:dim(hu48350_scores)[2])
colnames(tvalue.lm) <- rownames(summary.fit$coefficients)
rownames(PctExp) <- paste0('PC.', 1:dim(hu48350_scores)[2])
colnames(PctExp) <- c("disease", "sex", "age", "residual")
rownames(Pvalue) <- paste0('PC.', 1:dim(hu48350_scores)[2])
colnames(Pvalue) <- c("disease", "sex", "age", "residual")

write.table(tvalue.lm, file = "GSE48350_tvalueLM.txt", quote = FALSE, sep = "\t")
write.table(Pvalue.lm, file = "GSE48350_pvalueLM.txt", quote = FALSE, sep = "\t")
write.table(PctExp, file = "GSE48350_PctExp.txt", quote = FALSE, sep = "\t")
write.table(Pvalue, file = "GSE48350_Pvalue.txt", quote = FALSE, sep = "\t")

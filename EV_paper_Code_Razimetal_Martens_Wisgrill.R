#### EV data proteomics ####
#### 13.12.2023 ####

# Stimulation of differentiated air-liquid interface cell cultures, generated from primary nasal epithelial cells of adults and term infants.
# Stimulation was conducted with with EV, EC, and MOCK.
# EV = Stimulation with vesicles
# EC = Stimulation with E.coli
# MOCK = unstimulated
# Following stimulation, proteins were extracted, TMT-labled and RPLC-MS/MS analysis was conducted.
# Raw data were processed and analyzed using Proteome Discoverer v 2.4 and the SwissProt protein database (human) with Mascot v 2.5.1 search engine. 

###############
### Library ###
###############

library(tidyverse)
library(edgeR)
library(ggfortify)
library(impute)
library(sva)
library(Biobase)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(org.Hs.eg.db)

#################
### Functions ###
#################

apply_tmm_factors <- function(y, color = NULL, plot = TRUE) {
  lib_facs <- mean(y$samples$lib.size) / y$samples$lib.size
  norm_facs <- lib_facs / y$samples$norm.factors
  cat("Overall Factors (lib.size+TMM):\n", sprintf("%-5s -> %f\n", colnames(y$counts), norm_facs))
  tmt_tmm <- as.data.frame(sweep(y$counts, 2, norm_facs, FUN = "*"))
  colnames(tmt_tmm) <- str_c(colnames(y$counts), "_tmm")
  if(plot == TRUE) {
    boxplot(log10(tmt_tmm), col = color, notch = TRUE, main = "TMM Normalized data")
  }
  tmt_tmm
}
CV <- function(df) {
  ave <- rowMeans(df)    # compute averages
  sd <- apply(df, 1, sd) # compute standard deviations
  cv <- 100 * sd / ave   # compute CVs in percent (last thing gets returned)
}
labeled_boxplot <- function(df, ylim, title) {
  cv = CV(df)
  boxplot(cv, ylim = c(0, ylim), notch = TRUE, main = title)
  text(x = 0.65, y = boxplot.stats(cv)$stats[3], 
       labels = round(boxplot.stats(cv)$stats[3], 1))
}

##############################
### Data Import & Trimming ###
##############################

#Import raw proteome matrix
mat <- read.csv("matrix_prot_ec.csv", header=TRUE, sep=";")
colnames(mat)[1] <- "Protein.IDs"

#Remove rows containing only NAs
mat <- na.omit(mat)
rownames(mat) <- mat$Protein.IDs
mat$Protein.IDs <- NULL

mat_corr <- mat

#Check if any residuals are left
summary(is.infinite(as.matrix(log10(mat_corr))))
summary(is.na(as.matrix(log10(mat_corr))))

#Create Metadata table
meta <- as.data.frame(c(rep("Adult_MOCK", 3), rep("Adult_EV", 3), rep("Adult_EC", 3),
                        rep("Term_MOCK", 3), rep("Term_EV", 3), rep("Term_EC", 3)))
colnames(meta)[1] = "group"
rownames(meta) = colnames(mat[,c(1:18)])
meta$group <- factor(meta$group, levels=c("Adult_MOCK", "Adult_EC", "Adult_EV", "Term_MOCK", "Term_EC", "Term_EV"))

#let's see what the starting data look like
color = c(rep("gray", 3), rep("red", 3), rep("blue", 3), rep("gray", 3), rep("red", 3), rep("blue", 3), rep("black",2))
boxplot(log2(mat_corr), col = color, notch = TRUE, main = "Starting MaxQuant data")
plotDensities(log2(mat_corr), col = c(rep("gray", 10),rep("red", 10)), main = 'Raw data', legend = F)
format(round(colSums(mat_corr), digits = 0), big.mark = ",")

##### SL Normalization #####
#with raw files
F11_raw <- mat_corr[,c(4,6,7,8,10,12,13,15,16,19)]
F12_raw <- mat_corr[,c(1,2,3,5,9,11,14,17,18,20)]

#with corr files norm_raw_pro
#F11_raw <- mat_corr[,c(4,5,7,9,10,12,13,16,18,19)]
#F12_raw <- mat_corr[,c(1,2,3,6,8,11,14,15,17,20)]
target <- mean(c(colSums(F11_raw), colSums(F12_raw)))
norm_facs <- target / colSums(F11_raw)
F11_sl <- sweep(F11_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(F12_raw)
F12_sl <- sweep(F12_raw, 2, norm_facs, FUN = "*")

# make a pre-IRS data frame after sample loading normalizations
data_sl <- cbind(F11_sl, F12_sl)

# see what the SL normalized data look like
boxplot(log2(data_sl), col = rep(c("gray", "red"), each = 10), 
        notch = TRUE, main = "Sample Loading (SL) normalized data: \nF8 (gray), F9 (red)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity')

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_sl), col = c(rep(c("gray", "red"), each = 10)), main = "SL normalization", legend=F)
format(round(colSums(data_sl), digits = 0), big.mark = ",")

##### TMM Normalization ######
# do TMM on raw data - we use it later for CV analyses
raw_tmm <- calcNormFactors(mat_corr)
data_raw_tmm <- sweep(mat_corr, 2, raw_tmm, FUN = "/") # raw data after TMM on original scale

# perform TMM on the SL-normed data and visualize resulting distributions
sl_tmm <- calcNormFactors(data_sl)
data_sl_tmm <- sweep(data_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale

boxplot(log2(data_sl_tmm), notch = TRUE, col = rep(c("gray", "red"), each = 10), 
        main = "TMM normalization of SL data\nF8 (gray), F9 (red)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity')

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_sl_tmm), col = rep(c("gray", "red"), each = 10), main = "SL/TMM normalization", legend = F)
format(round(colSums(data_sl_tmm), digits = 0), big.mark = ",")

# see how things cluster now that we have nice boxplots and density plots
plotMDS(log2(data_sl_tmm), col = rep(c("gray", "red"), each = 10), 
        main = "SL/TMM clusters group by TMT experiment")


# make new data frame with row sums from each frame
irs <- tibble(F11_sl$F11_pooled, F12_sl$F12_pooled)
colnames(irs) <- c("sum1", "sum2")

# get the geometric average intensity for each protein
irs$average <- apply(irs, 1, function(x) exp(mean(log(x))))

# compute the scaling factor vectors
irs$facF11 <- irs$average / irs$sum1
irs$facF12 <- irs$average / irs$sum2

# make new data frame with IRS normalized data
data_irs <- F11_sl * irs$facF11
data_irs <- cbind(data_irs, F12_sl * irs$facF12)

# see what the IRS data look like
boxplot(log2(data_irs), col = rep(c("gray", "red"), each = 10), 
        main = "Internal Reference Scaling (IRS) normalized data: \nF8 (gray), F9 (red)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity', notch = TRUE)

# can also look at density plots (like a distribution histogram)    
plotDensities(log2(data_irs), col = rep(c("gray", "red"), each = 10), main = "IRS data", legend = FALSE)


# this is data after SL, IRS, and TMM normalized on original scale
irs_tmm <- calcNormFactors(data_irs)
data_irs_tmm <- sweep(data_irs, 2, irs_tmm, FUN = "/") 

# see if box plots are aligned
boxplot(log2(data_irs_tmm), notch = TRUE, col = rep(c("gray", "red"), each = 10), 
        main = "TMM normalization of IRS data\nF8 (gray), F9 (red)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity')

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_irs_tmm), col = rep(c("gray", "red"), each = 10), main = "IRS/TMM data", legend=FALSE)

plotMDS(log10(data_irs_tmm[,-c(10,20)]), col = rep(c("gray", "red", "blue", "gray", "red", "blue"), each = 3), 
        main = "IRS/TMM clusters group by TMT experiment")

####################################
### Principal component analysis ###
####################################

meta$probands = c("Adult", "Adult", "Adult",
                  "Adult", "Adult", "Adult",
                  "Adult", "Adult", "Adult",
                  "Term", "Term", "Term",
                  "Term", "Term", "Term",
                  "Term", "Term", "Term")

meta$condition = c("MOCK", "MOCK", "MOCK",
                   "EV", "EV", "EV",
                   "EC", "EC", "EC",
                   "MOCK", "MOCK", "MOCK",
                   "EV", "EV", "EV",
                   "EC", "EC", "EC")

# index variable defines paired samples (those cell cultures that were derived from the same proband)
meta$index = c("1", "2", "3",
               "1", "2", "3",
               "1", "2", "3",
               "4", "5", "6",
               "4", "5", "6",
               "4", "5", "6")

meta$index <- as.factor(meta$index)
summary(meta)

PCA_s <- prcomp(t(log2(data_irs_tmm[,-c(10,20)])))
summary(PCA_s)

autoplot(PCA_s, data = meta, 
         frame = TRUE,
         frame.colour = 'group',
         colour = 'group',
         shape = 'probands',
         label = TRUE, 
         label.label = 'index',
         size = 6) +
  theme_classic() +
  scale_fill_manual(name = "Group",
                    values = c('#253494', '#41B6C4', '#EDF8B1', '#253494', '#41B6C4', '#EDF8B1'),
                    breaks=c('Adult_MOCK', 'Adult_EC', 'Adult_EV', 'Term_MOCK', 'Term_EC', 'Term_EV'),
                    labels = c('MOCK', 'EC', 'EV', 'MOCK', 'EC', 'EV')) +
  scale_colour_manual(name = "Group",
                    values = c('#253494', '#41B6C4', '#EDF8B1', '#253494', '#41B6C4', '#EDF8B1'),
                    breaks=c('Adult_MOCK', 'Adult_EC', 'Adult_EV', 'Term_MOCK', 'Term_EC', 'Term_EV'),
                    labels = c('MOCK', 'EC', 'EV', 'MOCK', 'EC', 'EV')) +
  scale_shape_manual(name = "Probands",
                     breaks = c('Adult', 'Term'),
                     values=c(1, 2),
                     labels = c('Adult', 'Term')) +
  guides(colour = FALSE)

######### PCA Adult only #########

data_adult <- data_irs_tmm[,c(1:4, 11:15)]
meta_adult <- meta%>%
  dplyr::filter(probands == "Adult")

PCA_Adult <- prcomp(t(log2(data_irs_tmm[,c(1:4, 11:15)])))
summary(PCA_Adult)

autoplot(PCA_Adult, data = meta_adult, 
         frame = TRUE,
         frame.colour = 'group',
         colour = 'group',
         label = TRUE, 
         label.label = 'index',
         shape = FALSE,
         size = 6
         ) +
  theme_classic() +
  scale_fill_manual(name = "Group",
                    values = c('#253494', '#41B6C4', '#EDF8B1'),
                    breaks=c('Adult_MOCK', 'Adult_EC', 'Adult_EV'),
                    labels = c('MOCK', 'EC', 'EV')) +
  scale_colour_manual(name = "Group",
                      values = c('#253494', '#41B6C4', '#EDF8B1'),
                      breaks=c('Adult_MOCK', 'Adult_EC', 'Adult_EV'),
                      labels = c('MOCK', 'EC', 'EV')) +
  guides(colour = FALSE)



######################
### Final datasets ###
######################

data_prot_nona = data_irs_tmm[,-c(10,20)]
data_prot_nona_log = log2(data_irs_tmm[,-c(10,20)])



########################################
### Differentially abundant proteins ###
########################################

##### Create Eset #####

all(rownames(meta)==colnames(data_prot_nona_log))
# FALSE
meta_new = meta[match(colnames(data_prot_nona_log), rownames(meta)),] 
all(rownames(meta_new)==colnames(data_prot_nona_log))
# TRUE

pset <- ExpressionSet(assayData = as.matrix(data_prot_nona_log), phenoData = AnnotatedDataFrame(meta_new))
dim(pset)

# adjustment for paired samples is implemented within the design matrix
design <- model.matrix(~0 + group + index, meta_new)

colnames(design) <- c("Adult_MOCK", "Adult_EC", "Adult_EV", "Term_MOCK", "Term_EC", "Term_EV",
                      "index2", "index3", "index4", "index5", "index6")
design

##### contrasts ####
contrast <- makeContrasts(contrasts = c("Adult_EC-Adult_MOCK",
                                        "Adult_EV-Adult_MOCK",
                                        "Adult_EV-Adult_EC",
                                        "Term_EC-Term_MOCK",
                                        "Term_EV-Term_MOCK",
                                        "Term_EV-Term_EC",
                                        "Adult_MOCK-Term_MOCK",
                                        "Adult_EV-Term_EV",
                                        "Adult_EC-Term_EC"), 
                          levels = design) 

contrast

######################  limma  ###########################
fit.trend <- lmFit(data_prot_nona_log, design)
fit.trend2 <- eBayes(contrasts.fit(fit.trend, contrast), trend = T)

summary(decideTests(fit.trend2, p.value = 0.05))

coefs <- colnames(contrast) %>%
  set_names(.,.)

result <- coefs %>%
  purrr::map(function(coef, fit){
    topTable(fit = fit,
             coef = coef,
             number = Inf)
  },
  fit = fit.trend2)


########### Differentially abundant proteins (DAPs) of every contrast #############

##### Adult_EC-Adult_MOCK #####
DAP_Adult_EC <- result$`Adult_EC-Adult_MOCK` %>%
  filter(adj.P.Val < 0.05)

# 0 DAPs with p < 0.05  


##### Adult_EV-Adult_MOCK #####
DAP_Adult_EV <- result$`Adult_EV-Adult_MOCK` %>%
    filter(adj.P.Val < 0.05) %>%
    rownames_to_column(var = "UniProt") %>%
    mutate(entrez = mapIds(org.Hs.eg.db, keys=UniProt, column="ENTREZID", keytype="UNIPROT", multiVals="first")) %>%
    mutate(symbol = mapIds(org.Hs.eg.db, keys=UniProt, column="SYMBOL", keytype="UNIPROT", multiVals="first")) %>%
    mutate(genename = mapIds(org.Hs.eg.db, keys=UniProt, column="GENENAME", keytype="UNIPROT", multiVals="first"))
    
sum(is.na(DAP_Adult_EV$entrez))

DAP_Adult_EV <- DAP_Adult_EV %>%
  mutate(entrez = ifelse(UniProt == "Q58FG1", "3323", entrez))

sum(is.na(DAP_Adult_EV$entrez))


##### Adult_EV-Adult_EC #####
DAP_Adult_EV_EC <- result$`Adult_EV-Adult_EC` %>%
  filter(adj.P.Val < 0.05) %>%
  rownames_to_column(var = "UniProt") %>%
  mutate(entrez = mapIds(org.Hs.eg.db, keys=UniProt, column="ENTREZID", keytype="UNIPROT", multiVals="first")) %>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys=UniProt, column="SYMBOL", keytype="UNIPROT", multiVals="first")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys=UniProt, column="GENENAME", keytype="UNIPROT", multiVals="first"))

sum(is.na(DAP_Adult_EV_EC$entrez))

DAP_Adult_EV_EC <- DAP_Adult_EV_EC %>%
  mutate(entrez = ifelse(UniProt == "Q58FG1", "3323", entrez)) %>%
  drop_na(entrez)

sum(is.na(DAP_Adult_EV_EC$entrez))
sum(is.na(DAP_Adult_EV_EC$symbol))


##### Term_EC-Term_MOCK #####
DAP_Term_EC <- result$`Term_EC-Term_MOCK` %>%
  filter(adj.P.Val < 0.05)

# 0 DAPs with p < 0.05


##### Term_EV-Term_MOCK #####
DAP_Term_EV <- result$`Term_EV-Term_MOCK` %>%
  filter(adj.P.Val < 0.05) %>%
  rownames_to_column(var = "UniProt") %>%
  mutate(entrez = mapIds(org.Hs.eg.db, keys=UniProt, column="ENTREZID", keytype="UNIPROT", multiVals="first")) %>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys=UniProt, column="SYMBOL", keytype="UNIPROT", multiVals="first")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys=UniProt, column="GENENAME", keytype="UNIPROT", multiVals="first"))

sum(is.na(DAP_Term_EV$entrez))

DAP_Term_EV <- DAP_Term_EV %>%
  mutate(entrez = ifelse(UniProt == "Q58FG1", "3323", entrez)) %>%
  mutate(entrez = ifelse(UniProt == "Q13670", "107161145", entrez)) %>%
  drop_na(entrez)

sum(is.na(DAP_Term_EV$entrez))
sum(is.na(DAP_Term_EV$symbol))


##### Term_EV-Term_EC #####
DAP_Term_EV_EC <- result$`Term_EV-Term_EC` %>%
  filter(adj.P.Val < 0.05) %>%
  rownames_to_column(var = "UniProt") %>%
  mutate(entrez = mapIds(org.Hs.eg.db, keys=UniProt, column="ENTREZID", keytype="UNIPROT", multiVals="first")) %>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys=UniProt, column="SYMBOL", keytype="UNIPROT", multiVals="first")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys=UniProt, column="GENENAME", keytype="UNIPROT", multiVals="first"))

sum(is.na(DAP_Term_EV_EC$entrez))

DAP_Term_EV_EC <- DAP_Term_EV_EC %>%
  mutate(entrez = ifelse(UniProt == "Q58FG1", "3323", entrez)) %>%
  drop_na(entrez)

sum(is.na(DAP_Term_EV_EC$entrez))
sum(is.na(DAP_Term_EV_EC$symbol))


##### Adult_MOCK-Term_MOCK #####
DAP_MOCK <- result$`Adult_MOCK-Term_MOCK` %>%
  filter(adj.P.Val < 0.05) %>%
  rownames_to_column(var = "UniProt") %>%
  mutate(entrez = mapIds(org.Hs.eg.db, keys=UniProt, column="ENTREZID", keytype="UNIPROT", multiVals="first")) %>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys=UniProt, column="SYMBOL", keytype="UNIPROT", multiVals="first")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys=UniProt, column="GENENAME", keytype="UNIPROT", multiVals="first"))

sum(is.na(DAP_MOCK$entrez))

DAP_MOCK <- DAP_MOCK %>%
  mutate(entrez = ifelse(UniProt == "Q9H853", "80086", entrez)) %>%
  mutate(entrez = ifelse(UniProt == "H7BZ55", "728763", entrez)) %>%
  drop_na(entrez)

sum(is.na(DAP_MOCK$entrez))
sum(is.na(DAP_MOCK$symbol))


##### Adult_EC-Term_EC #####
DAP_EC <- result$`Adult_EC-Term_EC` %>%
  filter(adj.P.Val < 0.05) %>%
  rownames_to_column(var = "UniProt") %>%
  mutate(entrez = mapIds(org.Hs.eg.db, keys=UniProt, column="ENTREZID", keytype="UNIPROT", multiVals="first")) %>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys=UniProt, column="SYMBOL", keytype="UNIPROT", multiVals="first")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys=UniProt, column="GENENAME", keytype="UNIPROT", multiVals="first"))

sum(is.na(DAP_EC$entrez))

DAP_EC <- DAP_EC %>%
  mutate(entrez = ifelse(UniProt == "Q9H853", "80086", entrez)) %>%
  mutate(entrez = ifelse(UniProt == "H7BZ55", "728763", entrez)) %>%
  mutate(entrez = ifelse(UniProt == "P19013", "3851", entrez)) %>%
  drop_na(entrez)

sum(is.na(DAP_EC$entrez))
sum(is.na(DAP_EC$symbol))


##### Adult_EV-Term_EV #####
DAP_EV <- result$`Adult_EV-Term_EV` %>%
  filter(adj.P.Val < 0.05) %>%
  rownames_to_column(var = "UniProt") %>%
  mutate(entrez = mapIds(org.Hs.eg.db, keys=UniProt, column="ENTREZID", keytype="UNIPROT", multiVals="first")) %>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys=UniProt, column="SYMBOL", keytype="UNIPROT", multiVals="first")) %>%
  mutate(genename = mapIds(org.Hs.eg.db, keys=UniProt, column="GENENAME", keytype="UNIPROT", multiVals="first"))

sum(is.na(DAP_EV$entrez))

DAP_EV <- DAP_EV %>%
  mutate(entrez = ifelse(UniProt == "H7BZ55", "728763", entrez)) %>%
  mutate(entrez = ifelse(UniProt == "P19013", "3851", entrez)) %>%
  drop_na(entrez)

sum(is.na(DAP_EV$entrez))
sum(is.na(DAP_EV$symbol))


####################
### Volcano Plot ###
####################

##### Volcano Plot Term_EV-Term_MOCK
DAP_Term_EV_entrez <- DAP_Term_EV$entrez

write.table(DAP_Term_EV_entrez, file = "DAP_Term_EV_entrez.txt", quote = FALSE, row.names = FALSE)

EnhancedVolcano(DAP_Term_EV,
                lab = DAP_Term_EV$symbol,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1)


##### Volcano Plot Adult_EV-Adult_MOCK 
DAP_Adult_EV_entrez <- DAP_Adult_EV$entrez

write.table(DAP_Adult_EV_entrez, file = "DAP_Adult_EV_entrez.txt", quote = FALSE, row.names = FALSE)

EnhancedVolcano(DAP_Adult_EV,
                lab = DAP_Adult_EV$symbol,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1)


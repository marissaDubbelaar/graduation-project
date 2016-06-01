####################################################################
# Author    : M. Dubbelaar
# Date      : 08-april-2016
# File Name : CalculateTPMbtGTF.R
# Purpose   : Obtains the TPM values of the genes within a 
####################################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("biomaRt")
#biocLite("Biostrings")
#install.packages("plyr")
#install.packages("stringr")
library("edgeR")
library("biomaRt")
library(Biostrings)
library(plyr)
library(stringr)
####################################################################
#                              Functions                           #
####################################################################
setRowname <- function(data) {
  # This function checks for the column with the name "probe"
  # This column will be used as the rownames and the column will be
  # put to null (it will remove the colum).
  colNumber <- grep("probes", colnames(data))
  if (colNumber != 0) {
    row.names(data) <- data[, colNumber]
    data[,colNumber] <- NULL
    return (data) 
  }
} 

# The functions below are used to split a string and to obtain a part
# of when that is necessary.
splitOrganism <- function(s) strsplit(s, "G")[[1]][1]
obtainTranscriptInfo <- function(s) strsplit(s, "; ")[[1]]
splitSampleName <- function(s) strsplit(s, "_")[[1]]
####################################################################
# A connection is made with molgenis to obtain the needed information.
source("http://localhost:8080/molgenis.R")
molgenis.login("admin", "admin")
####################################################################
#                            Load Dataset                          #
####################################################################
# A Study id is give to obtain, the count file and target file of that
# given study. It is also used to write the TPM file.
studyID <- "EMTAB1030"
# The count file is obtained from MOLGENIS.
M1 <- molgenis.get(entity = studyID)
# The column probe is used to define the rownames.
M1 <- setRowname(M1)

# The targets file is obtained from MOLGENIS.
targets <- molgenis.get(entity = paste(studyID, "targets", sep = "_"))
# The replicate number and the '_' are removed from the sample name.
targets$Description <- gsub('(_)+([0-9]+)?$', '', targets$Description)
# The DGE list is made.
dge <- DGEList(counts=M1, group = factor(targets$Description))
####################################################################
#                    Create design and normalize                   #
####################################################################
# The data is filtered with a cut-off of 1 count per million (cpm).
# The data with 2 replicated samples will be saves into the vector.
isExpr <- rowSums(cpm(dge)>1) >= 2
# These expressed genes will be saved within the dge and BioM dataset.
dge <- dge[isExpr, ]
# The raw library sizes are calculated with the use of calNormFactors
# This is a step before calculating the estimates.
dge <- calcNormFactors(dge)
data <- matrix(nrow = length(rownames(dge$counts)), ncol = 0)
colName <- NULL

# Obtains the different conditions within the dge list.
for (item in as.character(unique(dge$samples[,1]))) {
  # Creates a shorter name for the condition.
  for (i in paste(splitSampleName(item)[2] , splitSampleName(item)[3], sep="_" )){
    # Checks if there are any wild type conditions
    if (grepl("wt|wildtype|wild type|health|healthy|control|normal", gsub('_|-', '', paste(tolower(as.character(unique(dge$samples[,1]))), collapse = "")))) {
      # Makes sure that the mutation or type is added to the name
      for (i in paste(splitSampleName(item)[2], splitSampleName(item)[3] , splitSampleName(item)[4], sep="_" )){
        # Only the Wildtype samples are used further.
        if (grepl("wt|wildtype|wild type|health|healthy|control|normal", gsub('_|-', '', tolower(i)))){
          # This if-else statement creates the unique name for the condition.
          if (splitSampleName(item)[2] != "") {
            colName <- c(colName, splitSampleName(item)[2])
          } else if (splitSampleName(item)[3] != "") {
            colName <- c(colName, splitSampleName(item)[3])
          } else {
            colName <- c(colName, splitSampleName(item)[1])
          }
          # Data of this sample is added into the data matrix.
          data <- cbind(data, apply(as.matrix(dge$counts[,as.numeric(grep(tolower(i), tolower(dge$samples[,1])))]), 1, mean))
          }
      }
    # The else statement is given when there is no wild type found. 
    } else {
      # This if-else statement creates the unique name for the condition.
      if (splitSampleName(item)[2] != "") {
        colName <- c(colName, splitSampleName(item)[2])
      } else if (splitSampleName(item)[3] != "") {
        colName <- c(colName, splitSampleName(item)[3])
      } else {
        colName <- c(colName, splitSampleName(item)[1])
      }
      # Data of this sample is added into the data matrix.
      data <- cbind(data, apply(as.matrix(dge$counts[,as.numeric(grep(tolower(i), tolower(dge$samples[,1])))]), 1, mean))
      }
    }
}
# The information saved in calName is used as column name for the data matrix.
colnames(data) <- colName
# Replicate coloms are filtered out of the matrix.
data <- t(unique(t(data)))

# Makes sure that the first part of the ensemble name is obtained
# to determine the organism.
ensemblOrganism <- splitOrganism(rownames(data)[1])
####################################################################
#                    Obtain the gene information                   #
####################################################################
# Checks if the genes are from mice, if yes obtain the mice
# ensemble data, otherwise the human ensemble data.
if ("ENSMUS" %in% ensemblOrganism) {
  # The transcript files created from the GTF file are saved.
  transcriptInfo <- read.csv("/Users/marissa/Desktop/goadFiles/Transcript_Files/Mus_musculus_transcripts.txt", sep=";", header=F, col.names = paste0("V",seq_len(20)), fill = T, na.strings=c("","NA"), stringsAsFactors = F)
  # Information as the gene_id and gene_name are removed.
  filteredInfo <- cbind(gsub("gene_id ", "", as.matrix(transcriptInfo[,1])), gsub("gene_name ", "", as.matrix(transcriptInfo[,5])), apply(transcriptInfo, 1, function(x) tail(na.omit(x), 1)))
} else {
  # The transcript files created from the GTF file are saved.
  transcriptInfo <- read.csv("/Users/marissa/Desktop/goadFiles/Transcript_Files/Homo_sapiens_transcripts.txt", sep=";", header=F, col.names = paste0("V",seq_len(12)), fill = T, na.strings=c("","NA"), stringsAsFactors = F)
  # Information as the gene_id and gene_name are removed.
  filteredInfo <- cbind(gsub("gene_id ", "", as.matrix(transcriptInfo[,1])), gsub("gene_name ", "", as.matrix(transcriptInfo[,3])), apply(transcriptInfo, 1, function(x) tail(na.omit(x), 1)))
}

# The min, mean and max transcript lengths are determined.
minTranscript <- as.matrix(sapply(split(as.integer(filteredInfo[,3]), filteredInfo[,1]), min))
meanTranscript <- as.matrix(sapply(split(as.integer(filteredInfo[,3]), filteredInfo[,1]), mean))
maxTranscript <- as.matrix(sapply(split(as.integer(filteredInfo[,3]), filteredInfo[,1]), max))

# The information from the transcript file is marged with the mean, min and max 
# transcript information.
filteredInfo <- filteredInfo[match(rownames(meanTranscript), filteredInfo[,1]),]
filteredInfo <- cbind(filteredInfo[,1:2], meanTranscript, minTranscript, maxTranscript)
# Rows that start with ENS (that contain a ensemble annotation) are kept.
rowsToKeep <- which(grepl("^ENS", rownames(filteredInfo)))
filteredInfo <- unname(filteredInfo[rowsToKeep,])

filteredInfo <- as.data.frame(filteredInfo[match(rownames(as.matrix(data)), filteredInfo[,1]),], stringsAsFactors=F)
# Data is written as integer (to solve difficulties when calculating the TPM values).
filteredInfo[,3] <- as.integer(filteredInfo[,3])
filteredInfo[,4] <- as.integer(filteredInfo[,4])
filteredInfo[,5] <- as.integer(filteredInfo[,5])
# The colom names are added to the data frame filteredInfo
colnames(filteredInfo) <- c("ensemble_id", "external_gene_name", "transcript_length", "transcript_length_min", "transcript_length_max")
####################################################################
#                       Calculate the TPM value                    #
####################################################################
# Source = https://github.com/johnstantongeddes/sim-transcriptome/blob/master/TPM.R
####################################################################
# TPM = (Rg * 10^6) / (Tn * Lg)
# Rg = number of reads mapped to a particular transcript g = count
# Tn = sum of all length normalized transcript counts
# Lg = length of transcript g (kilobases)
tpmMatrix <- matrix(nrow = length(filteredInfo[,1]), ncol = 0)
tpmMatrix <- cbind(tpmMatrix, as.data.frame(filteredInfo[,2]))
colName <- NULL

for (item in colnames(data)){
  for (i in 1:3) {
    # Obtains the counts of the matrix data.
    count <- as.matrix(data[,which(as.matrix(colnames(data) %in% item))])
    tempMat <- cbind(filteredInfo[,c(1,2,2+i)], count)
    if (i == 2) {
      # The tpm low is calculated, and the tpm low name is saved into colName
      Tn <- sum(ddply(tempMat, .(external_gene_name), summarize, Tn = count/(transcript_length_min))[2], na.rm = T)
      tpmMatrix <- cbind(tpmMatrix, ddply(tempMat, .(external_gene_name), summarize, tpm_low = (count*1e6)/(Tn*transcript_length_min))[2])
      colName <- c(colName, paste(item, "low", sep="_"))
    } else if (i == 3) {
      # The tpm high is calculated, and the tpm high name is saved into colName
      Tn <- sum(ddply(tempMat, .(external_gene_name), summarize, Tn = count/(transcript_length_max))[2], na.rm = T)
      tpmMatrix <- cbind(tpmMatrix, ddply(tempMat, .(external_gene_name), summarize, tpm_high = (count*1e6)/(Tn*transcript_length_max))[2])
      colName <- c(colName, paste(item, "high", sep="_"))
    } else {
      # The tpm is calculated, and the tpm name is saved into colName
      Tn <- sum(ddply(tempMat, .(external_gene_name), summarize, Tn = count/(transcript_length))[2], na.rm = T)
      tpmMatrix <- cbind(tpmMatrix, ddply(tempMat, .(external_gene_name), summarize, tpm = (count*1e6)/(Tn*transcript_length))[2])
      colName <- c(colName, item)
    }
  }
}

# The external gene name is given as column name for the gene symbols and the 
# found information (saved in colName) is used to identify the other columns.
colnames(tpmMatrix) <- c("external_gene_name", colName)
# The '-' in colum names are replaced by '_'.
colnames(tpmMatrix) <- gsub('-', '_', colnames(tpmMatrix))[1:length(colnames(tpmMatrix))]

# Abbreviations and other difficulties with the colom names can be added below.
colnames(tpmMatrix) <- gsub("oligodendrocyte_precursor_cell", "OPC", colnames(tpmMatrix))
colnames(tpmMatrix) <- gsub("newly_formed_oligodendrocyte", "NFO", colnames(tpmMatrix))
colnames(tpmMatrix) <- gsub("myelinating_oligodendrocyte", "MO", colnames(tpmMatrix))
colnames(tpmMatrix) <- gsub("Acutely_isolated_microglia", "microglia", colnames(tpmMatrix))
colnames(tpmMatrix) <- gsub("Acutely_isolated_microglia", "microglia", colnames(tpmMatrix))
colnames(tpmMatrix) <- gsub("superior_temporal_gyrus", "STG", colnames(tpmMatrix))

# A synonym of the PISD gene is used since mice data contains two PISD genes (one uper case and
# one camelcase).
tpmMatrix$external_gene_name <- gsub("PISD", "PSDC", tpmMatrix$external_gene_name)
tpmMatrix <- tpmMatrix[!duplicated(tpmMatrix$external_gene_name),]

# This for loop is used since the tpm low and high values swap some of the time.
for (column in 2:length(colnames(tpmMatrix))) {
  if (column%%3 == 0) {
    # Making sure that the tpm low contains the low values
    tpmLow <- as.matrix(apply(tpmMatrix[c(column,column+1)], 1, min))
    # And the tpm high, contains the high values.
    tpmHigh <- as.matrix(apply(tpmMatrix[c(column,column+1)], 1, max))
    # This information is rewritten within the tpmMatrix
    tpmMatrix[c(column,column+1)] <- cbind(tpmLow, tpmHigh)
  }
}

# Removes rows with NA's from the tpmMatrix
tpmMatrix <- tpmMatrix[complete.cases(tpmMatrix),]
# Removes whitespaces within the columns of the tpmMatrix
tpmMatrix <- apply(tpmMatrix,2,function(x)gsub('\\s+', '',x))
tpmMatrix <- as.data.frame(tpmMatrix, stringsAsFactors = F)
#################################################################
#                      Obtaining Percentiles                    #
#################################################################
tpmMatrix[,2:length(colnames(tpmMatrix))] <- lapply(tpmMatrix[,2:length(colnames(tpmMatrix))], as.numeric)
Expressed <- tpmMatrix[,1] > 0 

count = 2
for (i in 1:(length(colnames(tpmMatrix)) + ((length(colnames(tpmMatrix))-1)/3))) {
   # if the column is 2 after %% 4 (also known as the table with the normal TPM values).
   if (i %% 4 == 2) {
      # The percentile column is added into the data.frame (only containing columns with the text "TEXT")
      tpmMatrix <- as.data.frame(append(tpmMatrix, list(percentile=rep("TEXT", length(rownames(tpmMatrix)))), match(as.character(colnames(tpmMatrix)[i+2]), names(tpmMatrix))),  stringsAsFactors = F)
      # The column name is renamed after the cell type with "_percentile" behind it.
      colnames(tpmMatrix)[i+3] <- paste(colnames(tpmMatrix[i]), "percentile", sep="_")
      
      # The significance of the genes is determined and a percentile indication is given.
      tpmMatrix[as.numeric(tpmMatrix[,i+1]) == 0,i+3] = c("Not significantly expressed")
      tpmMatrix[as.numeric(tpmMatrix[,i]) == 0,i+3] = c("Not expressed")
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0)] <- "90-100 percentile: very low expression"
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0.1)] <- "80-90 percentile: low expression"
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0.2)] <- "70-80 percentile: low expression"
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0.3)] <- "60-70 percentile: low expression"
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0.4)] <- "50-60 percentile: moderate expression"
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0.5)] <- "40-50 percentile: moderate expression"
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0.6)] <- "30-40 percentile: moderate expression"
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0.7)] <- "20-30 percentile: moderately high expression"
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0.8)] <- "10-20 percentile: moderately high expression"
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0.9)] <- "5-10 percentile: high expression"
      tpmMatrix[Expressed, i+3][tpmMatrix[Expressed, i] > quantile(as.numeric(tpmMatrix[Expressed, i]), probs=0.95)] <- "0-5 percentile: very high expression"
      count <- count + 3
    }
}

# Information from the tpmMatrix is written into a CSV file.
write.table(tpmMatrix, file = paste("/Users/marissa/Desktop/goadFiles/TPM/TPM_", studyID, ".csv", sep = ""), eol="\n", quote = F, row.names = F, sep=",")

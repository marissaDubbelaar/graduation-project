##############################################################
# Author    : M. Dubbelaar
# Date      : 08-april-2016
# File Name : databaseRAPI.R
# Purpose   : Performs the normalization, the pairwise
#             DE analysis and the creation of the matrix
##############################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library("edgeR")
library("RCurl")
##############################################################
#                          Functions                         #
##############################################################
setRowname <- function(data) {
  # This function checks for the column with the name "probe"
  # This column will be used as the rownames and the column will be
  # put to null (it will remove the colum).
  colNumber <- grep("probes", colnames(data))
  row.names(data) <- data[, colNumber]
  data[,colNumber] <- NULL
  return (data)
} 

obtainCondition <- function(conditions, dataTab) {
  # This function is used to add information from one sample together
  # It makes sure that these conditions can be used to specify the 
  # group to create the model.matrix for the DE analysis.
  # A is returned containing each condition found in the dge dataset.
  count <- 1
  dataName <- matrix(NA, nrow = nrow(t(dataTab)), ncol = 1)
  for (x in 1:conditions) {
    dataName[count] <- apply(as.matrix(dataTab[,x][4:length(dataTab[,x])-2]), 2, paste, collapse= " ")
    count <- count + 1
  } 
  return (dataName)
}

createDGEGroups <- function(conditions, dataTab) {
  # This function is used to add information from one sample together
  # It makes sure that these conditions can be used to specify the 
  # group to create the model.matrix for the DE analysis.
  # A is returned containing each condition found in the dge dataset.
  count <- 1
  dataName <- matrix(NA, nrow = ncol(t(dataTab)), ncol = 1)
  for (x in 1:conditions) {
    dataName[count] <- dataTab[x]
    count <- count + 1
  }
  return(dataName)
}

# The function separate is used to split the vector on "_"
separate <- function(s) strsplit(s, "_")[[1]]
# The function findGroup is used to split the vector on ")"
# and is used to remove the factor in front of the unique condition.
findGroup <- function(s) strsplit(s, ")")[[1]][2]
# Splits the ensemble name on the 'G' making sure that the organism
# can be validated.
splitOrganism <- function(s) strsplit(s, "G")[[1]][1]
##############################################################
# A connection is made with molgenis to obtain the needed information.
source("http://localhost:8080/molgenis.R")
molgenis.login("admin", "admin")
conditie1 = "${condition1}"
conditie2 = "${condition2}"  
##############################################################
#                       Load Dataset                         #
##############################################################
# Reads the count file
M1 <- molgenis.get(entity = "${entityName}")

# The column probe is used to define the rownames.
M1 <- setRowname(M1)

# Reads the targets file.
targets <- molgenis.get(entity = "${targetFile}")

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

# The design and all of the needed calculations are made.
design <- model.matrix(~0+factor(targets$Description), data = dge$samples )
# lmfit fits the linear model.
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
##############################################################
#            Check which condition is specified              #
#                   and create the contrast                  #
##############################################################
# Returns the names of the chosen conditions by the user
######################
con1index <- which(targets$SRA==conditie1)
con2index <- which(targets$SRA==conditie2)

# Makes sure that the format of the colnames within the design matrix 
# are the same as the format of the given conditions.
# After this the index of the design with the given condition is obtained.
simplefiedDesignColName <- sapply(colnames(design), findGroup)
indexCon1 <- which(simplefiedDesignColName %in% targets$Description[con1index])
indexCon2 <- which(simplefiedDesignColName %in% targets$Description[con2index])

# Given the indexes obtained above, the contrast matrix is filled in.
# Where condition 1 is compared against condition 2.
madeContrast <- rep(0, length(colnames(design)))
######################
# Adjust the strength of the contrasts
######################
madeContrast[indexCon1] = 1
madeContrast[indexCon2] = -1
##############################################################
#                Obtain the significant DEG                  #
##############################################################
# The DEG are calculated and are filtered on a significant FDR
# value of 0.05 (Benjamini Corrected values).
lrt <- glmLRT(fit, contrast = madeContrast)
toptable <- topTags(lrt, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
uniqueGenes <- which(toptable[[1]]$FDR < 0.05)
toptable <- toptable[uniqueGenes,]

if (length(toptable) != 0) {
  ##############################################################
  #                    Obtain genesymbols                      #
  ##############################################################
  # A connection is made with the ensembl database to retreive the 
  # genesymbols that match with the already known ensemble names.
  ######################
  # Checks if the genes are from mice, if yes obtain the mice
  # ensemble data, otherwise the human ensemble data.
  if ("${organism}" == "Mus musculus") {
    geneList <- molgenis.get(entity = "miceGenes")
  } else if ("${organism}" == "Homo sapiens") {
    geneList <- molgenis.get(entity = "humanGenes")
  }
  
  # The getBM function retrieves the given attributes. 
  genes <- geneList[match(rownames(toptable), geneList[,2]), ]
  # The FDR values are rounded and the genesymbols are added into the matrix.
  toptable <- cbind(toptable, genes[1])
  toptable[6] <- lapply(toptable[6], as.character)
  # NA's within the genesymbols are captured and replaced by there ensemble name.
  toptable$Associated_Gene_Name[is.na(toptable$Associated_Gene_Name)] <- as.character(rownames(toptable)[is.na(toptable$Associated_Gene_Name)])
  toptable[,5] <- round(toptable[,5], 9)
  toptable <- toptable[!duplicated(toptable[,6]),]
  toptable <- toptable[!grepl("^20", rownames(toptable)),]
  rownames(toptable) <- toptable[,6]
  print(as.matrix(toptable[c(1,5)]))
} else {
  print(paste("No differentially expressed genes where found with comparison:", sub("^\\s+|\\s+$|  ", "", gsub("_", " ", targets$Description[which(targets$SRA==conditie1)])), "vs.", sub("\\s+$", "", gsub("_", " ", targets$Description[which(targets$SRA==conditie1)]))))
}

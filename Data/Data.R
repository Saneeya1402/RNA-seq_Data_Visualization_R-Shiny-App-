###>>>>> DATA GENERATION FOR R FINAL PROJECT BF591 <<<<<###


########>>>>>>>  1.    METADATA

# Install GEOquery if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

library(tidyverse)
library(GEOquery)

# Load GEO data
geo_data <- getGEO("GSE64810", GSEMatrix = TRUE)
# Check if geo_data is a list and extract the first element (first platform)
if (length(geo_data) > 1) {
  cat("Multiple platforms detected. Selecting the first platform.\n")
}
# Extract metadata from the first platform
metadata <- pData(geo_data[[1]])
# View the metadata
head(metadata)

# Loading counts data
counts <- read.delim("Counts data.csv", stringsAsFactors = FALSE)

# Selecting specific colus of the metadata
mod.metadata <- metadata%>%select(characteristics_ch1.1:characteristics_ch1.11)
mod.metadata[mod.metadata==""] <- NA
# Save the cleaned data as a CSV file
write.csv(mod.metadata, "mod_metadata_cleaned.csv", row.names = FALSE)






#######>>>>>   2.    FGSEA DATA


# Load necessary libraries
library(tidyverse)

# Load DESeq results
deseq_results <- read.csv("DESeq_results.csv")

# Remove rows with NA in log2FoldChange and arrange by descending log2FoldChange
modified_deseq <- deseq_results %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

# Save the modified DESeq results
write.csv(modified_deseq, "modified_deseq_results.csv", row.names = FALSE)


library(fgsea)
library(dplyr)
library(data.table)

# Read the ordered DESeq results
deseq_results <- read.csv("modified_deseq_results.csv")

# Read the GMT file
pathways <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt")

# Prepare ranking vector (using log2FoldChange)
ranks <- deseq_results$log2FoldChange
names(ranks) <- deseq_results$Gene.symbol  # Ensure column name matches your file

# Perform FGSEA
set.seed(123)  # for reproducibility
fgsea_results <- fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 15,   # minimum gene set size
  maxSize = 500  # maximum gene set size
)

# Add more information and sort
fgsea_results_detailed <- fgsea_results %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = ","),
    pathway = gsub("HALLMARK_", "", pathway)  # Optional: clean pathway names
  )

# Write detailed results
write.csv(fgsea_results_detailed, "fgsea_results.csv", row.names = FALSE)



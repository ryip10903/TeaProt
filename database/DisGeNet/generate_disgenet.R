# Script to edit DisGeNet
# First dowload All gene disease associations from https://www.disgenet.org/downloads
# the file is called: all_gene_disease_associations.tsv.gz
library(tidyverse)

# Load the downloaded datafile
# Data is filtered using the significance cutoff as defined by the IMPC (1e-4)
# https://www.mousephenotype.org/understand/start-using-the-impc/how-to-use-gene-pages/

df <- readr::read_delim("curated_gene_disease_associations.tsv.gz") %>% select(1,2,6)

df %>% select(geneSymbol, diseaseName) %>% distinct() %>% group_by(geneSymbol) %>% 
  summarise(diseaseName = paste(diseaseName, collapse = "|")) %>% 
  write.csv(., file = "disgenet_genesymbol.csv", row.names = FALSE)

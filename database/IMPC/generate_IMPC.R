# Script to download IMPC genotype-phenotype associations
# IMPC data is filtered for significance and summarized for TeaProt

library(tidyverse)

# Download the lastest version of IMPC data
file = "./genotype-phenotype-assertions-IMPC.csv.gz"
url <- "http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-15.1/results/genotype-phenotype-assertions-IMPC.csv.gz"

curl::curl_download(url, file)

# Load the downloaded datafile
# Data is filtered using the significance cutoff as defined by the IMPC (1e-4)
# https://www.mousephenotype.org/understand/start-using-the-impc/how-to-use-gene-pages/
df <- readr::read_csv("genotype-phenotype-assertions-IMPC.csv.gz") %>% filter(p_value < 1e-4)
mapping <- as.data.frame(org.Mm.eg.db::org.Mm.egALIAS2EG)

df %>% select(marker_symbol, procedure_name) %>% distinct() %>% group_by(marker_symbol) %>% 
  summarise(impc_significant_procedure_name = paste(procedure_name, collapse = "|")) %>% left_join(., mapping, by = c("marker_symbol" =  "alias_symbol")) %>% 
  write.csv(., file = "impc_procedure.csv", row.names = FALSE)

df %>% select(marker_symbol, parameter_name) %>% distinct() %>% group_by(marker_symbol) %>% 
  summarise(impc_significant_parameter_name = paste(parameter_name, collapse = "|")) %>% left_join(., mapping, by = c("marker_symbol" =  "alias_symbol")) %>% 
  write.csv(., file = "impc_parameter.csv", row.names = FALSE)


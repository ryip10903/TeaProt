# Download BRENDA from https://www.brenda-enzymes.org/ and place the file in the BRENDA folder
library(tidyverse)

# Load the data and extract chunks of data for each enzymatic reaction
lines <- readr::read_lines(file = "brenda_download.tar.gz", skip = 16, skip_empty_rows = FALSE)

lines_start <- grep("///", lines)
lines_ID <- lines[grep("ID\\\t", lines)]

lines_range <- data.frame(start = lines_start, end = c(lines_start[2:length(lines_start)]-1, length(lines)))
lines_range <- lines_range[1:nrow(lines_range)-1,]

# Make a list containing all Uniprot accessions for each enzymatic reaction
lines_list <- list()

for(i in 1:(nrow(lines_range))){

  lines_list[[i]] <- lines[lines_range[i,1]:lines_range[i,2]]
  lines_list[[i]] <- lines_list[[i]][lines_list[[i]] %>% grep("UniProt|Swissprot", .)] %>% str_extract(., "[A-z,0-9]* (UniProt|Swissprot)") %>% sub(" (UniProt|Swissprot)", "", .)
  names(lines_list)[i] <- lines_ID[i]

}

# Remove empty enzymatic reactions (no uniprot)
lines_list <- lines_list[lapply(lines_list, length) %>% unlist != 0]

# Retain uniprot IDs for each reaction
for(i in 1:length(lines_list)){
  
  lines_list[[i]] <- data.frame(uniprot = lines_list[i], reaction = names(lines_list)[i]) %>% `colnames<-`(c("uniprot", "reaction"))
  
}

# Combine reactions into a table
df <- do.call(rbind, lines_list) %>% group_by(uniprot) %>% summarize(reaction = paste(reaction, collapse = "|"))
df <-  df[9:nrow(df),]

mapping <- read.delim(file = "mapping.tab.txt")

# Add entrez mapping and filter
df <- left_join(df, mapping, by = c("uniprot" = "From")) %>% filter(!is.na(To)) %>% `colnames<-`(c("uniprot", "reaction", "entrez"))

write.csv(df, file = "brenda_processed.csv")

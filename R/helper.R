#' Test whether ID column contains UNIPROT IDs, ENSEMBL IDs or neither
#'
#' \code{cp_idtype} Returns a string indicating the ID type
#'
#' @param id_col A vector containing gene/transcript/protein identifiers
#' 
#' @return A character string indicating the ID type ("uniprot", "ensembl", "neither")
#'
#' @examples
#' cp_idtype(id_colcol)
#' 
cp_idtype <- function(id_col){
  
  if((grepl("^[A-z][0-9][0-9,A-z][0-9,A-z][0-9,A-z][0-9]", id_col) %>% sum / length(id_col) * 100) > 50){
    
    idtype <- "uniprot"
    
  } else if((grepl("^ENS", id_col) %>% sum / length(id_col) * 100) > 50) {
    
    idtype <- "ensembl"
    
  } else {
    
    idtype <- "neither"
    
  }
  
  return(idtype)
  
}


#' Converts the ID column of a data frame to gene names
#'
#' \code{cp_idconvert} Converts the identifiers in the supplied data frame to gene names
#'
#' @param data A vector containing gene/transcript/protein identifiers in a column named "ID"
#' @param id_type A string containing the ID type ("uniprot", "ensembl", "neither")
#' 
#' @return A character string indicating the ID type ("uniprot", "ensembl", "neither")
#'
#' @examples
#' cp_idconvert(df, "uniprot")
#' 
cp_idconvert <- function(data, id_type){
  
  if(id_type == "uniprot"){
    
    uniKeys <- (AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype="SYMBOL")) %>%  c(., AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype="SYMBOL")) #Take all gene symbols from DB 
    Hs_g <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=uniKeys, columns="UNIPROT", keytype="SYMBOL") %>% bind_rows(., AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=uniKeys, columns="UNIPROT", keytype="SYMBOL"))
    
    df <- left_join(data, Hs_g, by = c("ID" = "UNIPROT")) %>% select(-ID) %>% select(1, SYMBOL, everything()) %>% rename(., ID = SYMBOL) %>% filter(is.na(ID) == FALSE)
    
    message("Status: Converted Uniprot to Gene names")
    
  } else if(id_type == "ensembl") {
    
    uniKeys <- (AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype="SYMBOL")) %>%  c(., AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype="SYMBOL")) #Take all gene symbols from DB 
    Hs_g <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=uniKeys, columns="ENSEMBL", keytype="SYMBOL") %>% bind_rows(., AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=uniKeys, columns="ENSEMBL", keytype="SYMBOL"))
    
    df <- left_join(data, Hs_g, by = c("ID" = "ENSEMBL")) %>% select(-ID) %>% select(1, SYMBOL, everything()) %>% rename(., ID = SYMBOL) %>% filter(is.na(ID) == FALSE)
    
    message("Status: Converted Ensembl Gene to Gene names")
    
  } else {
    
    df <- data
    
    message("Status: No gene name conversion needed")
    
  }
  
  return(df)
  
}

cp_dl_table_csv <- function(table, filename){
  downloadHandler(
    filename = function() { filename },
    content = function(file) {write.csv(table, file, row.names = FALSE)}
  )
}


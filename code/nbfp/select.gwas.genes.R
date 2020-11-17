#Read in GWAS gene lists from MAGMA and extract a set of significant genes

select.gwas.genes <- function(file_path = NULL, pval_thresh = 0.01) {
  
  require(tidyverse)
  #read in the magma gene conversion file
  #Make sure the GENE column is read in as character
  #P column should also be read in as a double
  mag_table <- readr::read_table(file_path, 
                                 col_names = TRUE, 
                                 col_types = cols(GENE = col_character(), 
                                                  CHR = col_double(), 
                                                  START = col_double(), 
                                                  STOP = col_double(), 
                                                  NSNPS = col_double(), 
                                                  NPARAM = col_double(), 
                                                  N = col_double(), 
                                                  ZSTAT = col_double(), 
                                                  P = col_double()))
  
  sig_genes <- mag_table %>% 
    dplyr::filter(P < pval_thresh) %>% 
    dplyr::select(GENE)
  sig_genes <- as.character(t(sig_genes))
  
  return(sig_genes)
}
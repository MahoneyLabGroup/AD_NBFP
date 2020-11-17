get_biomart_gene_position <- function(score_tab = NULL,
                                      genomic_pos_start_col = "start",
                                      genomic_pos_stop_col = "stop",
                                      chromosome_col = "chr",
                                      mart = "hsapiens") {
  
  gene_position_biomart <- getBM(attributes = c("external_gene_name","entrezgene_id", "start_position", "end_position", "chromosome_name"), 
                                 filters = "entrezgene_id",
                                 values = score_tab$gene,
                                 mart = {{ hsapiens }}) %>% 
    mutate(entrezgene_id = as.character(entrezgene_id))
  
  score_tab_updated_positions <- score_tab %>% 
    inner_join(gene_position_biomart, by = c("gene" = "entrezgene_id")) %>% 
    dplyr::filter(!duplicated(external_gene_name), !duplicated(gene), !duplicated(name)) %>% 
    dplyr::select(-external_gene_name, -{{ genomic_pos_start_col }}, -{{ genomic_pos_stop_col }}, -{{ chromosome_col }})
  
  score_tab_updated_positions %<>% mutate(mean_gene_position = rowMeans(score_tab_updated_positions[c("start_position", "end_position")]))
  
  return(score_tab_updated_positions)
}
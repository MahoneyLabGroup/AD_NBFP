go_comparison_plotr <- function(score_tab = NULL, gene_col = "gene", pval_col = "log_p_val", score_comp_col = "logFDR", num_nets = "4") {
  
  score_tab_named <- score_tab %>% 
    rename(gene = {{ gene_col }}, log_pval = {{ pval_col }}, comp_col = {{ score_comp_col}}) %>% 
    mutate(across(gene, as.character)) %>% 
    rowwise() %>% 
    mutate(comb_score = sum(log_pval > .$log_pval & comp_col > .$comp_col)) %>%
    ungroup() %>%
    mutate(comb_score = comb_score/n()) 
    
  score_tab_named %<>% arrange(desc(log_pval))
  
  score_tab_pval_gost <- gost(score_tab_named$gene,
                                    organism = "hsapiens",
                                    ordered_query = TRUE,
                                    numeric_ns = "ENTREZGENE_ACC",
                                    sources = c("GO:BP", "REAC", "KEGG")) 
  
  score_tab_pval_gost_res <- score_tab_pval_gost$result %>% 
    as_tibble() %>% 
    select(term_name, pval_score_pval = p_value) %>% 
    mutate(log_pval_score_pval = -log10(pval_score_pval))
  
  score_tab_named %<>% arrange(desc(comp_col))
  
  score_tab_comp_col_gost <- gost(score_tab_named$gene,
                                  organism = "hsapiens",
                                  ordered_query = TRUE,
                                  numeric_ns = "ENTREZGENE_ACC",
                                  sources = c("GO:BP", "REAC", "KEGG")) 
  
  score_tab_comp_col_gost_res <- score_tab_comp_col_gost$result %>% 
    as_tibble() %>% 
    select(term_name, comp_col_pval = p_value) %>% 
    mutate(log_comp_col_pval = -log10(comp_col_pval))
  
  score_tab_named %<>% arrange(desc(comb_score))
  
  score_tab_comb_score_gost <- gost(score_tab_named$gene,
                                    organism = "hsapiens",
                                    ordered_query = TRUE,
                                    numeric_ns = "ENTREZGENE_ACC",
                                    sources = c("GO:BP", "REAC", "KEGG")) 
  
  score_tab_comb_score_gost_res <- score_tab_comb_score_gost$result %>% 
    as_tibble() %>% 
    select(term_name, comb_score_pval = p_value) %>% 
    mutate(log_comb_score_pval = -log10(comb_score_pval))
    
  pval_comp_col_comparison <- score_tab_pval_gost_res %>% 
    full_join(score_tab_comp_col_gost_res, by = "term_name") %>% 
    select(term_name, log_pval_score_pval, log_comp_col_pval) %>% 
    mutate(across(starts_with("log_"), replace_na, 0))
    
  
  pval_comb_score_comparison <- score_tab_pval_gost_res %>% 
    full_join(score_tab_comb_score_gost_res, by = "term_name") %>% 
    select(term_name, log_pval_score_pval, log_comb_score_pval)
  
  x <- list(
    title = "-log 10 P-value GWAS Pval GSEA"
    )
  y1 <- list(
    title = paste("-log 10 P-value of", score_comp_col, "GSEA")
  )
  
  y2 <- list(
    title = "-log10 P-value of Combined Score GSEA"
  )
  
  pval_comp_col_comparison_plot <- plot_ly(pval_comp_col_comparison,
                                           x = ~log_pval_score_pval,
                                           y = ~log_comp_col_pval,
                                           text = ~term_name,
                                           hoverinfo = 'text',
                                           type = "scatter") %>% 
    layout(xaxis = x, yaxis = y1, showlegend = FALSE)
    
  pval_comb_score_comparison_plot <- plot_ly(pval_comb_score_comparison,
                                             x = ~log_pval_score_pval,
                                             y = ~log_comb_score_pval,
                                             text = ~term_name,
                                             hoverinfo = 'text',
                                             type = "scatter") %>% 
    layout(xaxis = x, yaxis = y2, showlegend = FALSE)
  
  side_by_side <- subplot(pval_comp_col_comparison_plot, pval_comb_score_comparison_plot, 
                          titleX = TRUE, 
                          titleY = TRUE) %>% 
    layout(title = paste("Functional Enrichment Comparisions of SIGNET Results with", num_nets,"Brain Networks"))
  return(side_by_side)
}



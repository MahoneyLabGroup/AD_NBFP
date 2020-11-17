gene_rank_plotr <- function(score_tab = NULL,
                            sig_gene_tab = NULL,
                            gene_col = "gene", 
                            log_pval_col = "log_p_val", 
                            log_score_comp_col = "logFDR",
                            set_ceiling = FALSE
) {
  
  score_tab_renamed <- score_tab %>% 
    rename(gene = {{ gene_col }}, log_pval = {{ log_pval_col }}, log_comp_col = {{ log_score_comp_col}}) %>% 
    rowwise() %>%
    mutate(comb_score = sum(log_pval > .$log_pval & log_comp_col > .$log_comp_col)) %>%
    ungroup() %>%
    mutate(comb_score = comb_score/n()) %>% 
    mutate(p_ceiling = ifelse(log_comp_col > 20, 20, log_comp_col))
  
  top_10_comb_score <- score_tab_renamed %>% 
    slice_max(comb_score, n = 10)
  
  if(set_ceiling){
    score_tab_plot <-
      ggplot(score_tab_renamed,
             aes(x = log_pval, y = p_ceiling, label = gene)) +
      geom_point(aes(color = comb_score), size = 5) +
      geom_hline(yintercept = 19, linetype = "dashed", color  = "lightgrey") +
      geom_text(aes(x = 4, y = 18, label = "Genes with -log10 P > 50 transformed to 20 for redability", size = 1), color = "grey") +
      theme(legend.position = "none")
  }else{
    score_tab_plot <-
      ggplot(score_tab_renamed,
             aes(x = log_pval, y = log_comp_col, label = gene)) +
      geom_point(aes(color = comb_score), size = 5) +
      scale_color_viridis()
  }
  
  
  
  if(!(is.null(sig_gene_tab))){
    sig_gene_names <- score_tab_renamed %>% 
      filter(gene %in% sig_gene_tab$name)
    score_tab_plot <- score_tab_plot + geom_point(data = sig_gene_names, color = "lightblue2") 
  }
  
  score_tab_plot <- score_tab_plot + geom_text_repel(data = top_10_comb_score, aes(label = gene, fontface = "italic"), 
                                                     size = 17,
                                                     box.padding = 1,
                                                     #nudge_x = 4, for amygdala
                                                     nudge_x = 5, #for hippocampus
                                                     direction = "y",
                                                     xlim = c(2, Inf))
  return(score_tab_plot)
}
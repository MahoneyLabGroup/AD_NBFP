gene_rank_plotlyr <- function(score_tab = NULL,
                              gene_col = "gene_name",
                              log_p_col = "log_p",
                              log_comp_col = "log_fpr",
                              add_col_to_hover = NULL,
                              title = "Intractive Ranking Plot"
                              ) {
  
  score_tab_renamed <- score_tab %>% 
    rename(gene = {{ gene_col }}, 
           log_pval = {{ log_p_col }}, 
           log_comp_score = {{ log_comp_col}},
           for_hover = {{ add_col_to_hover }}) %>% 
    rowwise() %>%
    mutate(comb_score = sum(log_pval > .$log_pval & log_comp_score > .$log_comp_score)) %>%
    ungroup() %>%
    mutate(comb_score = comb_score/n())
  
  if(!is.null(add_col_to_hover)) {
    score_tab_renamed %<>% 
      rowwise() %>% 
      mutate(hover = paste(gene, for_hover, sep = " | ")) %>% 
      ungroup()
  }else{
    score_tab_renamed %<>% mutate(hover = gene)
  }
  
  x <- list(
    title = "-log10 GWAS P-Value"
  )
  
  y <- list(
    title = "-log10 FPR"
  )
  
  plotly_plot <- plot_ly(score_tab_renamed,
                         x = ~log_pval,
                         y = ~log_comp_score,
                         text = ~hover,
                         color = ~comb_score,
                         hoverinfo = 'text',
                         type = "scatter") %>%
    plotly::layout(xaxis = x,
                   yaxis = y,
                   legend = list(title = list(text ='<b> Combined Score </b>')),
                   title = title
    )
  return(plotly_plot)
}
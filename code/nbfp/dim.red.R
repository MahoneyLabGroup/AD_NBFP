dim.red <- function(net = NULL, rank_opt = 200, gene_fgn_file = NULL) {
  t.bmp <- t(net)
  scale.bmp <- scale(t.bmp)
  
  svd.bmp <- svd(t.bmp, nu = rank_opt)$u
  rownames(svd.bmp) <- rownames(t.bmp)
  
  svd.bmp <- t(svd.bmp)
  
  #Save dimension reduced object for use with generate.triage.models
  saveRDS(svd.bmp, gene_fgn_file)
  
  return(svd.bmp)
}
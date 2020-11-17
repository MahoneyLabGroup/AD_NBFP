make.adj.mat <- function(net_path = NULL, tissue = NULL, sig_genes = NULL, remove_unconnected = FALSE){
  net <- readRDS(net_path)
  
  adj <- tissue.adj.mat(net,
                        gene.list = sig_genes,
                        remove.unconnected = remove_unconnected)
  
  saveRDS(adj, paste(net_path, tissue, "_adj.rds"))
  
  return(adj) 
}
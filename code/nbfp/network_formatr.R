network_formatr <- function(network.dir = NULL, net_extension = c("_top.txt"), edge_threshold = 0.2) {
  nets <- list.files(here("data"), pattern = "_top.txt")
  
  for(i in 1:length(brain_nets)){
    net <- read_delim(here("data", nets[i]), delim = "\t", col_names = FALSE) #not great practice but okay for now.
    
    net_filt <- net %>% 
      filter(X3 > edge_threshold)
    
    write_delim(net_filt, here("hb_brain_nets", paste("filt", nets[i], sep = "_")), delim = "\t", col_names = FALSE) #should be able to specify output dir
    
    rm(net)
  }
}
multi.net.fusion <- function(adj_pattern = c("_adj", "_adj_hv", "_adj_dx")){
  
  #Get list of adjacency matrices from the global environment
  adjs <- ls(.GlobalEnv, pattern = adj_pattern)
  
  #Check the lengths of each adj mat to find the smallest
  length_check <- matrix(0, nrow = 1, ncol = 4)
  for(i in 1:length(adjs)){
    adj <- get(adjs[i], .GlobalEnv)
    length_check[i] <- ncol(adj)
  }
  
  colnames(length_check) <- adjs
  
  #Find which network has the least number of connections to the significant gene set
  idx_choice <- which(length_check == min(length_check))
  
  #Get the string of the smallest network
  small_adj <- colnames(length_check)[idx_choice]
  
  #Make a subset of the remaining networks
  adjs_sub <- adjs[which(adjs != small_adj)]
  
  adj_idx <- colnames(get(small_adj, .GlobalEnv))
  
  #Trim the other networks to match the smallest
  net_list <- list()
  for(i in 1:length(adjs_sub)){
    adj <- get(adjs_sub[i], .GlobalEnv)
    net_list[[i]] <- adj[,which(colnames(adj) %in% adj_idx)]
  }
  
  names(net_list) <- adjs_sub
  net_list$small_adj <- get(small_adj, .GlobalEnv)
  
  if(any(lengths(net_list[1:length(adjs_sub)]) != length(net_list$small_adj))){
    still_to_trim <- net_list[which(lengths(net_list[1:length(adjs_sub)]) == length(net_list$small_adj))]
    still_to_trim$small_adj <- net_list$small_adj
    
    new_small_adj <- net_list[which(lengths(net_list[1:length(adjs_sub)]) != length(net_list$small_adj))]
    new_adj_idx <- colnames(new_small_adj[[1]])
    
    trimmed_list <- list()
    for(i in 1:length(still_to_trim)){
      adj <- still_to_trim[[i]]
      trimmed_list[[i]] <- adj[,which(colnames(adj) %in% new_adj_idx)]
    }
    
    #Check to ensure all networks are the same size now
    if(any(lengths(trimmed_list) != length(new_small_adj[[1]]))){
      stop("Adjacency matrices are different sizes!")
    }else{
      print("Adjacency matrices are ready to be combined!")
    }
  }
  
  trim_order_list <- list()
  for (i in 1:length(trimmed_list)) {
    adj <- trimmed_list[[i]]
    trim_order_list[[i]] <- adj[,order(match(colnames(adj), colnames(new_small_adj[[1]])))]
  }
  
  names(trim_order_list) <- names(still_to_trim)
  trim_order_list$new_small_adj <- new_small_adj[[1]]
  
  for (i in 1:length(trim_order_list)) {
    trim_order_list[[i]] <- as.data.frame(trim_order_list[[i]])
  }
  
  #Now rbind them together!
  brain_mush <- bind_rows(trim_order_list)
  
}
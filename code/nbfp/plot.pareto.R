#This function is to process the files necessary to format the plot to get the pareto front
#Then to plot the data. The options x_lim and y_lim and rat [ratio] are SOLELY for use when is_normalized = FALSE
#is_normalized should be TRUE only if normalize = TRUE for make.table.for.pareto()
#file_name is the name to save the file as in the directory here() points to
#fig_title should change based on the data being processed
plot.pareto <- function(fin_tab = NULL, bet_than = 10, is_normalized = FALSE, x_lim = c(0,14), y_lim = c(0, 0.4), rat = 25, file_name = NULL, fig_title = NULL) {

top10 <- top_n(fin_tab, 10, fin_tab$BetterThan)
#will need to take out any intersection with pref below

# LABELS AND POINTS NEEDED FOR PLOT  
pref <- rPref::psel(fin_tab, rPref::high(logFP) * rPref::high(logPV))
subset_log <- fin_tab$GENE_SYMB %in% as.character(pref$GENE_SYMB) #need this to subset labels on pareto front

#need to take out any intersection with pref from subset top 10 = top10 names minus pref
top10naMinuspref <- top10[[10]][!(top10$GENE_SYMB %in% as.character(pref[[10]]))] #what are these 10s subsetting?

#full logical array taking out labels for the top 10 minus pref - will use this in ggplot for labels
subset_top10 <- fin_tab$GENE_SYMB %in% as.character(top10naMinuspref) 

#Maybe make main ggplot object then have conditional if() sections that will alter the ggplot settings
#dependent upon whether or not the data has been normalized.

# PLOT FOR NON-NORMALIZED THINGY
p <- ggplot(fin_tab, aes(x = logPV, y = logFP, label = GENE_SYMB)) + 
  geom_point(shape = 21, color = "gray30") + #all other gene points
  geom_step(data = pref, direction = "vh", color = "lightblue", size = 1.5) +  # pareto line
  geom_point(data = pref, shape = 20, size = 10, color = "steelblue1") + #pareto line points 
  geom_point(data = top10, shape = 20, size = 5.5, color = "orange") + #top 10 ranked (better than) points
  geom_text(aes(x = logPV, y = logFP, label= GENE_SYMB), subset(fin_tab, subset = subset_log), hjust = -0.2, vjust = -0.5, angle = 15, size = 5) + #pareto labels
  geom_text(aes(x = logPV, y = logFP, label= GENE_SYMB), subset(fin_tab, subset = subset_top10), hjust = -0.2, vjust = -0.5, angle = 15, position = "jitter", size = 5) + #top 10 labels
  ggtitle(fig_title) + 
  labs(y="-log10 SVM FP Rate", x = "-log10 Hip. Vol. p value") +
  theme(axis.title = element_text(size=30), axis.text = element_text(size = 15))
if(!(is_normalized)){
  p_coords <- p + coord_fixed(ratio = rat, xlim = x_lim, ylim = y_lim)
}else{
  p_coords <- p + coord_equal() +
    xlim(0, 1.05) + 
    ylim(0, 1.05)
  }

ggsave(file_name, plot = last_plot(), height = 10, width = 10, units = "in")


}
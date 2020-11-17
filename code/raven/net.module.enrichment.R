net.module.enrichment <- function(gene.fgn, mods, path = "."){

	require(gProfileR)
	
	u_mods <- unique(mods)
	all.enrich <- vector(mode = "list", length = length(u_mods))
	names(all.enrich) <- paste0("Module", u_mods)
	if(!file.exists("Module.Enrichment.pdf")){
		for(i in 1:length(u_mods)){
			mod.locale <- which(mods == u_mods[i])
			mod.genes <- colnames(gene.fgn)[mod.locale]
			cat("Module ", i, ": ", length(mod.genes), " Genes\n", sep = "")
			gene.locale <- match(as.numeric(mod.genes), as.numeric(gene.info[,"entrezgene"]))
			enrichment <- gprofiler(gene.info[gene.locale,1], "mmusculus", max_p_value = 0.05)
			all.enrich[[i]] <- enrichment
			}
		term.list <- lapply(all.enrich, function(x) x[,"term.name"])
		}

	pdf(paste0(path, "/Module.Enrichment.pdf"), width = 20, height = 10)
	par(mfrow = c(1,2))
	for(i in 1:length(all.enrich)){
		mod.locale <- which(mods == u_mods[i])
		mod.genes <- colnames(annotated.pos.mat)[mod.locale]
		gene.locale <- match(as.numeric(mod.genes), as.numeric(gene.info[,"entrezgene"]))
		plot.enrichment(all.enrich[[i]], 30, plot.label = paste("Module", i, ":", length(mod.genes), "Genes\nDefault Order"), order.by = "gprofiler")
		plot.enrichment(all.enrich[[i]], 30, plot.label = paste("Module", i, ":", length(mod.genes), "Genes\nOrdered by p value"), order.by = "p.value")
		}
	dev.off()

	}
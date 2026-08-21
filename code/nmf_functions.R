# These are functions to run NMF and create associated plots and tests

# Load packages

packages <- c("NMF", "nnls", "mclust", "RColorBrewer", "ggplotify",
              "pheatmap", "factoextra", "ggplot2", "ggpubr", "cowplot",
              "grid", "gridExtra")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

# Assess model performance

nmf_performance <- function(data.list, rank = 3, seed = 3, meta.data){
  # data = list of data from all labs
  # rank = rank for determining basis (3 for this study design because A, B, C groups)
  # meta.data = meta data with sample info
  
  # Function for one lab
  
  one.lab.perf <- function(data, rank, seed, meta.data){
    
    model <- nmf(data, rank, seed = seed)
    
    fitted.model <- fitted(model)
    
    # Residuals 
    res <- residuals(model)
    
    ## Normalize residuals to show fraction explained by model
    data_mat <- as.matrix(data)
    storage.mode(data_mat) <- "numeric"
    
    frac_unexplained <- res / sum((data_mat - mean(data_mat))^2)
    frac_explained <- format(1 - frac_unexplained, scientific = T)
    
    
    return(data.frame("residuals" = res, "fraction_explained_variance" = frac_explained))
  }
  
  # Apply across labs
  
  all.perf <- lapply(data.list, function(x) one.lab.perf(x, rank = rank, seed = seed, meta.data = meta.data))
  
  return(do.call(rbind, all.perf))
}

# Estimated variance (see if same/similar to residuals method, above)

nmf_evar <- function(data.list, rank = 3, seed = 3, meta.data){
  # data = list of data from all labs
  # rank = rank for determining basis (3 for this study design because A, B, C groups)
  # meta.data = meta data with sample info
  
  # Function for one lab
  
  one.lab.var <- function(data, rank, seed, meta.data){
    
    model <- nmf(data, rank, seed = seed)
    
    # Estimated variance 
    est.var <- evar(model, as.matrix(data))
    
    return(data.frame("explained_variance" = est.var))
    
  }
  
  # Apply across labs
  
  all.var <- lapply(data.list, function(x) one.lab.var(x, rank = rank, seed = seed, meta.data = meta.data))
  
  return(do.call(rbind, all.var))
  
}

# Test different rank values and output plots

nmf_rank <- function(data.list, rank.range = c(seq(2,6)), nrun = 30, randomize = FALSE,
                     plot.list = c("evar", "silhouette")){
  # data = list of data from all labs
  # rank.range = range of ranks to test
  # nrun = number of runs to do on each rank tested
  # randomize = logical indicating whether to randomize data before input
  # plot.list = list of metrics to plot. Includes "cophenetic", "dispersion", "evar", "residuals", "rss", "silhouette", "sparseness"
  
  # randomize data
  
  if(randomize == TRUE){
    
    data.random <- permute.data(data.list)
    
    # sanity check
    
    if(all(mapply(function(x, y) all.equal(colMeans(x), colMeans(y)), data.random, data.list)) == TRUE){
      
      data.list <- data.random
      
    } else {
      
      print("randomized data means do not match original data - check permutation")
      
    }
    
  }
  
  # calculate rank estimates
  
  ## function for one lab
  rank.one <- function(data, rank.range, nrun, lab){
    
    rank.est <-  nmfEstimateRank(as.matrix(data), rank.range, 
                                 method = nmf.getOption("default.algorithm"), nrun = nrun,
                                 model = NULL, verbose = FALSE, stop = FALSE)
    
    rank.plot <- plot(rank.est, y = NULL,
                      what = plot.list, 
                      na.rm = FALSE, xname = "x", yname = "y",
                      xlab = "Factorization rank", ylab = "",
                      main = paste("NMF Rank Evaluation", lab))
    
  }

  ## apply to all labs
  
  mapply(function(x, y) rank.one(x, rank.range, nrun, y), data.list, names(data.list))
  
}


# Calculate W matrix and plot heatmaps

nmf_basis <- function(data, rank = 3, meta.data, meta.col = c("group", "concentration", "both"), print.heatmap = TRUE, seed = 3){
  # data = list of data from all labs
  # rank = rank for determining basis (3 for this study design because A, B, C groups)
  # meta.data = meta data with sample info
  # meta.col = column to use for coloring heatmap
  # print.heat = logical saying whether to print out the heatmaps
  
  if(!is.null(seed)){
    
    set.seed(seed)
    
  }
  
  # subset data into labs and groups if plotting by sample concentration
  
  groups <- unique(meta.data$group)
  
  if(meta.col == "concentration"){
    
    data.new <- list()
    length(data.new) <- length(data) * length(groups)
    names(data.new) <- apply(expand.grid(names(data), groups), 1, paste, collapse="_")
    
    
    for(i in 1:length(data.new)){
      for(j in 1:length(data)){
        
        group.rows <- gsub("_.*", "", row.names(data[[j]]))
        
        for(k in 1:length(groups)){
          
          if(grepl(names(data)[[j]], names(data.new)[[i]]) == T &
             groups[k] == gsub(".*_", "", names(data.new)[[i]])){
            
            data.new[[i]] <- data[[j]][which(group.rows == groups[k]),]  
            
          }
          
        }
      }
    }
    data <- data.new
  }
  
  w <- list(length(data))
  
  if(print.heatmap == T){
    
    heats <- list(length(data))
    
    for(i in 1:length(data)){
      
      # calculate NMF
      res <- nmf(data[[i]], rank)
      
      # extract basis
      w[[i]] <- res@fit@W
      
      # plot heatmap of basis
      
      ## scale basis
      scale.w <- w[[i]] / max(w[[i]])
      colnames(scale.w) <- seq(1, ncol(scale.w), by = 1)
      
      ## set annotation row
      if(meta.col == "group" |
         meta.col == "both"){
        rows <- data.frame("group" = gsub("_.*", "", row.names(scale.w)))
      }
      
      if(meta.col == "concentration"){
        conc <- gsub("_R.*", "", row.names(scale.w))
        conc <- gsub(".*_", "", conc)
        rows <- data.frame("concentration" = conc)
      }
      
      if(meta.col == "both"){
        
        conc <- gsub("_R.*", "", row.names(scale.w))
        conc <- gsub(".*_", "", conc)
        
        rows <- data.frame(
          group = meta.data$group,
          concentration = conc,
          row.names = rownames(scale.w)
        )
      }
      
      row.names(rows) <- row.names(scale.w)
      
      # Determine rank labels
      
      if(meta.col == "group" | 
         meta.col == "both"){
        
        groups <- as.character(rows$group)
        unique.groups <- sort(unique(groups))
        n.basis <- ncol(scale.w)
        
        # mean loading per group, per basis column
        loading.matrix <- sapply(unique.groups, function(g){
          colMeans(scale.w[groups == g, , drop = FALSE])
        })
        # loading.matrix: rows = basis columns, cols = groups
        
        # for each basis column, pick the group with the highest mean loading
        basis.to.group <- unique.groups[apply(loading.matrix, 1, which.max)]
        
        # relabel heatmap columns with matched group names
        colnames(scale.w) <- basis.to.group
      }
      
      
      ## set annotation column
      cols <- data.frame("basis" = seq(1, ncol(scale.w), by = 1))
      row.names(cols) <- colnames(scale.w)
      
      ## set annotation colors
      if(meta.col == "group"){
        
        group.colors <- c("red","orange","brown","purple","gold","darkgreen", "blue")
        names(group.colors) <- levels(meta.data$group)
        
        basis.colors <- brewer.pal(n = ncol(scale.w), "Paired")
        annotation_colors <- list(group = group.colors, basis = basis.colors)
        
      }
      
      if(meta.col == "concentration"){
        
        conc.colors <- brewer.pal(length(unique(conc)), "Blues")
        names(conc.colors) <- unique(conc)
        
        basis.colors <- brewer.pal(n = ncol(scale.w), "Paired")
        annotation_colors <- list(concentration = conc.colors, basis = basis.colors)
        
      }

      if(meta.col == "both"){
        
        group.colors <- c("red","orange","brown","purple","gold","darkgreen", "blue")
        names(group.colors) <- levels(meta.data$group)
        
        conc <- gsub("_R.*", "", row.names(scale.w))
        conc <- gsub(".*_", "", conc)
        conc.colors <- brewer.pal(length(unique(conc)), "Blues")
        names(conc.colors) <- unique(conc)
        
        basis.colors <- brewer.pal(n = ncol(scale.w), "Paired")
        
        annotation_colors <- list(group = group.colors,
                                  concentration = conc.colors,
                                  basis = basis.colors)
        
      }
      
      ## plot
      
      heats[[i]] <- as.ggplot(pheatmap::pheatmap(scale.w, 
                                       color = colorRampPalette(brewer.pal(n = 9, name = "Oranges"))(100),
                                       annotation_row = rows,
                                       cluster_cols = F,
                                       annotation_colors = annotation_colors,
                                       annotation_legend = F,
                                       show_rownames = F,
                                       treeheight_col = 0,
                                       treeheight_row = 10,
                                       cellheight = 25,
                                       fontsize = 16,
                                       fontsize_row = 14,
                                       angle_col = 0,
                                       legend = F,
                                       main = names(data)[[i]])) +
                                    theme(plot.margin = margin(2, 10, 2, 10))
      
    }
      
      ## plot legends
    # make a throwaway pheatmap purely to harvest its legend grobs
    legend.src <- pheatmap::pheatmap(scale.w,
                           color = colorRampPalette(brewer.pal(n = 9, name = "Oranges"))(100),
                           annotation_row = rows,
                           cluster_cols = F,
                           cluster_rows = F,
                           fontsize = 20,
                           annotation_colors = annotation_colors,
                           annotation_legend = TRUE,
                           legend = TRUE,
                           silent = TRUE)
    
    layout.names <- legend.src$gtable$layout$name
    
    ## color-scale legend grob (usually named exactly "legend")
    color.legend.grob <- legend.src$gtable$grobs[[which(layout.names == "legend")]]
    
    ## annotation legend grob(s) (named "annotation_legend" or "annotation_legend_1", "_2", etc.
    ## if group + basis annotations each get their own)
    annot.idx <- grep("^annotation_legend", layout.names)
    annot.legend.grobs <- legend.src$gtable$grobs[annot.idx]
    
    ## stack any multiple annotation legend grobs vertically into one combined grob
    annot.legend.combined <- gtable::gtable(
      widths = unit(1, "npc"),
      heights = unit(rep(1 / length(annot.legend.grobs), length(annot.legend.grobs)), "npc")
    )
    for(j in seq_along(annot.legend.grobs)){
      annot.legend.combined <- gtable::gtable_add_grob(
        annot.legend.combined, annot.legend.grobs[[j]],
        t = j, l = 1
      )
    }
    
    # place color-scale legend TOP-RIGHT, annotation legend BOTTOM-RIGHT
    ## each occupies its own half of the legend canvas -> no overlap
    legend.plot <- ggplot() +
      theme_void() +
      annotate("text", x = 0.3, y = 1, label = "scaled basis",
               size = 7, fontface = "bold", hjust = 0, vjust = 1) +
      annotation_custom(color.legend.grob,
                        xmin = 0.3, xmax = 1.9, ymin = 0.5, ymax = 0.97) +
      annotation_custom(annot.legend.combined,
                        xmin = 0.3, xmax = 1, ymin = 0, ymax = 0.45) +
      xlim(0, 1) + ylim(0, 1)
    
    #combine 8 equal-size panels + 1 external legend column 
    panels <- plot_grid(plotlist = heats, ncol = 4, nrow = 2,
                        align = "hv", axis = "tblr",
                        rel_heights = c(1, 1))
    
    final.plot <- plot_grid(panels, legend.plot, ncol = 2, rel_widths = c(8.5, 1.5))
    
    print(final.plot)
    
  }else{
    
    for(i in 1:length(data)){
      
      # calculate NMF
      res <- nmf(data[[i]], rank)
      
      # extract basis
      w[[i]] <- res@fit@W
      
    }
    
  }
  
  # Calculate clustering performance
  
  cluster.cor <- data.frame(matrix(nrow = length(w), ncol = 1))
  row.names(cluster.cor) <- names(data)
  colnames(cluster.cor) <- "correlation"
  
  for(i in 1:length(w)){
    
    w.hcut <- hcut(w[[i]], 7)
    sample.classes <- data.frame("cluster" = w.hcut$cluster, "group" = meta.data$group)
    
    cluster.cor.test <- cor.test(as.numeric(as.factor(meta.data$group)), w.hcut$cluster, method = "pearson")
    cluster.cor[i,c("correlation")] <- cluster.cor.test$estimate
    
  }
  
  names(w) <- names(data)
  
  return(list("w" = w, "correlation" = cluster.cor))
  
}

# Validate clusters using correlations

nmf_cluster <- function(w, meta.data){
  # w = basis from nmf_basis()
  # meta.data = dataframe with sample groups
  
  # set group colors
  group.colors <- c("red","orange","brown","purple","gold","darkgreen", "blue")
  names(group.colors) <- levels(meta.data$group)
  
  group.colors.df <- data.frame(group = names(group.colors), group.color = group.colors)
  
  clusterplot <- list(length(w))
  
  for(i in 1:length(w)){
    
    w.hcut <- hcut(w[[i]], 7)
    sample.classes <- data.frame("cluster" = w.hcut$cluster, "group" = meta.data$group)
    
    clusterplot[[i]] <- ggplot(sample.classes, aes(x = group, y = cluster, color = group)) +
      geom_point(position = position_jitter(w = 0.25, h = 0), size = 2) +
      scale_color_manual(values = group.colors) +
      theme_classic() +
      scale_y_continuous(breaks = seq(min(sample.classes$cluster), max(sample.classes$cluster), by = 1)) +
      ggtitle(names(w)[[i]]) +
      ylab("Basis Cluster") +
      xlab("Sample Group")
    
  }
  
  ggarrange(plotlist = clusterplot, ncol = 4, nrow = 2, common.legend = T)
  
  
}

# Extract H matrix and important features

nmf_feats <- function(data, rank = 3, meta.data, meta.col = c("group", "concentration"), comp.data){
  # data = list of data from all labs
  # rank = rank for determining basis (3 for this study design because A, B, C groups)
  # meta.data = meta data with sample info
  # meta.col = column to use for coloring heatmap
  # comp.data = list of metabolite meta-data from all labs
  
  set.seed(3)
  
  # subset data into labs and groups if plotting by sample concentration
  
  groups <- unique(meta.data$group)
  
  if(meta.col == "concentration"){
    
    data.new <- list()
    length(data.new) <- length(data) * length(groups)
    names(data.new) <- apply(expand.grid(names(data), groups), 1, paste, collapse="_")
    
    
    for(i in 1:length(data.new)){
      for(j in 1:length(data)){
        
        group.rows <- gsub("_.*", "", row.names(data[[j]]))
        
        for(k in 1:length(groups)){
          
          if(grepl(names(data)[[j]], names(data.new)[[i]]) == T &
             groups[k] == gsub(".*_", "", names(data.new)[[i]])){
            
            data.new[[i]] <- data[[j]][which(group.rows == groups[k]),]  
            
          }
          
        }
      }
    }
    data <- data.new
  }
  
  # extract coefficients and basis
  
  h <- list(length(data))
  w <- list(length(data))
  
  for(i in 1:length(data)){
    
    # calculate NMF
    res <- nmf(data[[i]], rank)
    
    # extract coefficients
    h[[i]] <- res@fit@H
    
    # extract basis
    w[[i]] <- res@fit@W
    
  }
  
  names(h) <- names(data)
  names(w) <- names(data)
  
  # determine labels for ranks
  rank.labels <- list(length(w))
  
  for(i in 1:length(w)){
    
    ## subset data to only include a, b, and c groups and calculate sums across ranks for each group
    groups <- gsub("_.*", "", row.names(w[[i]])) 
    
    a.sum <- colSums(w[[i]][which(groups == "a"),])
    b.sum <- colSums(w[[i]][which(groups == "b"),])
    c.sum <- colSums(w[[i]][which(groups == "c"),])
    
    # assign labels to the ranks
    a.rank <- which(a.sum == max(a.sum))
    b.rank <- which(b.sum == max(b.sum))
    c.rank <- which(c.sum == max(c.sum))
    
    rank.labels[[i]] <- c("a" = a.rank, "b" = b.rank, "c" = c.rank)
    
  }    
  
  # convert coefficients to percent contribution across ranks and add bin names to results
  
  h.perc <- h
  
  for(i in 1:length(h.perc)){
    
    h.totals <- colSums(h.perc[[i]])
    
    h.perc[[i]] <- sweep(h.perc[[i]], 2, h.totals, FUN = "/")
    
    # rank names
    
    rows <- names(rank.labels[[i]][order(rank.labels[[i]])])
    
    row.names(h.perc[[i]]) <- rows
    
    # feature bins
    
    colnames(h.perc[[i]]) <- row.names(comp.data[[i]])
    
  }
  
  names(h.perc) <- names(h)
  
  return(h.perc)
  
}

# Adjusted Rand Index: calculate clustering accuracy

ari.nmf <- function(w, meta.data){
  # w = calculated basis matrix
  # meta.data = sample metadata with true groups
  
  # Prepare results
  result <- data.frame(matrix(nrow = length(w$w), ncol = 1))
  colnames(result) <- c("ARI")
  row.names(result) <- names(w$w)
  
  for(lab in 1:length(w$w)){
    
    # Determine w clusters
    
    w.hcut <- hcut(w$w[[lab]], 7)
    sample.classes <- data.frame("cluster" = w.hcut$cluster, "group" = meta.data$group)
    
    # Calculate ARI
    
    ari <- adjustedRandIndex(sample.classes$cluster, meta.data$group)
    
    # Add to results
    
    result$ARI[lab] <- ari
    
  }

  return(result)
    
}


# Permutation test: randomize metabolite abundances across samples for input into NMF

## Prepare data
permute.data <- function(data.list){
  # data.list = metabolite abundances for all labs
  
  # function for one lab
  one.lab.perm <- function(data){
    
    perm_mat <- apply(data, 2, function(col) sample(col))
    dimnames(perm_mat) <- dimnames(data)
    return(perm_mat)
    
  }
  
  # apply to all labs

  fin.list <- lapply(data.list, one.lab.perm)
  
  # sanity check
  
  sanity <- mapply(function(x, y) all.equal(colMeans(x), colMeans(y)), fin.list, data.list)
  
  if(!all(sanity == T)){
    
    print("permuted and original data don't match - check inputs")
    
  }
  
  return(fin.list)
  
}

## Calculate basis matrices and accuracy results across permutations

permute.w <- function(data.list, iterations, meta.data){
  # data.list = metabolite abundances for all labs
  # iterations = number of iterations to run
  # meta.data = sample metadata
  
  # Prepare results table
  results <- data.frame(matrix(nrow = iterations, ncol = length(data.list)*2))
  colnames(results) <- c(paste0(names(data.list), "_ARI"), paste0(names(data.list), "_Correlation"))

  for(i in 1:iterations){
    print(i)
    
    # permute data
    set.seed(3000 + i)
    perm.list <- permute.data(data.list)
    
    # calculate w
    w.i <- nmf_basis(data = perm.list, rank = 3, meta.data = meta.data, meta.col = "group", print.heatmap = F, seed = 3)
    
    # calculate ARI
    ari.i <- ari.nmf(w.i, meta.data)
    
    # add to results
    results[i, which(grepl("ARI", colnames(results)))] <- t(ari.i)
    results[i, which(grepl("Correlation", colnames(results)))] <- t(w.i$correlation)
        
  }
  
  return(results)
  
}

## Calculate summary statistics

permute.sum <- function(perm.results, metric = c("ARI", "Correlation")){
  # perm.results = results from permute.w()
  # metric = whether evaluating ARI or correlation results 
    
    sub.results <- perm.results[,which(grepl(metric, colnames(perm.results)))]
  
    means <- format(colMeans(sub.results), scientific = T)
    sds <- apply(sub.results, 2, function(x) round(sd(x), 2))
  
    labs <- gsub("_.*", "", colnames(sub.results))
    
    results <- data.frame("mean" = means, "sd" = sds, row.names = labs)    
  
    return(results)
}

## Test significance

perm.sig <- function(perm.results, real.results, metric = c("ARI", "Correlation")){
  # perm.results = results from permute.w()
  # real.results = ARI or correlation (w$correlation) results for real data
  # metric = whether evaluating ARI or correlation results
  
  # identify the appropriate columns in the permuted results
  cols <- grep(metric, colnames(perm.results), value = TRUE)
  
  # extract lab names from column names (strip "_ARI" suffix)
  lab_names <- sub(metric, "", cols)
  lab_names <- sub("_", "", lab_names)

  # calculate p-value for each lab
  p_values <- sapply(seq_along(cols), function(i) {
    
    perm_dist <- perm.results[[cols[i]]]
    real_val  <- real.results[which(row.names(real.results) == lab_names[i]),]
    
    prob <- sum(perm_dist >= real_val) / length(perm_dist)
    
    if(prob == 0){
      
      prob <- "< 1 E -3"
       
    }
    
  })
  
  results <- data.frame(Lab = lab_names, Real = real.results, P_value = p_values)
  
  return(results)
  
}

# Use training/testing split to test prediction accuracy

nmf.tt <- function(data.list, split.reps = T, split.n = 16, rank = 3, sample.seed = NULL, nmf.seed = 3, meta.data){
	# data.list = list of metabolite abundances for all labs
  # split.reps = logical indicating whether to split training/testing based on rep1/rep2 samples
	# split.n = number of samples to split into training data if not based on reps
	# rank = nmf rank
	# seed = random seed
	# meta.data = sample metadata
	
	# Split data
	
	if(split.reps == TRUE){
		
		train_idx <- which(meta.data$rep == 1)
		test_idx <- which(meta.data$rep == 2)
		
	} else {
		
	  set.seed(sample.seed)
	  
		all_idx <- seq(1:nrow(meta.data))
		train_idx <- sample(all_idx, split.n)
		test_idx <- all_idx[!(all_idx %in% train_idx)]
		
	}
	
	train.list <- lapply(data.list, function(x) x[train_idx,])
	test.list <- lapply(data.list, function(x) x[test_idx,])
	
	train.meta <- meta.data[train_idx,]
	test.meta <- meta.data[test_idx,]	
	
	# NMF on training

	nmf_train <- lapply(train.list, function(x) nmf(x, rank = rank, seed = nmf.seed))

	w_train <- lapply(nmf_train, basis)
	h_train <- lapply(nmf_train, coef)
	
	## assign training sample to dominant factor
	
	train_factors <- lapply(w_train, function(x) apply(x, 1, which.max))
	
	### apply group names to factors
	
	factor_to_type <- lapply(train_factors, function(lab) sapply(split(train.meta$group, lab), function(x) { names(sort(table(x), decreasing = T))[1] }))
	
	# Project onto testing
	
	A <- lapply(h_train, t)
	
	w_test <- Map(function(test_mat, A_mat) {
  				w <- t(apply(test_mat, 1, function(sample_row) {
    			nnls(A_mat, sample_row)$x
  				}))
  				rownames(w) <- rownames(test_mat)
  				w
				}, test.list, A)
	
	# Map to dominant factor
	
	test_factor <- lapply(w_test, function(lab) apply(lab, 1, which.max))
	predicted_type <- Map(function(f2t, tf) {
  						f2t[as.character(tf)]
						}, factor_to_type, test_factor)

	# Evaluate
								
	results.list <- Map(function(test_mat, pred_type) {
  						df <- data.frame(
    							sample = rownames(test_mat),
    							true_type = test.meta$group,
    							predicted_type = pred_type,
    							row.names = NULL
  								)
  						df$correct <- ifelse(mapply(grepl, df$predicted_type, df$true_type), 1, 0)
  						df
						}, test.list, predicted_type)
	
	accuracy <- lapply(results.list, function(lab) sum(lab$correct) / nrow(lab) * 100)
	
	mapply(function(lab, name) print(paste0("The accuracy of predictions for ", name, " is ", lab, "%" )), accuracy, names(accuracy))
	
	return(results.list)
	
}






















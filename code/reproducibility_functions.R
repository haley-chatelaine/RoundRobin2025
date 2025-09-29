# This script contains functions used in evaluating inter-lab reproducibility for the RoundRobin project

# i. concatenate list of data

concat_results <- function(list.of.results){
  # list.of.results = list of lists by lab for any results output
  
  # method if list of results is a list of vectors: 
  
  if(is.vector(list.of.results[[1]])){
    
    # label results from labs and concatenate
    
    all.result <- data.frame()
    
    ## make all results the same length (NA's will fill spots)
    
    lengths <- lapply(list.of.results, length)
    longest.index <- which.max(lengths)
    longest <- length(list.of.results[[longest.index]])
    
    for(i in 1:length(list.of.results)){
      
      length(list.of.results[[i]]) <- longest
      
      if(i == 1){
        
        all.result <- list.of.results[[i]]
        
      }else{
        
        all.result <- cbind(all.result, list.of.results[[i]])
        
      }
      
    }
    
    colnames(all.result) <- names(list.of.results)
    
  }else{
 
  # method is list is a list of data frames
    
    # label results from labs and concatenate
    
    all.result <- data.frame()
    
    for(i in 1:length(list.of.results)){
      
      list.of.results[[i]][nrow(list.of.results[[i]]) + 1,] <- names(list.of.results)[[i]]
      row.names(list.of.results[[i]])[nrow(list.of.results[[i]])] <- "lab"
      
      if(i == 1){
        
        all.result <- list.of.results[[i]]
        
      }else{
        
        all.result <- cbind(all.result, list.of.results[[i]])
        
      }
      
    }

  }
  
  return(all.result)
  
}

# ii. Only use data for features that ARE NOT duplicated in labs (for writing code for now)

no.dupes <- function(comp.found, results.list = NULL, column){
  # results.list = list of metabolite data (normalized, sc, or raw) from all labs *NOT with metabolite names mapped yet*
  # comp.found = list of compounds found results from all labs
  # column = column name used to search for duplicates (usually Name)
  
    fin.comp <- list()
    fin.results <- list()
      
    if(!is.null(results.list)){
      
      for(i in 1:length(comp.found)){
        
        # make all metabolite names upper case (especially if using "Name" as column to filter)
        
        comp.found[[i]]$Name <- toupper(comp.found[[i]]$Name)
        
        # check that metabolite order in comp.found and results are the same
        
        comp.found[[i]] <- comp.found[[i]][order(row.names(comp.found[[i]])),]
        results.list[[i]] <- results.list[[i]][,order(colnames(results.list[[i]]))]
        uniq.metabs <- comp.found[[i]]
        
        # remove any metabolites that have more than 1 value 
        
        if(all.equal(row.names(comp.found[[i]]), colnames(results.list[[i]])) == F)
          
        {print("Feature names don't match - ensure that results list uses bin names")
          
        }else{
          
          for(j in 1:length(comp.found[[i]][,c(column)])){
            if(length(which(comp.found[[i]][,c(column)] == comp.found[[i]][j,c(column)])) > 1){
              uniq.metabs.index <- which(uniq.metabs[,c(column)] == comp.found[[i]][j,c(column)])
            }else{next}
            
            if(length(uniq.metabs.index > 0)){
              uniq.metabs <- uniq.metabs[-uniq.metabs.index,]
            }
          }
          
          fin.results.list <- results.list[[i]][,which(colnames(results.list[[i]]) %in% row.names(uniq.metabs))]
          
          if(i == 1){

            fin.results[[i]] <- fin.results.list
            
          }else{

            fin.results[[length(fin.results) + 1]] <- fin.results.list          
            
          }
          
        }    
        
      }

      names(fin.results) <- names(results.list) 
      
      return(fin.results)
    }  
  
  if(is.null(results.list)){
    
    for(i in 1:length(comp.found)){
      
      # make all metabolite names upper case (especially if using "Name" as column to filter)
      
      comp.found[[i]]$Name <- toupper(comp.found[[i]]$Name)
      
      # check that metabolite order in comp.found and results are the same
      
      comp.found[[i]] <- comp.found[[i]][order(row.names(comp.found[[i]])),]
      uniq.metabs <- comp.found[[i]]
      
      # remove any metabolites that have more than 1 value 
      
      for(j in 1:length(comp.found[[i]][,c(column)])){
        if(length(which(comp.found[[i]][,c(column)] == comp.found[[i]][j,c(column)])) > 1){
          uniq.metabs.index <- which(uniq.metabs[,c(column)] == comp.found[[i]][j,c(column)])
        }else{next}
        
        if(length(uniq.metabs.index > 0)){
          uniq.metabs <- uniq.metabs[-uniq.metabs.index,]
        }
      }
        
        fin.comp.list <- uniq.metabs
        
        if(i == 1){
          
          fin.comp[[i]] <- fin.comp.list
          
        }else{
          
          fin.comp[[length(fin.comp) + 1]] <- fin.comp.list          
          
        }
        
      }    
    
    names(fin.comp) <- names(comp.found)
    
    return(fin.comp) 
  }
}

# 1. Calculate RSD of features

feat.rsd <- function(results.list){
  # results.list = list of metabolite data (normalized, sc, or raw) from all labs
    ## *metabolite data has samples in rows and metabolites in columns
  
  rsd.list <- list()
  
  for(i in 1:length(results.list)){
    
    means <- colMeans(results.list[[i]])
    sds <- apply(results.list[[i]], 2, sd)
    rsds <- sds/means * 100
    names(rsds) <- colnames(results.list[[i]])
    
    if(i == 1){
      
      rsd.list[[i]] <- rsds
    
      }else{
      
      rsd.list[[length(rsd.list) + 1]] <- rsds
    
    }
  }
  
  names(rsd.list) <- names(results.list)
  
  return(rsd.list)
  
}

# 2. Plot heatmap of RSDs of features in each lab

heat.rsd <- function(rsd.list, comp.class = NULL){
  
  # label results from labs and concatenate
  
  all.rsd <- concat_results(rsd.list)
 
  # plot heatmap with labs as columns *and super classes as annotations, if provided*
  
  if(is.null(comp.class)){
    
    pheatmap(all.rsd, 
             color = colorRampPalette(brewer.pal(n = 7, name =
                                                   "Blues"))(100),
             cluster_rows = T, 
             show_rownames = F,
             na_col = "red")
    
  }else{
    
    
    
    
  }

   
}



# 3. Create list of only metabolites that overlap (either by name or formula) across all labs to be used for parsing

overlap.metabs <- function(comp.found, column){
  # comp.found = list of compounds found results from all labs
  # results.list = list of metabolite data (normalized, sc, or raw) from all labs
  # column = column name to use for matching (molecular formula or compound name)
  
  # create list of metabolites from all labs 
  for (i in 1:length(comp.found)) {
    
    if(i == 1){
      metabs <- data.frame("metab" = comp.found[[i]][,c(column)], "lab" = names(comp.found)[[i]])
    }else{
      metabs <- rbind(data.frame("metab" = comp.found[[i]][,c(column)], "lab" = names(comp.found)[[i]]), metabs)
    }
    
  }
  
  unique.metabs <- unique(metabs$metab)
  unique.labs <- names(comp.found)
  retain.metabs <- c()
  
  for(i in 1:length(unique.metabs)){
    
    check.metab <- metabs[which(metabs$metab == unique.metabs[i]),]
    
    check.metab <- check.metab[order(check.metab$lab),]
    
    if(length(names(table(check.metab$lab))) != length(unique.labs)){next}
    
    if(all.equal(names(table(check.metab$lab)), unique.labs)){
      
      retain.metabs <- c(retain.metabs, check.metab$metab)
      
    }
    
  }
  
  retain.metabs <- unique(retain.metabs)
  
  return(retain.metabs) 
}

# 4. Calculate associations of features with lab  (linear model)

lab.assoc <- function(results.list){
  # results.list = list of metabolite data (normalized, sc, or raw) from all labs
  
  # concatenate data
  
  all.data <- concat_results(results.list)

  
  
}

# 5. Calculate associations of shared features, controlling for method variables



# 6. Calculate pairwise correlations of features (lab pairs)

lab.cor <- function(results.list, cortype){
  # results.list = list of metabolite data (normalized, sc, or raw) from all labs
  # cortype = character string of method for correlation
  
  # remove labs with no features (i.e., if features were pre-filtered out)
  cols <- sapply(results.list, ncol)
  
  if(any(cols == 0)){
    results.list <- results.list[-which(cols == 0)]
  }
  
  
  # create index of lab pairs
  
  lab.names <- names(results.list)
  lab.pairs <- combn(lab.names, 2)
  
  # create empty results list
  
  results <- list()
  
  # run correlations using indexes
  
  for(i in 1:ncol(lab.pairs)){
    
    #print(paste("comparing", lab.pairs[1,i], "and", lab.pairs[2,i]))
    
    # pull out data from lab pairs
    
    lab.a <- results.list[[lab.pairs[1,i]]]
    
    lab.b <- results.list[[lab.pairs[2,i]]]
    
    # subset data from lab pairs to only include overlapped metabolites
    
    overlap.names <- intersect(colnames(lab.a), colnames(lab.b))
    
    prep.a <- as.data.frame(lab.a[,which(colnames(lab.a) %in% overlap.names)])
    
    prep.b <- as.data.frame(lab.b[,which(colnames(lab.b) %in% overlap.names)])
    
    if(ncol(prep.a) != ncol(prep.b)){print(paste("dimensions don't match - ensure duplicated named features are filtered from datasets:", lab.pairs[1,i], lab.pairs[2,i]))}
    
    prep.a <- prep.a[,order(colnames(prep.a))]
    prep.b <- prep.b[,order(colnames(prep.b))]
    
    if(all.equal(colnames(prep.a), colnames(prep.b)) == F){print(paste("column names don't match - ensure names are consistent in:", lab.pairs[1,i], lab.pairs[2,i]))}
    
    if(all.equal(colnames(prep.a), colnames(prep.b)) == T){
      
      cor.ab <- cor(prep.a, prep.b, method = cortype)
      
    }
    
    if(i == 1){
      
      results[[i]] <- diag(cor.ab)
      names(results[[i]]) <- row.names(cor.ab)
      names(results)[[i]] <- paste0(lab.pairs[1,i], " | ", lab.pairs[2,i])
      
    }else{
      
      results[[length(results) + 1]] <- diag(cor.ab)
      names(results[[i]]) <- row.names(cor.ab)
      names(results)[[i]] <- paste0(lab.pairs[1,i], " | ", lab.pairs[2,i])
      
    }
  }
  
  # list results in a dataframe
  
    ## find unique metabs compared
    metab.names <- c()
    for(i in 1:length(results)){
      
      metab.names <- c(metab.names, unlist(names(results[[i]])))
      
    }
    metab.names <- unique(metab.names)
 
  results.df <- data.frame(matrix(ncol = length(results), nrow = length(metab.names)))
  colnames(results.df) <- names(results)
  row.names(results.df) <- metab.names
  
  for(i in 1:length(results)){
    
    if(length(results[[i]]) == 0){next}
    
    for(j in 1:nrow(results.df)){
      for(k in 1:length(results[[i]]))
        
        if(row.names(results.df)[j] == names(results[[i]])[k]){
          
          results.df[j,i] <- results[[i]][k]
          
      }
    }
  }

  
  return(results.df) 
  
}

# 7. Correlation visualization

heat.cor <- function(labcor.results, comp.class = NULL, order.by.class = FALSE, show.rownames = FALSE){
  # labcor.results = results from finding pairwise correlations of features between labs
  # comp.class = metabolite metadata with chemical classes mapped
  # order.by.class = logical stating whether to use chemical class to order the rows (as opposed to hierarchical clustering)
  # show.rownames = logical saying whether to include row names in the heatmap
  
  plot.data <- labcor.results
  
  #plot.data[is.na(plot.data)] <- 0 # make color for 0 white
  fin.data <- as.matrix(as.data.frame(apply(plot.data, 2, as.numeric)))
  row.names(fin.data) <- row.names(plot.data)
  
  # prepare colors for gradient
  
      cols <- colorRampPalette(c("darkred", "white", "darkblue"))(100)
      
  # add row annotations based on comp.class, if provided
  
  if(!is.null(comp.class)){
    
    # create list of metabolites and classes from all data
    all.class <- data.frame()
    
    i = 1
    while(i <= length(comp.class)){
      
      if(i == 1){
        comp.class[[i]]$ramp.class[is.na(comp.class[[i]]$ramp.class)] <- "unknown"
        
        all.class <- subset(comp.class[[i]], select = c(Name, retain.name, ramp.class))
      }else{
        comp.class[[i]]$ramp.class[is.na(comp.class[[i]]$ramp.class)] <- "unknown"
        
        all.class <- rbind(all.class, subset(comp.class[[i]], select = c(Name, retain.name, ramp.class)))
      }
      
      i <- i+1 
    }
    all.class <- all.class[which(duplicated(all.class$retain.name) == F),]
    
    # filter to only include metabolites in data
    all.class <- all.class[which(all.class$retain.name %in% row.names(fin.data)),]
    
    # create data frame for annotations
    row.annot <- data.frame(class = all.class$ramp.class)
    row.names(row.annot) <- all.class$retain.name
    
  }else{
    
    row.annot <- data.frame(name = row.names(fin.data))
    row.names(row.annot) <- row.names(fin.data)
    
  }
  
  # plot heatmap
    
  if(order.by.class == T){
    
    fin.data <- fin.data[order(row.names(fin.data)),]
    row.annot <- rownames_to_column(row.annot, "name")
    row.annot <- row.annot[order(row.annot$name),]
    row.names(row.annot) <- row.annot$name
    row.annot <- subset(row.annot, select = c("class"))
    
    if(!all.equal(row.names(row.annot), row.names(fin.data))){
      print("row annotations and metabolite names don't match")
    }else{
     
      hmap <- pheatmap(fin.data[order(row.annot$class),], 
                       color = cols,
                       breaks = seq(-1, 1, length.out = 101),
                       cluster_rows = F,
                       cluster_cols = F,
                       show_rownames = show.rownames,
                       annotation_row = row.annot)
       
    }
    
  }else{
    
    hmap <- pheatmap(fin.data, 
                     color = cols,
                     breaks = seq(-1, 1, length.out = 101),
                     cluster_rows = F,
                     cluster_cols = F,
                     show_rownames = show.rownames)
    
  }    
      
  return(hmap)
  
}

# Plot distributions of correlations as violin plots

violin.cor <- function(labcor.results, data.type = NULL){
  # labcor.results = results from finding pairwise correlations of features between labs
  # data.type = character string of data type for title
  
  # melt data for ggplot format
  
    labcor.results <- rownames_to_column(labcor.results, var = "feature")  
    plot.data <- reshape2::melt(labcor.results, id = "feature")
    colnames(plot.data) <- c("feature", "lab_pair", "correlation")
  
    if(any(is.na(plot.data$correlation))){
      
      plot.data <- plot.data[-which(is.na(plot.data$correlation)),]     
      
    }

    
    ggplot(plot.data, aes(x = lab_pair, y = correlation)) + 
      geom_violin() +
      geom_boxplot(width = 0.1) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90)) + 
      ggtitle(data.type)
  
}

# 8. Calculate how many formulas have isobars (duplicates) within each lab 

form.isobars <- function(comp.found){

  form.dups.list <- list()
  
  for(i in 1:length(comp.found)){
    
    formulas <- comp.found[[i]]$MF
    form.dups <- c()
    
    for(j in 1:length(formulas)){
        
        form.dups.i <- length(which(formulas == formulas[j]))
        names(form.dups.i) <- formulas[j]
        
        form.dups <- c(form.dups, form.dups.i)
      
      }
    
    if(i == 1){
      
      form.dups.list[[i]] <- form.dups
      
    }else{
      
      form.dups.list[[length(form.dups.list) + 1]] <- form.dups
      
    }
    
  }
  
  names(form.dups.list) <- names(comp.found)
  
  return(form.dups.list)
  
}

# 9. Calculate CV of each metabolite for each sample across labs

cv.bylab <- function(results.list, comp.class = NULL){
  # results.list = list of metabolite data (normalized, sc, or raw) from all labs
  # comp.class = optionally include mapping of compounds to classes for figure later
  
  # create sub-function for calculating CV
  cv.func <- function(x){
    mean.x <- mean(x)
    sd.x <- sd(x)
    return(mean.x / sd.x)
  }
  
  # create list of all metabolite names included in the data
  
  metab.names <- unique(unlist(lapply(results.list, function(x) colnames(x))))
  
  # create list of all sample names included in the data without replicate number
  
  sample.names <- unique(gsub("_rep.*", "", row.names(results.list[[1]])))
  
  # concatenate all data (for subsetting by metabolite and sample later)
  
  all.data <- concat_results(results.list)
  
  # create empty results list
  
  all.results <- data.frame(matrix(nrow = length(metab.names), ncol = 2))
  colnames(all.results) <- c("CV", "n_lab")
  row.names(all.results) <- metab.names
  
  # iterate through all metabolite names and sample names to find CV of metabolite across labs
  
  for(i in 1:length(metab.names)){
    #print(paste("metab", metab.names[i]))
    for(j in 1:length(sample.names)){
      #print(paste("sample", sample.names[j]))
      # subset data for a given metabolite from all labs that observe it
      
      subset.metab <- all.data[, which(colnames(all.data) == metab.names[i])]
      
      # calculate CV only for metabolites that were detected across multiple labs
      
      if(is.vector(subset.metab) == F){
      
      # subset data for above metabolite to only include the sample group of interest
      
      subset.sample <- subset.metab[which(grepl(sample.names[j], row.names(subset.metab))),]
      
      # convert data to a vector
      
      subset.vector <- as.numeric(unlist(subset.sample, use.names = F))
      
      # calculate the CV of the subset data
      
      subset.cv <- cv.func(subset.vector)
      
      # Add results to list
      
      all.results$CV[i] <- subset.cv
      all.results$n_lab[i] <- ncol(subset.metab)
      
      }else{
        
        # if metabolite was only detected in one lab, report NA and n_lab = 1
        
        all.results$CV[i] <- NA
        all.results$n_lab[i] <- length(which(colnames(all.data) == metab.names[i]))
      }
    }
  }
  
  # add compound classes, if provided
  
  if(is.null(comp.class) == FALSE){
    
    # extract all metabolite names and classes from comp.class list 
    
    for(k in 1:length(comp.class)){
      #print(k)
      if(k == 1){
        
        all.comp <- data.frame("metab" = comp.class[[k]]$retain.name, "class" = comp.class[[k]]$ramp.class)
      
        }else{
        
        all.comp <- rbind(all.comp, data.frame("metab" = comp.class[[k]]$retain.name, "class" = comp.class[[k]]$ramp.class))
          
      }
      
    }
    
    all.comp <- all.comp[which(duplicated(all.comp$metab) == F),]
    
    # add class to results
    
    temp.results <- rownames_to_column(all.results, var = "metab")
    
    all.results <- merge(temp.results, all.comp, by = "metab")
    
  }
  
  return(all.results)
}

# Visualize CV results

cv.vis <- function(abundance.cv, ratio.cv, all.labs = T){
  # abundance.cv = results of cv.bylab() for norm, sc, or raw data
  # ratio.cv = results of cv.bylab() for ratio data
  # all.labs = logical indicating whether to subset to only include metabolites detected in all 8 labs
  
  # subset data to only include metabolites in all 8 labs
  
  if(all.labs == T){
    
    abund.subset <- abundance.cv[which(abundance.cv$n_lab == 8),]
    abund.subset$type <- "abundance"
    abund.subset$class[is.na(abund.subset$class)] <- "unknown"
    
    ratio.subset <- ratio.cv[which(ratio.cv$n_lab == 8),]
    ratio.subset$type <- "ratio"
    ratio.subset$class[is.na(ratio.subset$class)] <- "unknown"
    
  }else{
    
    abund.subset <- abundance.cv
    abund.subset$type <- "abundance"
    abund.subset$class[is.na(abund.subset$class)] <- "unknown"
    
    ratio.subset <- ratio.cv
    ratio.subset$type <- "ratio"
    ratio.subset$class[is.na(ratio.subset$class)] <- "unknown"
    
  }
  
  # combine data to include CVs of all metabolites from each data type
  cv.data <- rbind(abund.subset, ratio.subset)
  
  # add means of CVs for each class and each data type
  
  cv.means <- aggregate(cv.data$CV, FUN = mean, by = list(class = cv.data$class, type = cv.data$type))
  
  for(i in 1:nrow(cv.data)){
    for(j in 1:nrow(cv.means)){
      
      if(cv.data$class[i] == cv.means$class[j] &
         cv.data$type[i] == cv.means$type[j]){
        
        cv.data$mean_CV[i] <- cv.means$x[j]
        
      }
    }  
  }
  

          # plot bar graph
          
          bar_plot <- ggplot(cv.data, aes(x = class, y = mean, fill = class)) + 
                          geom_bar(stat = "identity") + 
                          geom_errorbar(aes(x = class, ymin = mean - sd, ymax = mean + sd)) +
                          facet_grid(.~type) +
                          theme(axis.text.x = element_blank())
  
  # interactive bar plot with bars the size of the mean CV, include boxes for individual metabolite CVs
  
  ## create color palettes by class with gradient for CV
  ### set class colors
  class.colors <- Polychrome::createPalette(length(unique(cv.data$class)), seedcolors = "#009FF9")
  names(class.colors) <- unique(cv.data$class)
  
  ### apply gradient to lighten/darken
    cv.data$color <- NA  
  
    mean.cv <- mean(cv.data$CV)
    
    for(metab in 1:nrow(cv.data)){
      for(class in 1:length(class.colors)){
        
        if(cv.data$CV[metab] < mean.cv &
           cv.data$class[metab] == names(class.colors)[class]){
          
            cv.data$color[metab] <- lighten(palette(class.colors)[class], 
                                            amount = cv.data$CV[metab] / mean.cv,
                                            method = "relative")
        }
        
        if(cv.data$CV[metab] > mean.cv &
           cv.data$class[metab] == names(class.colors)[class]){
          
          cv.data$color[metab] <- darken(palette(class.colors)[class], 
                                          amount = mean.cv/cv.data$CV[metab],
                                         method = "relative")
        }
        
        if(cv.data$CV[metab] == mean.cv &
           cv.data$class[metab] == names(class.colors)[class]){
          
          cv.data$color[metab] <- class.colors[class]
          
          }
        }
      }
  
  color <- cv.data$color
    
  new_plot <- ggplot(cv.data, aes(x = class, y = CV, col = color)) + 
    geom_point(stat = "identity") + 
    scale_color_manual(values = c(color)) +
    facet_wrap(.~type) + 
    theme_classic() + 
    theme(legend.position = "none", axis.text.x = element_blank()) 
  
  ggplotly(new_plot)
  
  return(bar_plot)
                 
}

# function to subset lab's metabolite data to only those with standards

match.subset <- function(data, comp.list, match.type, is.data = T){
  # data = abundance results matrix for 1 lab
  # comp.list = compound meta-data for 1 lab
  # match.type = level of match type to use for subset (character)
  # is.data = logical indicating whether to return data or comp
  
  if(is.data == T){
    
    new.data <- as.data.frame(data[,which(comp.list$match.type == match.type)])
    
    return(new.data)    
    
  }else{
    
    new.comp <- comp.list[which(comp.list$match.type == match.type),]
    
    return(new.comp)
    
  }
  
}

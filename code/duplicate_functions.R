# This script contains functions for filtering duplicate annotated features into the most likely correct annotation

# Create RT indexes and relative RT indexes

rt.index <- function(comp.found){
  # comp.found = list of annotated metabolites for each lab
  
  for(i in 1:length(comp.found)){
    
    orig.order <- order(row.names(comp.found[[i]]))
    
    max.rt <- max(comp.found[[i]]$BinRT)
    
    comp.found[[i]] <- comp.found[[i]][order(comp.found[[i]]$BinRT),]
    
    comp.found[[i]]$rt.index <- order(comp.found[[i]]$BinRT)
    
    comp.found[[i]]$rt.index.perc <- comp.found[[i]]$rt.index / nrow(comp.found[[i]])
    
    comp.found[[i]]$rt.fraction <- comp.found[[i]]$BinRT/max.rt
    
    comp.found[[i]] <- comp.found[[i]][order(row.names(comp.found[[i]])),]
    
  }
  
  return(comp.found)
  
}

# Plot distributions of RT indexes for compounds with same InChiKey

rt.distrib <- function(comp.found, subset.variable, feat.retain, rt.type, lab.rts = NULL, append.name = FALSE, add.flags = FALSE){
  # comp.found = list of annotated metabolites for each lab with rt indexes calculated
  # subset.variable = column/identifier used to subset plots: Name (most specific), ID, InChiKey2D, or Formula (least specific, not recommended) 
  # feat.retain = list of features with highest pairwise correlations, calculated using cor.filter(), or highest intensities using int.filter()
  # rt.type = character string indicating "rt.index", "rt.index.perc", or "rt.fraction"
  # lab.rts = optional list of standard retention times from each lab **SAME ORDER as comp.found**
  # append.name = logical stating whether to add metabolite common name to plot titles
  # add.flags = logical stating whether to add a flag if a metabolite (a) has isobars and/or (b) has features that represent in-source fragments
  
  # add label for labs and concatenate data by rows
  
  # label results from labs and concatenate
  
  all.result <- data.frame()
  
  for(i in 1:length(comp.found)){
    
    comp.found[[i]][,ncol(comp.found[[i]]) + 1] <- names(comp.found)[[i]]
    colnames(comp.found[[i]])[ncol(comp.found[[i]])] <- "lab"
    
    if(i == 1){
      
      all.result <- comp.found[[i]]
      
    }else{
      
      all.result <- rbind(all.result, comp.found[[i]])
      
    }
    
  }
  
  # subset separate data frames for each metabolite, based on subset.variable
  
  subset.metabs <- list()
  unique.metabs <- unique(all.result[,c(subset.variable)])
  unique.metabs.names <- gsub("\\s*\\([^\\)]+\\)", "", all.result$Name[which(duplicated(all.result[,c(subset.variable)]) == F)])
  
  for(j in 1:length(unique.metabs)){
    
    subset.metabs[[j]] <- all.result[which(all.result[,c(subset.variable)] == unique.metabs[j]),]
    
  }
  
  # append metabolite names to identifiers, if desired
  
  if(append.name == F){
    
    names(subset.metabs) <- unique.metabs
    
  }else{
    
    names(subset.metabs) <- paste0(unique.metabs, "_", unique.metabs.names)
    
  }
  
  
  # remove elements without names
  
  name.metabs <- subset.metabs[!is.na(unique.metabs)]
  
  # if adding flags, create flags for metabolites with features that represent in-source fragments
  
  if(add.flags == T){
    
    # create a list of adducts that reflect NOT in-source fragments
    adducts <- c("(M-H)", "(M+H)", "(M+Na)", "(M+2H)", "(M-e)", "(M+e)", "(M+NH4)")
    
    for(m in 1:length(name.metabs)){
      
      # subset adducts from feature names
      feat.adducts <- gsub(".*? ", "", name.metabs[[m]]$Name)
      
      # create operator for finding elements that aren't in a list
      '%!in%' <- function(x,y)!('%in%'(x,y))
      
      # if any of the adducts of a feature are NOT in the adducts list (are in-source fragments), flag it
      
      if(any(feat.adducts %!in% adducts)){
        
        names(name.metabs)[[m]] <- paste0("FLAG:contains-fragments_", names(name.metabs)[[m]])
        
      }
      
    }
    
  }
  
  # only plot metabolites that have duplicate values in at least 1 lab
  
  fin.metabs <- list()
  
  for(k in 1:length(name.metabs)){
    
    if(any(duplicated(name.metabs[[k]]$lab)) == T){
      
      fin.metabs[[length(fin.metabs) + 1]] <- name.metabs[[k]]
      
      names(fin.metabs)[[length(fin.metabs)]] <- names(name.metabs)[[k]]
        
      }else{
      next
    }
    
  }
  
  # Add column for whether a feature was listed in the feat.retain list (from cor.filter() or int.filter() function)
    
    for(m in 1:length(fin.metabs)){
      #print(paste0("fin.metabs", m))
      if(is.na(names(fin.metabs)[[m]])){next}
      
      fin.metabs[[m]]$retained <- "not retained"
      
      for(n in 1:length(feat.retain)){
        #print(paste0("feat.retain", n))
        if(is.na(names(feat.retain)[[n]])){next}
        
        if(grepl(names(feat.retain)[[n]], names(fin.metabs)[[m]], ignore.case = T)){
          
          for(p in 1:nrow(fin.metabs[[m]])){
            #print(paste0("fin.metabs row", p))
            for(q in 1:nrow(feat.retain[[n]])){
              #print(paste0("feat.retain row", q))
              if(toupper(fin.metabs[[m]]$Name[p]) == toupper(feat.retain[[n]]$Feature[q]) &
                fin.metabs[[m]][p,c(rt.type)] == feat.retain[[n]][q,c(rt.type)] & 
                 fin.metabs[[m]]$lab[p] == row.names(feat.retain[[n]])[q]){
                
                fin.metabs[[m]]$retained[p] <- "retained"
                
              }
            } 
          }
        }
      }
    }
  
  # if list of standards was supplied, add row to fin.metabs for each metabolite that had a standard for each lab
  
    if(!is.null(lab.rts)){
      
      for(r in 1:length(fin.metabs)){
      
        if(is.na(names(fin.metabs)[[r]])){next}
        
        # match InChiKey2D and lab name between fin.metabs and lab.rts for features with standards
        # add data from lab.rts to fin.metabs and mark as a "standard"
        
        for(s in 1:length(lab.rts)){
          
          lab.name <- names(lab.rts)[[s]]
          
          for(q in 1:nrow(lab.rts[[s]])){
            
            if(grepl(lab.rts[[s]]$InChiKey2D[q], names(fin.metabs)[[r]]) &
               !is.na(lab.rts[[s]]$Retention[q])){
                
                  #print(r)
                
                  fin.metabs[[r]] <- add_row(fin.metabs[[r]],
                                             # replace the adduct name from lab.rts with the adduct found by the lab
                                             "Name" = gsub("\\([^\\)]+\\)", paste0("(", lab.rts[[s]]$`Adduct type`[q], ")"), lab.rts[[s]]$Name[q]),
                                             "ID" = fin.metabs[[r]]$ID[1],
                                             "MF" = lab.rts[[s]]$Formula[q],
                                             "BinC12mz" = lab.rts[[s]]$`M/Z`[q],
                                             "BinRT" = lab.rts[[s]]$Retention[q],
                                             "InChiKey2D" = lab.rts[[s]]$InChiKey2D[q],
                                             "IDMechanism" = "lab_standard",
                                             
                                             # use the rt.index in comp.found that's closest to the RT of lab.rts
                                             "rt.index" = min(comp.found[[s]]$rt.index[which(abs(comp.found[[s]]$BinRT - lab.rts[[s]]$Retention[q]) == min(abs(comp.found[[s]]$BinRT - lab.rts[[s]]$Retention[q])))]),
                                             
                                             # use the rt.index.perc in comp.found that's closest to the RT of lab.rts
                                             "rt.index.perc" = min(comp.found[[s]]$rt.index.perc[which(abs(comp.found[[s]]$BinRT - lab.rts[[s]]$Retention[q]) == min(abs(comp.found[[s]]$BinRT - lab.rts[[s]]$Retention[q])))]),
                                             
                                             # use the rt.fraction in comp.found that's closest to the RT of lab.rts *IS THIS RIGHT?
                                             #"rt.fraction" = min(comp.found[[s]]$rt.fraction[which(abs(comp.found[[s]]$BinRT - lab.rts[[s]]$Retention[q]) == min(abs(comp.found[[s]]$BinRT - lab.rts[[s]]$Retention[q])))]),
                                             
                                             "rt.fraction" = lab.rts[[s]]$Retention[q] / max(comp.found[[s]]$BinRT),
                                             
                                             "lab" = lab.name,
                                             
                                             "retained" = "standard"
                  )
          }
        }
      }
    }
  }
  
  # plot scatter plots for each named metabolite

  plot.fun <- function(x, y) { x <- arrange(x, retained)
  
                              colors <- c("grey", "violet", "blue")
                              names(colors) <- c("not retained", "retained", "standard")
                              
                              mean.line <- mean(x[,c(rt.type)])
                              min.line <- mean.line - 1 * sd(x[,c(rt.type)])
                              max.line <- mean.line + 1 * sd(x[,c(rt.type)])
    
                              ggplot(x, aes(x = lab, y = x[,c(rt.type)], col = retained, label = Name)) + 
                                              #geom_point() +
                                              geom_point(position = position_jitter(w = 0.25, h = 0)) +
                                              scale_color_manual(values = c(colors)) +
                                              theme_classic() +
                                              theme(axis.text.x = element_text(angle = 90)) +
                                              geom_hline(yintercept = mean.line, alpha = 0.5) +
                                              geom_hline(yintercept = min.line, linetype = "dotted", alpha = 0.3) +
                                              geom_hline(yintercept = max.line, linetype = "dotted", alpha = 0.3) +
                                              ggtitle(y) }
  
  plotly.fun <- function(x) { ggplotly(x) }
  
  plots <- mapply(plot.fun, fin.metabs, names(fin.metabs), SIMPLIFY = FALSE)
  
  plotlys <- mapply(plotly.fun, plots, SIMPLIFY = FALSE)
  
  return(plotlys)
}


# Calculate pairwise correlations of each feature, grouped by an umbrella term*, across labs
  ## *umbrella term: some taxon higher than Name (ID, InChi, or Formula)

lab.assoc <- function(results.list, comp.found, umbrella){
  # results.list = list of CF results from all labs
  # comp.found = list of annotations fo CF results from all labs
  # umbrella = some taxon higher than Name (ID, InChi, or Formula)
  
  # label results from labs and concatenate ### NEED TO FIX so that all.comp and all.results can be merged
  
  ## comp.found
  
    all.comp <- data.frame()
    
    for(i in 1:length(comp.found)){
      
      comp.found[[i]][,ncol(comp.found[[i]]) + 1] <- names(comp.found)[[i]]
      colnames(comp.found[[i]])[ncol(comp.found[[i]])] <- "lab"
      
      if(i == 1){
        
        all.comp <- comp.found[[i]]
        
      }else{
        
        all.comp <- rbind(all.comp, comp.found[[i]])
        
      }
      
    }

  ## results.list
    
    all.results <- data.frame()
    
    for(i in 1:length(results.list)){
      
      results.list[[i]] <- t(results.list[[i]])
      
      if(i == 1){
      
        all.results <- results.list[[i]]
        
      }else{
        
        all.results <- rbind(all.results, results.list[[i]])
        
      }
      
    }
 
  #  all.comp and all.results have metabolites in the same order by definition
    # because names for results.list files were assigned based on order of comp.found in map.metabs() summary function
    
  # For each ID group of features, find association of feature with lab
    
    id.groups <- unique(all.comp[,c(umbrella)])
    
    assoc.list <- list()
    
    for(i in 1:length(id.groups)){
      #print(id.groups[i])
      
      if(is.na(id.groups[i])){next}
      
      # extract data for features of a given umbrella term
      
      feats.results <- t(all.results[which(all.comp[,c(umbrella)] == id.groups[i]),])
      comp.results <- all.comp[which(all.comp[,c(umbrella)] == id.groups[i]),]
      
      # if multiple species for a given umbrella are found, find correlations
    
      if(nrow(feats.results) > 1){
        
        # set up assoc.list for results
        
        assoc.list[[i]] <- data.frame(matrix(nrow = ncol(feats.results), ncol = ncol(feats.results)))
        .rowNamesDF(assoc.list[[i]], make.names = T) <- paste0(comp.results$Name, "_", comp.results$lab)
        colnames(assoc.list[[i]]) <- paste0(comp.results$Name, "_", comp.results$lab) 
        
        # calculate pairwise correlations for each feature
        
        pairs <- expand.grid(rep(1:ncol(feats.results)), rep(1:ncol(feats.results)))
        
        for(j in 1:nrow(pairs)){
          
          assoc.list[[i]][pairs[j,1], pairs[j,2]] <- cor(feats.results[,pairs[j,1]], feats.results[,pairs[j,2]])
          
        }
        
      }else{
        
      # if only one species is available for an umbrella term, give it a correlation value of 1
        
        assoc.list[[i]] <- data.frame(matrix(nrow = 1, ncol = 1))
        .rowNamesDF(assoc.list[[i]], make.names = T) <- paste0(comp.results$Name, "_", comp.results$lab)
        colnames(assoc.list[[i]]) <- paste0(comp.results$Name, "_", comp.results$lab)
        
      assoc.list[[i]][1, 1] <- 1
      
      }
      
    }
    
    names(assoc.list) <- id.groups
    
    return(assoc.list)
     
}

# Select feature with highest correlations for duplicates found in a given lab

cor.filter <- function(feat.cor){
  # feat.cor = list of pairwise correlations generated using lab.assoc()
  
  retain.metabs <- list()
  
  for(i in 1:length(feat.cor)){
    #print(i)
    
    # skip features that don't have an ID
    
    if(is.na(names(feat.cor)[[i]])){next}
    
    # extract lab names from feature correlations
    
    lab.names <- unique(gsub(".*_", "", colnames(feat.cor[[i]])))
    
    # create results matrix
    
    retain.metabs[[i]] <- data.frame(matrix(nrow = length(lab.names), ncol = 2))
    row.names(retain.metabs[[i]]) <- lab.names
    colnames(retain.metabs[[i]]) <- c("Feature", "Mean Correlation")
    
    for(j in 1:length(lab.names)){
      
      # subset list of features for each lab under a given ID
      
      compare.metabs <- feat.cor[[i]][,which(grepl(lab.names[j], colnames(feat.cor[[i]])))]
      
      # if only one feature is listed for a lab, extract its name and mean correlation
      
      if(is.vector(compare.metabs)){
 
        retain.metabs[[i]][j,1] <- gsub("_.*", "", colnames(feat.cor[[i]])[which(grepl(lab.names[j], colnames(feat.cor[[i]])))])
        retain.metabs[[i]][j,2] <- mean(compare.metabs)      
        
      }else{
      
      # if multiple features were listed for a lab, select the feature with the highest average correlation
        
        compare.means <- colMeans(compare.metabs)
        
        # if there are ties (multiple features with the exact same correlation value), select the first instance
        
        if(length(which(compare.means == max(compare.means))) > 1){

          retain.metabs[[i]][j,1] <- gsub("_.*", "", colnames(compare.metabs)[which(compare.means == max(compare.means))[1]])
          retain.metabs[[i]][j,2] <- compare.means[which(compare.means == max(compare.means))[1]]          
          
        }else{
          
          retain.metabs[[i]][j,1] <- gsub("_.*", "", colnames(compare.metabs)[which(compare.means == max(compare.means))])
          retain.metabs[[i]][j,2] <- compare.means[which(compare.means == max(compare.means))]
          
        }
      }
    }
  }
  
  names(retain.metabs) <- names(feat.cor)
  
  return(retain.metabs)
  
}

# Calculate mean intensities of all features under an umbrella term and mark whether they should be retained

int.calculate <- function(results.list, comp.found.rt, umbrella, rt.type){
  # results.list = list of CF results from all labs
  # comp.found.rt = list of annotations fo CF results from all labs with rt.index mapped
  # umbrella = some taxon higher than Name (ID, InChi, or Formula)
  # rt.type = style of rt to use (rt.index, rt.index.perc, rt.fraction)
  
  # label results from labs and concatenate
  
    ## comp.found
    
    all.comp <- data.frame()
    
    for(i in 1:length(comp.found.rt)){
      
      comp.found.rt[[i]][,ncol(comp.found.rt[[i]]) + 1] <- names(comp.found.rt)[[i]]
      colnames(comp.found.rt[[i]])[ncol(comp.found.rt[[i]])] <- "lab"
      
      if(i == 1){
        
        all.comp <- comp.found.rt[[i]]
        
      }else{
        
        all.comp <- rbind(all.comp, comp.found.rt[[i]])
        
      }
      
    }
  
    ## results.list
    
    all.results <- data.frame()
    
    for(i in 1:length(results.list)){
      
      results.list[[i]] <- t(results.list[[i]])
      
      if(i == 1){
        
        all.results <- results.list[[i]]
        
      }else{
        
        all.results <- rbind(all.results, results.list[[i]])
        
      }
      
    }
  
  # all.comp and all.results have metabolites in the same order by definition
  # because names for results.list files were assigned based on order of comp.found.rt in map.metabs() summary function

  # For each ID group of features, find feature with highest mean intensity
  
  id.groups <- unique(all.comp[,c(umbrella)])
  
  int.list <- list()
  
  for(i in 1:length(id.groups)){
    #print(i)
    
    if(is.na(id.groups[i])){next}
    
    # extract data for features of a given umbrella term
    
    feats.results <- t(all.results[which(all.comp[,c(umbrella)] == id.groups[i]),])
    comp.results <- all.comp[which(all.comp[,c(umbrella)] == id.groups[i]),]
    
    # for multiple species of a given umbrella term
    
    if(nrow(feats.results) > 1){
      # set up int.list for results
      
      int.list[[i]] <- data.frame(matrix(ncol = ncol(feats.results), nrow = 5))
      colnames(int.list[[i]]) <- paste0(comp.results$Name, "_", comp.results$lab)
      row.names(int.list[[i]]) <- c("mean_intensity", "retain", rt.type, "feature", "lab")
      
      # calculate mean intensities for each feature
      
      for(j in 1:ncol(int.list[[i]])){
        
        int.list[[i]][c("mean_intensity"),j] <- mean(feats.results[,j])
        
      }
      
      # mark the highest mean as retained for each lab
      
      labs <- unique(gsub(".*_", "", colnames(int.list[[i]])))
      
      for(k in 1:length(labs)){
        
        compare_index <- which(grepl(labs[k], colnames(int.list[[i]])))
        
        int.list[[i]][c("retain"), which(int.list[[i]][c("mean_intensity"), ] == max(int.list[[i]][c("mean_intensity"), compare_index]))] <- "X"
        # try creating a separate data frame for each lab, mark the highest mean, then match based on column name in the int.list
      }
      
      # add rt.indexes for all features
      
      int.list[[i]][c(rt.type),] <- comp.results[,c(rt.type)]
      
      # subset feature name
      
      int.list[[i]][c("feature"),] <- gsub("_.*", "", colnames(int.list[[i]]))
      
      # subset lab name
      
      int.list[[i]][c("lab"),] <- gsub(".*_", "", colnames(int.list[[i]]))
      
    }else{
      # set up int.list for results
      
      int.list[[i]] <- data.frame(matrix(ncol = 1, nrow = 3))
      colnames(int.list[[i]]) <- paste0(comp.results$Name, "_", comp.results$lab)
      row.names(int.list[[i]]) <- c("mean_intensity", "retain", rt.type)
      
      # calculate mean intensity
      
      int.list[[i]][c("mean_intensity"),] <- mean(feats.results)
      
      # mark retained
      
      int.list[[i]][c("retain"),] <- "X"
      
      # add rt.indexes for all features
      
      int.list[[i]][c(rt.type),] <- comp.results[,c(rt.type)]
      
      # subset feature name
      
      int.list[[i]][c("feature"),] <- gsub("_.*", "", colnames(int.list[[i]]))
      
      # subset lab name
      
      int.list[[i]][c("lab"),] <- gsub(".*_", "", colnames(int.list[[i]]))
      
    }
  }
  
  names(int.list) <- id.groups
  
  return(int.list)

}

# Filter list of intensity features to only retain highest mean

int.filter <- function(int.list, rt.type){
  # int.list = list of feature intensity means and "X" retain from int.calculate()
  # rt.type = style of rt to use (rt.index, rt.index.perc, rt.fraction)
  
  retain.metabs <- list()
  
  # loop through each InChiKey
  for(i in 1:length(int.list)){
    #print(i)
    
    if(is.na(names(int.list)[[i]])){next}
    
    # determine which labs are represented
    labs <- unique(gsub(".*_", "", colnames(int.list[[i]])))
    
    # create results matrix
    retain.metabs[[i]] <- data.frame(matrix(nrow = length(labs), ncol = 3))
    row.names(retain.metabs[[i]]) <- labs
    colnames(retain.metabs[[i]]) <- c("Feature", rt.type, "mean")
    
    # add retained metabs to results 
    for(j in 1:nrow(retain.metabs[[i]])){
      #print(j)
      
      retain.metabs[[i]][j, "Feature"] <- gsub("\\_.*", "", colnames(int.list[[i]])[which(grepl(row.names(retain.metabs[[i]])[j], colnames(int.list[[i]])) &
                                                                          int.list[[i]][c("retain"),] == "X")])
      
      retain.metabs[[i]][j, c(rt.type)] <- int.list[[i]][c(rt.type), which(grepl(row.names(retain.metabs[[i]])[j], colnames(int.list[[i]])) &
                                                                                      int.list[[i]][c("retain"),] == "X")]
      
      retain.metabs[[i]][j, "mean"] <- int.list[[i]][c("mean_intensity"), which(grepl(row.names(retain.metabs[[i]])[j], colnames(int.list[[i]])) &
                                                                          int.list[[i]][c("retain"),] == "X")]
      
    }
    
  }
  
  names(retain.metabs) <- names(int.list)
  
  return(retain.metabs)
  
}

# Do the filtering: highest mean intensity + RT index in majority

rt.int.filter.1 <- function(comp.found.rt, int.list, results.list, return.comp.found.fin = FALSE){
  # comp.found.rt = list of found metabolites for each lab with rt.index mapped using rt.index()
  # results.list = list of abundances from each lab (norm, sc, or raw)
  # int.list = list of features, organized by an umbrella term, with mean intensities from int.calculate()
  # return.comp.found.fin = logical indicating whether the object returned should be the filtered results list (if FALSE) or the filtered comp.found list (if TRUE)
  
  # change all metabolite names to all caps (for matching later)
      ## results and comp.found
      for(i in 1:length(comp.found.rt)){
        
        comp.found.rt[[i]][,c("Name")] <- toupper(comp.found.rt[[i]][,c("Name")])
        colnames(results.list[[i]]) <- toupper(colnames(results.list[[i]]))
        
      }
  
      ## int.list and int.retain
      for(i in 1:length(int.list)){
        
        if(length(int.list[[i]]) == 0){next}
        
        colnames(int.list[[i]]) <- toupper(colnames(int.list[[i]]))
        int.list[[i]][c("feature"),] <- toupper(int.list[[i]][c("feature"),])
        
      }

  # filter: see annotations for the steps
  
  for(i in 1:length(int.list)){
    
    # Determine majority RT index for each umbrella term
    ## "majority" = median of rt indexes
    
    rt.majority <- median(as.numeric(int.list[[i]][c("rt.index"),]))
    
    # mark the features to retain in comp.found, based off of intensity
    
    for(j in 1:length(comp.found.rt)){
      
      # if feature isn't found in a lab, move to the next lab
      
      if(length(which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])) == 0){next}
  
      # subset the feature of interest for the lab
        
        if(length(which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])) == 1){
          
          # if only one feature is found in a lab, mark to retain it by adding a column of the name without adduct to comp.found
          
          int.feats <- int.list[[i]][,which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])]
          names(int.feats) <- row.names(int.list[[i]])
          
          retain.feat <- int.feats[c("feature")]
          
          comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat), c("retain.name")] <- gsub(" \\(.*", "", retain.feat)
          
        }else{
          
        # if multiple features are found in a lab, subset the name of the retained (highest intensity) feature
          
          int.feats <- int.list[[i]][,which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])]
          
          retain.feat <- int.feats[c("feature"), which(int.feats[c("retain"),] == "X")]
          
        }
      
      # if only one species is listed as the feature name with the highest intensity, 
        # mark to retain it by adding a column of the name without adduct to comp.found
      
      if(length(which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])) > 1 &
        length(which(comp.found.rt[[j]][,c("Name")] == retain.feat)) == 1){
        
        comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat), c("retain.name")] <- gsub(" \\(.*", "", retain.feat) 
        
      }
      
      # if multiple species are listed as the feature name with the highest intensity, 
        # find the one closest to the rt.majority, then
        # mark to retain it by adding a column of the name without adduct 
      
      if(length(which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])) > 1 &
         length(which(comp.found.rt[[j]][,c("Name")] == retain.feat)) > 1){
        
        # find differences between rt.indexes and rt.majority 
        
        rt.difs <- apply(int.feats[c("rt.index"),], 2, function(x) abs(as.numeric(x) - rt.majority))
        
        # only look at features with the "retained" name (highest intensity)
        
        rt.retain.difs <- rt.difs[which(int.feats[c("feature"),] == retain.feat)]
        
        # find the feature with the minimum difference from rt.majority, and extract index from rt.difs
        
        rt.dif.retain.min <- which(int.feats[c("feature"),] == retain.feat & 
                                     rt.difs == min(rt.retain.difs))
        
        # if both features with the same retained name are equidistant from the rt.majority,
          # i.e., if one lab only observes two features, but they have the same name and different RTs
        # then select the feature with the lowest rt.index
        
        if(length(rt.dif.retain.min > 1)){
          
          rt.retain <- which(int.feats[c("rt.index"),] == min(as.numeric(int.feats[c("rt.index"),])))
          
          comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat &
                                     comp.found.rt[[j]][,c("rt.index")] == int.feats[c("rt.index"), rt.retain]), c("retain.name")] <- gsub(" \\(.*", "", retain.feat)
          
        }else{
          
          # otherwise, mark to retain the feature with the lowest rt.dif by adding a column of the name without adduct to comp.found
          
          comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat &
                                     comp.found.rt[[j]][,c("rt.index")] == int.feats[c("rt.index"), rt.dif.retain.min]), c("retain.name")] <- gsub(" \\(.*", "", retain.feat)
          
        }
        
      }
      
    }    

  }
  
  # select and change names of features in results/comp.found so that all labs only have features with names that would match
  ## **this assumes comp.found.rt and results.list have the same order--I can't check it in the function because of duplicate column names in results
  ### ***if you add metabolite names as columns using map.metabs() in the summary and add rt.index using rt.index() and do nothing else to adjust the order, this assumption holds
  #### ****CHECK THE OBJECTS if you do anything else
  
  if(return.comp.found.fin == FALSE){

    fin.results <- list()
    
    for (i in 1:length(results.list)) {
      
      fin.results[[i]] <- results.list[[i]]
      
      colnames(fin.results[[i]]) <- comp.found.rt[[i]]$retain.name
      
      fin.results[[i]] <- fin.results[[i]][,-which(is.na(comp.found.rt[[i]]$retain.name))] 
      
    }
    
    names(fin.results) <- names(results.list)
    
    return(fin.results)    
    
  }else{
    
    fin.comp <- list()
    
    for (i in 1:length(comp.found.rt)) {
      
      fin.comp[[i]] <- comp.found.rt[[i]]
      
      fin.comp[[i]] <- fin.comp[[i]][-which(is.na(comp.found.rt[[i]]$retain.name)),] 
      
    }
    
    names(fin.comp) <- names(comp.found.rt)
    
    return(fin.comp)
  }
  
}
  
# Do the filtering: highest mean intensity + RT index in majority (Ewy's logic)

rt.int.filter.2 <- function(comp.found.rt, int.list, results.list, rt.type, lab.rts = NULL, return.comp.found.fin = FALSE, ratio.list = NULL){
  # comp.found.rt = list of found metabolites for each lab with rt.index mapped using rt.index()
  # results.list = list of abundances from each lab (norm, sc, or raw)
  # int.list = list of features, organized by an umbrella term, with mean intensities from int.calculate()
  # rt.type = style of rt to use (rt.index, rt.index.perc, rt.fraction)
  # lab.rts = optional list of standard compounds and RTs for matching
  # return.comp.found.fin = logical indicating whether the object returned should be the filtered results list (if FALSE) or the filtered comp.found list (if TRUE)
  # ratio.list = list of ratio data to filter so that it includes the same metabolites as the results.list data
  
  # change all metabolite names to all caps (for matching later)
      ## results and comp.found
      for(i in 1:length(comp.found.rt)){
        
        comp.found.rt[[i]]$Name <- toupper(comp.found.rt[[i]]$Name)
        colnames(results.list[[i]]) <- toupper(colnames(results.list[[i]]))
        
      }
      
      ## int.list and int.retain
      for(i in 1:length(int.list)){
        
        if(length(int.list[[i]]) == 0){next}
        
        colnames(int.list[[i]]) <- toupper(colnames(int.list[[i]]))
        int.list[[i]][c("feature"),] <- toupper(int.list[[i]][c("feature"),])
        
      }
  
      ## lab.rts
      if(!is.null(lab.rts)){
        for(i in 1:length(lab.rts)){
          
          lab.rts[[i]]$Name <- toupper(lab.rts[[i]]$Name)
          
        }
      }
  
  # filter: see annotations for the steps
  
  for(i in 1:length(int.list)){
    
    # Determine majority and sd of RT index for each umbrella term
        ## "majority" = mean of rt indexes
        
        rt.majority <- mean(as.numeric(int.list[[i]][c(rt.type),]))
        
        ## find window +/- 1 standard deviation of rt indexes
        
        rt.window.min <- rt.majority - 1 * sd(as.numeric(int.list[[i]][c(rt.type),]))
        rt.window.max <- rt.majority + 1 * sd(as.numeric(int.list[[i]][c(rt.type),]))
    
    # mark the features to retain in comp.found, based off of intensity and RT
        # also mark criterion by which selected:
          # 1 = only feature in InChiKey
          # 2 = standard
          # 3a = highest intensity within 1 sd of rt
          # 3b = highest intensity not within 1 sd of rt
    
    for(j in 1:length(comp.found.rt)){
      
      # if feature isn't found in a lab, move to the next lab
      
      if(length(which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])) == 0){next}
      
      # subset the feature of interest for the lab
      
      if(length(which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])) == 1){
        
        # if only one feature is found in a lab, mark to retain it by adding a column of the name without adduct to comp.found
        
        int.feats <- int.list[[i]][,which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])]
        names(int.feats) <- row.names(int.list[[i]])
        
        retain.feat <- int.feats[c("feature")]
        
        comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat), c("retain.name")] <- gsub(" \\(.*", "", retain.feat) 
        
        comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat), c("match.type")] <- "1"
        
      }else{
        
        # if multiple features are found in a lab, subset them and mark name and rt.index of the retained feature for later
        
        int.feats <- int.list[[i]][,which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])]
        
        retain.feat <- int.feats[c("feature"), which(int.feats[c("retain"),] == "X")]
        retain.rt <- int.feats[c(rt.type), which(int.feats[c("retain"),] == "X")]
        
      }
      
      # if the species with the highest intensity is also within 1 sd of RT Majority for a lab, 
      # mark to retain it by adding a column of the name without adduct to comp.found
      
      if(length(which(int.list[[i]][c("lab"),] == names(comp.found.rt)[[j]])) > 1){
        
        # if no species fall within 1 sd of RT Majority, mark the highest intensity feature for retention
        
        if(all(rt.window.min > int.feats[c(rt.type),] | int.feats[c(rt.type),] > rt.window.max)){
          
          comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat &
                                     comp.found.rt[[j]][,c(rt.type)] == retain.rt), c("retain.name")] <- gsub(" \\(.*", "", retain.feat)
          
          comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat &
                                     comp.found.rt[[j]][,c(rt.type)] == retain.rt), c("match.type")] <- "3b"
          
        }
        
        # determine if feature with max intensity falls within 1 sd of RT Majority
        
        if(all(rt.window.min > int.feats[c(rt.type),] | int.feats[c(rt.type),] > rt.window.max) == FALSE &
          (rt.window.min < retain.rt & retain.rt < rt.window.max) == TRUE){
          
            comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat &
                                       comp.found.rt[[j]][,c(rt.type)] == retain.rt), c("retain.name")] <- gsub(" \\(.*", "", retain.feat)
            
            comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat &
                                       comp.found.rt[[j]][,c(rt.type)] == retain.rt), c("match.type")] <- "3a"
          
        }
        
        if(all(rt.window.min > int.feats[c(rt.type),] | int.feats[c(rt.type),] > rt.window.max) == FALSE &
          (rt.window.min < retain.rt & retain.rt < rt.window.max) == FALSE){
      # if the species with the highest intensity is NOT within 1 sd of the RT Majority for a lab,
      # mark the most intense feature that is within 1 sd of the RT Majority for retention
          
          rt.feats <- data.frame(int.feats[,which(rt.window.min < int.feats[c(rt.type),] & int.feats[c(rt.type),] < rt.window.max)])
          row.names(rt.feats) <- row.names(int.feats)
          
          max.int.rt.feats <- max(as.numeric(rt.feats[c("mean_intensity"),]))
         
          if(max.int.rt.feats == -Inf){print(paste("missing values in mean intensity", i, j))}
          
          retain.feat <- rt.feats[c("feature"), which(rt.feats[c("mean_intensity"),] == max.int.rt.feats)]
          retain.rt <- rt.feats[c(rt.type), which(rt.feats[c("mean_intensity"),] == max.int.rt.feats)]
          
          comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat &
                                     comp.found.rt[[j]][,c(rt.type)] == retain.rt), c("retain.name")] <- gsub(" \\(.*", "", retain.feat)
          
          comp.found.rt[[j]][which(comp.found.rt[[j]][,c("Name")] == retain.feat &
                                     comp.found.rt[[j]][,c(rt.type)] == retain.rt), c("match.type")] <- "3a_rt"
          
          }
        }
      } 
    }    
    

  # select and change names of features in results/comp.found so that all labs only have features with names that would match
  ## **this assumes comp.found.rt and results.list have the same order--I can't check it in the function because of duplicate column names in results
  ### ***if you add metabolite names as columns using map.metabs() in the summary and add rt.index using rt.index() and do nothing else to adjust the order, this assumption holds
  #### ****CHECK THE OBJECTS if you do anything else
  
  if(return.comp.found.fin == FALSE & !is.null(ratio.list)){
    
    ratio.results <- list()
    
    for (i in 1:length(ratio.list)) {
      
      ratio.results[[i]] <- ratio.list[[i]]
      
      colnames(ratio.results[[i]]) <- comp.found.rt[[i]]$retain.name
      
      ratio.results[[i]] <- ratio.results[[i]][,-which(is.na(comp.found.rt[[i]]$retain.name))] 
      
    }
    
    names(ratio.results) <- names(ratio.list)
    
    return(ratio.results)    
    
  }
  
  
    if(return.comp.found.fin == FALSE & is.null(ratio.list)){
    
    fin.results <- list()
    
    for (i in 1:length(results.list)) {
      
      fin.results[[i]] <- results.list[[i]]
      
      colnames(fin.results[[i]]) <- comp.found.rt[[i]]$retain.name
      
      fin.results[[i]] <- fin.results[[i]][,-which(is.na(comp.found.rt[[i]]$retain.name))] 
      
    }
    
    names(fin.results) <- names(results.list)
    
    return(fin.results)    
    
  }
  
  if(return.comp.found.fin == TRUE & is.null(ratio.list)){
    
    fin.comp <- list()
    
    for (i in 1:length(comp.found.rt)) {
      
      fin.comp[[i]] <- comp.found.rt[[i]]
      
      fin.comp[[i]] <- fin.comp[[i]][-which(is.na(comp.found.rt[[i]]$retain.name)),] 
      
    }
    
    names(fin.comp) <- names(comp.found.rt)
    
    return(fin.comp)
  }
  
}

# replace retained features with the appropriate feature, based on standard data, if needed

std.replace <- function(data.fin, data.orig, comp.fin, comp.orig, std.data, type.return){
  # data.fin = filtered data using rt.int.filter.2()
  # data.orig = original data used to filter using rt.int.filter.2()
  # comp.fin = compound meta data filtered using rt.int.filter.2()
  # comp.orig = compound meta data used to filter using rt.int.filter.2()
  # std.data = data frame of appropriate features to retain based on manual std inspection
  # type.return = data or comp.found - type of data to return
  
  
  for(i in 1:length(data.fin)){
    
    # replace with proper feature, if needed
    
    for(j in 1:nrow(std.data)){
      for(k in 1:ncol(data.fin[[i]])){
        
        if(names(data.fin)[[i]] == std.data$Lab[j] &
           colnames(data.fin[[i]])[k] == std.data$Metabolite[j]){
          
          real.retain.index <- which(comp.orig[[i]]$InChiKey2D == std.data$InChiKey2D[j] &
                                      round(comp.orig[[i]]$rt.fraction, 3) == std.data$rt.fraction[j] &
                                       grepl(std.data$Feature[j], comp.orig[[i]]$Name, fixed = T)) 
          
          print(paste("lab", i, "std metab", j, "orig metab", k, "new data:", comp.orig[[i]]$Name[real.retain.index]))
          
          data.fin[[i]][,k] <- data.orig[[i]][,real.retain.index]
          
          comp.fin[[i]][k,] <- comp.orig[[i]][real.retain.index,]
          comp.fin[[i]]$Name[k] <- toupper(comp.fin[[i]]$Name[k])
          comp.fin[[i]]$match.type[k] <- "2"                                 
          
        } 
      }
    }
  }
  
  if(type.return == "data"){
    return(data.fin)
  }
  
  if(type.return == "comp.found"){
    return(comp.fin)
  }
}

# Match features to standards, if available (replace an annotation with that closest to the standard)

std.match <- function(data.fin, data.orig, comp.fin, comp.orig, lab.rts, ppm = 20, rt = 0.1, type.return){
  # data.fin = filtered data using rt.int.filter.2()
  # data.orig = original data used to filter using rt.int.filter.2()
  # comp.fin = compound meta data filtered using rt.int.filter.2()
  # comp.orig = compound meta data used to filter using rt.int.filter.2()
  # lab.rts = list of rts of standards available for each lab
  # ppm = ppm bound for matching features with standard m/z
  # rt = window around rt.fraction, given as decimal
  # type.return = data or comp.found - type of data to return
  
  # make a list of the features closest to standards, if available for each lab
  
  for(i in 1:length(data.fin)){
    
    std.list <- data.frame(matrix(ncol = ncol(comp.orig[[i]]) + 1))
    colnames(std.list) <- c(colnames(comp.orig[[i]]), "meanInt")
    
    for(j in 1:nrow(lab.rts[[i]])){

      if(is.na(lab.rts[[i]]$Retention[j])){next}
      
      if(!is.na(lab.rts[[i]]$Retention[j])){
        
        for(k in 1:ncol(data.fin[[i]])){
          
          #print(paste0("Lab ", i, " std ", j, " metab ", k))
          
          # match InChiKey
          
          if(comp.fin[[i]]$InChiKey2D[k] == lab.rts[[i]]$InChiKey2D[j]){
            
            # find all of the features with m/z within ppm of the standard 
              # can't use adduct because would need to escape +H in grepl
            
            feats <- comp.orig[[i]][which(comp.orig[[i]]$InChiKey2D == lab.rts[[i]]$InChiKey2D[j]),]
            
            if(nrow(feats) == 1){
              
              feats$meanInt <- mean(data.orig[[i]][,which(comp.orig[[i]]$InChiKey2D == lab.rts[[i]]$InChiKey2D[j])])
              
            }else{
              
              feats$meanInt <- colMeans(data.orig[[i]][,which(comp.orig[[i]]$InChiKey2D == lab.rts[[i]]$InChiKey2D[j])])
              
            }
            
            ## calculate lower m/z bound
            min.mz <- (-(lab.rts[[i]]$`M/Z`[j]) * (ppm/1000000)) + lab.rts[[i]]$`M/Z`[j]
            
            ## calculate upper m/z bound
            max.mz <- (lab.rts[[i]]$`M/Z`[j] * (ppm/1000000)) + lab.rts[[i]]$`M/Z`[j]
            
            fin.feats <- feats[which(min.mz < feats$BinC12mz &
                                       feats$BinC12mz < max.mz),]
            
            if(nrow(fin.feats) == 0){
              
              print(paste0(comp.fin[[i]]$InChiKey2D[k], " has no m/z values within ", ppm, " ppm of standard for ", names(comp.fin)[[i]]))
              
            }
            
          # match the RT to the standards (closest and within rt% window)
            
            ## make rt.fraction for standards
            lab.rts[[i]]$rt.fraction <- lab.rts[[i]]$Retention / max(comp.orig[[i]]$BinRT)
            
            ## calculate RT lower bound 
            min.rt <- lab.rts[[i]]$rt.fraction[j] - rt
            
            ## calculate RT upper bound
            max.rt <- lab.rts[[i]]$rt.fraction[j] + rt
            
            std.difs <- abs(lab.rts[[i]]$Retention[j] - fin.feats$BinRT)
            
            ## find the feature corresponding to the RT match (closest to std and within percentage dif from standard)      
            index <- which(std.difs == min(std.difs) &
                             min.rt < fin.feats$rt.fraction &
                             fin.feats$rt.fraction < max.rt)
            
            if(length(index) > 0){
              
              retain.feat <- fin.feats[index,]
              
              std.list <- add_row(std.list, retain.feat)
            
            }else{
              
              print(paste0(comp.fin[[i]]$InChiKey2D[k], " has no rt.fraction values within ", rt, " of standard rt.fraction for ", names(comp.fin)[[i]]))
              
            }
        
            ## double-check logic with other labs (check that if std is missing it's bc adduct isn't matched in comp.fin)
            
          }
        } 
      }
    }
    
    ## if two standards for the same InChiKey, retain the feature with the highest intensity in data
    
    std.list$remove <- NA
    
    for(m in 1:nrow(std.list)){
      
      if(duplicated(std.list$InChiKey2D)[m] == T){
        
        std.list$remove[which(std.list$meanInt == min(std.list$meanInt[which(std.list$InChiKey2D == std.list$InChiKey2D[m])]))] <- "X"
        
      }
      
    }
    
    if(sum(is.na(std.list$remove)) < nrow(std.list)){
      
      std.list <- std.list[-(which(std.list$remove == "X")),]
      
    }

    std.list <- std.list[-which(is.na(std.list$Name)),]
    
    # replace data.fin and comp.fin with the retained features (replace 3a/3b matchy type with 2)
   
    if(nrow(std.list) > 0){
      
      for(n in 1:nrow(std.list)){
        for(p in 1:nrow(comp.fin[[i]])){
          
          if(std.list$InChiKey2D[n] == comp.fin[[i]]$InChiKey2D[p]){
            
            # find index of correct metab in original data
            
            retain.feat.index <- which(comp.orig[[i]]$BinC12mz == std.list$BinC12mz[n] &
                                         comp.orig[[i]]$BinRT == std.list$BinRT[n])
            
            comp.fin[[i]][p,] <- comp.orig[[i]][retain.feat.index,]
            
            row.names(comp.fin[[i]])[p] <- row.names(comp.orig[[i]])[retain.feat.index]
            
            comp.fin[[i]]$retain.name[p] <- toupper(gsub(" \\(.*", "", comp.fin[[i]]$retain.name[p]))
            
            comp.fin[[i]][p,c("match.type")] <- "2"
            
            data.fin[[i]][,p] <- data.orig[[i]][,retain.feat.index]
            
            colnames(data.fin[[i]])[p] <- comp.fin[[i]]$retain.name[p]
            
          }
        }
      }  
    }
  }
  
  if(type.return == "data"){
    return(data.fin)
  }
  
  if(type.return == "comp.found"){
    return(comp.fin)
  }
  
  
}

# plot the distribution of match types across labs

match.distrib <- function(comp.matched){
  # comp.matched = metabolite meta-data with match types mapped
  
  # add label for labs and concatenate data by rows
  
  # label results from labs and concatenate
  
    all.result <- data.frame()
    
    for(i in 1:length(comp.matched)){
      
      comp.matched[[i]][,ncol(comp.matched[[i]]) + 1] <- names(comp.matched)[[i]]
      colnames(comp.matched[[i]])[ncol(comp.matched[[i]])] <- "lab"
      
      if(i == 1){
        
        all.result <- comp.matched[[i]]
        
      }else{
        
        all.result <- rbind(all.result, comp.matched[[i]])
        
      }
      
    }
 
    all.result <- subset(all.result, select = c(retain.name, match.type, lab))  
    
  # plot heatmap  
  
  plot <- ggplot(data = all.result, aes(x = lab, y = retain.name, fill = match.type)) +
    geom_tile() +
    theme(axis.text.y = element_blank())

    return(ggplotly(plot))      
}


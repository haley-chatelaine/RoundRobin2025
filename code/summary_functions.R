# These are functions used to create summaries of CF results across labs

# Load packages

packages <- c("tidyverse", "dplyr", "ggplot2", "ggpattern", "ggvenn", "tidyr", "UpSetR", "plotly", "devtools",
              "RColorBrewer", "ggsci", "ggpubr", "pheatmap", "grafify", "pheatmap",
              "openxlsx", "data.table", "factoextra", "htmltools", "colorspace")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

if (!requireNamespace("ncats/RAMP-DB", quietly = TRUE)) {
  library(devtools)
  install_github("ncats/RAMP-DB")
}

if("pmp" %in% installed.packages() == F){
  BiocManager::install("pmp")  
}

library(pmp)
library(RaMP)
invisible(lapply(packages, library, character.only = TRUE))

# Read files

read.results <- function(path, file.text, sheet.text){
  # path = path containing folders with files from labs
  # file.text = character string of pattern for file names
  # sheet.text = character string of sheet name
  
  file.list <- list.files(path = path, "Lab", full.names = T, recursive = F)
  
  results <- list()
  
  for(i in 1:length(file.list)){
    
    iroa.file <- list.files(file.list[i], file.text, full.names = T, recursive = T)
    
    if( i == 1 ){
      
     results[[i]] <- read.xlsx(iroa.file, sheet = sheet.text, rowNames = T, na.strings = "0")
      
    }else{
      
      results[[length(results) + 1]] <- read.xlsx(iroa.file, sheet = sheet.text, rowNames = T, na.strings = "0")
      
    }
    
  }
  
  names(results) <- list.dirs(path = path, full.names = F, recursive = F)
 
  return(results) 
}

# Map sample names to data

map.samples <- function(path, file.text, sheet.text, master.map, results){
  # path = path where lab results folders are
  # file.text = character string pattern to grepl file
  # sheet.text = character string of experimental design sheet name
  # master.map = object with master sample mapping
  # results = list of results to map file names to
  
  lab.names <- list.files(path = path, "Lab", full.names = T, recursive = F)
  
  for(i in 1:length(lab.names)){
    print(paste("Lab", i))
    # file path
    iroa.file <- list.files(lab.names[i], file.text, full.names = T, recursive = T)
    
    # sample map data for dataset/files
    samp.map <- read.xlsx(iroa.file, sheet = sheet.text, rowNames = F)
    samp.map.data <- subset(samp.map, select = c(dataId, SampleId, ExperimentalGroup))
    samp.map.data <- samp.map.data[!(is.na(samp.map.data$dataId)),]
    samp.map.file <- subset(samp.map, select = -c(dataId, SampleId, ExperimentalGroup)) 
    
    # merge data and files
    samp.master <- merge(samp.map.data, samp.map.file, by = "SampleId", no.dups = T)
    
    # add description to map 
    samp.master$description <- NA
    for(j in 1:nrow(samp.master)){
      #print(j)
      for(k in 1:nrow(master.map)){
        #print(k)
        if(samp.master$sampleNum[j] == master.map$sampleNum[k]){
          
          samp.master$description[j] <- master.map$description[k]
          
        }
      }
    }
    
    # change sample names in data to fit descriptions
    for(m in 1:nrow(results[[i]])){
      #print(m)
      for(n in 1:nrow(samp.master)){
        #print(n)
        if(row.names(results[[i]])[m] == samp.master$dataId[n]){
          row.names(results[[i]])[m] <- samp.master$description[n]
        }
      }
      
      results[[i]] <- results[[i]][order(row.names(results[[i]])),]
      
    }
  }
  
  return(results)
  
}

# Map metabolite names to results

map.metabs <- function(comp.found, results.list, column){
  # comp.found = list of compound annotations for each lab
  # results.list = list of metabolite results to map metabolite names to
  # column = column used to map feature annotations (usually "Name")
  
  for(i in 1:length(results.list)){
    
    if(all.equal(colnames(results.list[[i]]), rownames(comp.found[[i]]))){
      
      colnames(results.list[[i]]) <- comp.found[[i]][,c(column)]
      
      # if any column names are repeated, append a counter so that the column name is fixed 
      ## this way, if a metabolite gets filtered out, its name stays the same, rather than getting a new automatic counter from R
      
      for(j in ncol(results.list[[i]]):1){
        
        if(duplicated(colnames(results.list[[i]]))[j]){
          
          counter <- length(which(colnames(results.list[[i]]) == colnames(results.list[[i]])[j]))
          #if(counter > 2){print(j)}
          
          colnames(results.list[[i]])[j] <- paste0(colnames(results.list[[i]])[j], "_", counter)
          
        }
        
      }
      
      
    }else{
      print("compounds found and results list don't match - check metabolite orders")
    }
    
  }
  
  return(results.list)
  
}

# Map sample classes to metabolites

## *update to include db = rampDB in all ramp calls

map.classes <- function(comp.found, class.level = "super_class"){
  # comp.found = list of results files with metabolites found by ClusterFinder
  # class.level = ClassyFire hierarchy level to map classes to 
  
  # Format metabolite names for RaMP queries
  
  for(i in 1:length(comp.found)){
    
    comp.found[[i]]$ramp.name <- NA
    
    for(j in 1:nrow(comp.found[[i]])){
      
      if(grepl("HMDB", comp.found[[i]]$ID[j]) == T){
        comp.found[[i]]$ramp.name[j] <- paste0("hmdb:", comp.found[[i]]$ID[j])
      }
      
      if(grepl("CHEBI", comp.found[[i]]$ID[j]) == T){
        comp.found[[i]]$ramp.name[j] <- gsub("CHEBI", "chebi", comp.found[[i]]$ID[j])
      }

      if(grepl("LM", comp.found[[i]]$ID[j]) == T){
        comp.found[[i]]$ramp.name[j] <- paste0("LIPIDMAPS:", comp.found[[i]]$ID[j])
      }
      
      if(grepl("HMDB", comp.found[[i]]$ID[j]) == F & 
         grepl("CHEBI", comp.found[[i]]$ID[j]) == F &
         grepl("LM", comp.found[[i]]$ID[j]) == F){
        comp.found[[i]]$ramp.name[j] <- "none"
      }
      
    }
    
  }
  
  
  
  # Run chemical class lookup and retrieve only the ClassyFire_super_class values
  
  ramp.col = paste0("ClassyFire_", class.level)
  
  for(i in 1:length(comp.found)){
    
    ramp.classes <- chemicalClassSurvey(mets = comp.found[[i]]$ramp.name)
    
    class.table <- as.data.frame(ramp.classes$met_classes)
    
    comp.found[[i]]$ramp.class <- NA
    
    for(j in 1:nrow(class.table)){
      for(k in 1:nrow(comp.found[[i]])){
        
        if(comp.found[[i]]$ramp.name[k] == class.table$sourceId[j] &
           class.table$class_level_name[j] == ramp.col){
          
          comp.found[[i]]$ramp.class[k] <- class.table$class_name[j]
          
        }      
        
      }
      
    }
    
  }
  
  comp.found[[i]]$ramp.class[is.na(comp.found[[i]]$ramp.class)] <- "unknown"
  
  return(comp.found)
  
}

# Change 1 to 0 so all counted as missing values

one.to.zero <- function(results.list){
  
  for(i in 1:length(results.list)){
    for(j in 1:ncol(results.list[[i]])){
      
      results.list[[i]][,j] <- ifelse(results.list[[i]][,j] == 1, 0, results.list[[i]][,j])
      
    }   
  }
  
  return(results.list)
}

# Find numbers of missing values per metabolite

n.miss <- function(results.list, report = F){
  # results.list = list of metabolite results, potentially with classes mapped
  
  miss.results <- list()
  
  for(i in 1:length(results.list)){
    
    miss.results[[i]] <- data.frame(matrix(nrow = ncol(results.list[[i]]), ncol = 4))
    colnames(miss.results[[i]]) <- c("miss_0", "miss_1", "total", "total_perc")
    row.names(miss.results[[i]]) <- paste0(colnames(results.list[[i]]), "_", order(colnames(results.list[[i]])))
    
    miss.results[[i]]$miss_0 <- apply(results.list[[i]], 2, function(x) length(which(x == 0)))
    miss.results[[i]]$miss_1 <- apply(results.list[[i]], 2, function(x) length(which(x == 1)))
    miss.results[[i]]$total <-  apply(results.list[[i]], 2, function(x) (length(which(x == 0)) + length(which(x == 1))))
    miss.results[[i]]$total_perc <- miss.results[[i]]$total / (nrow(results.list[[i]]) * ncol(results.list[[i]])) * 100
    
  }
  
  names(miss.results) <- names(results.list)
  
  if(report == T){
    
    results.table <- data.frame(matrix(nrow = length(miss.results), ncol = 5))
    colnames(results.table) <- c("miss_0", "miss_1", "total", "total_data", "total_perc_miss")
    row.names(results.table) <- names(miss.results)
    
    for(j in 1:length(miss.results)){
      
      results.table$miss_0[j] <- sum(miss.results[[j]]$miss_0)
      results.table$miss_1[j] <- sum(miss.results[[j]]$miss_1)
      results.table$total[j] <- sum(miss.results[[j]]$total)
      results.table$total_data[j] <- nrow(results.list[[j]]) * ncol(results.list[[j]])
      results.table$total_perc_miss[j] <- round(results.table$total[j] / results.table$total_data[j], 2)
      
    }
    
    print(results.table)
    
  }
  
  return(miss.results)
  
}

# Find number of missing values per sample and create lab-specific meta-data with # NAs

sample.miss <- function(results.list, meta.data){
  # results.list = list of metabolite abundances in samples WITHOUT missing values imputed
  # meta.data = general sample meta data dataframe created using create.meta() (below)
  
  lab.meta <- list()
  
  for(i in 1:length(results.list)){
    
    # prep lab meta data list
    lab.meta[[i]] <- meta.data
    
    # count missing values per sample
    miss <- apply(results.list[[i]], 1, function(x) length(which(x == 0)))
    
    # add numbers of misisng values per sample to lab meta
    lab.meta[[i]][,c("missing")] <- miss
    
  }
  
  names(lab.meta) <- names(results.list)
  
  return(lab.meta)
  
}

# Remove metabolites that aren't present in both replicate samples

rep.remove <- function(results.list, comp.found = NULL, ratio.list = NULL){
  # results.list = list of metabolite results with sample names mapped
  # comp.found = list of compound annotations (only include if filtering comp.found)
  # ratio.data = list of C12/C13 ratio data that matches the results list 
    # used to filter according to missing values in abundance data (as of June 12, 2024)
  
  filter.list <- list()
  filter.comp <- list()
  filter.ratio <- list()
  
  for(i in 1:length(results.list)){
    
    # set sample names
    samples <- gsub("_rep.*", "", row.names(results.list[[i]]))
    
    # create empty list for metabolites to remove
    
    metabs.remove <- c()
    
    for(j in 1:length(samples)){
      
      # select data for the replicates
      
      metabs.check <- results.list[[i]][which(grepl(samples[j], row.names(results.list[[i]]))), ]

      for(k in 1:ncol(results.list[[i]])){
        
        # get the column index for any metabolite that has any kind of missing value in at least one rep
        
        if(length(which(metabs.check[,k] == 0)) > 0 |
           length(which(metabs.check[,k] == 1)) > 0){
          
          metabs.remove <- c(metabs.remove, k)
          #metabs.remove <- c(metabs.remove, colnames(metabs.check)[k])
          
              }
            }
          }
      
      # remove duplicates and NAs from list
      metabs.remove <- unique(metabs.remove)
      
      # remove metabolites named in the list from data and comp.found
      
      filter.list[[i]] <- subset(results.list[[i]], select = -metabs.remove)
      #print(paste(length(metabs.remove), "metabolites removed for not being present in both", samples[j], "samples", "in", names(results.list)[[i]]))
      
      if(!is.null(comp.found)){
        
        filter.comp[[i]] <- comp.found[[i]][-metabs.remove,]
        
      }
      
      if(!is.null(ratio.list)){
        
        # check that columns of ratio.list and results.list match
          
          if(all.equal(colnames(ratio.list[[i]]), colnames(results.list[[i]])) == F){
            print(paste("ratio and results columns don't match for", i))
          }else{
          
            filter.ratio[[i]] <- subset(ratio.list[[i]], select = -metabs.remove) 
          
          }
        
      }
      
  }
  
  if(is.null(comp.found) & is.null(ratio.list)){
    
    names(filter.list) <- names(results.list)
    
    return(filter.list)
    
  }
  
  if(!is.null(comp.found)){
    
    names(filter.comp) <- names(comp.found)
    
    return(filter.comp)    
    
  }
  
  if(!is.null(ratio.list)){
    
    names(filter.ratio) <- names(ratio.list)
    
    return(filter.ratio)    
    
  }
  
  
}

# Filter by high percentage of NAs

na.perc.remove <- function(results.list, cutoff, comp.found = NULL){
  # results.list = list of metabolite results with sample names mapped and all missing values as 0
  # cutoff = threshold percentage of NAs as a decimal
  # comp.found = list of compound annotations (only include if filtering comp.found)
  
  filter.list <- list()
  comp.list <- list()
  
  for(i in 1:length(results.list)){
    
    # count NAs
    na.count <- apply(results.list[[i]], 2, function (x) length(which(x == 0)))
    na.perc <- na.count / nrow(results.list[[i]])
    
    # select which to remove
    metabs.remove <- which(na.perc > cutoff)
    print(paste(length(metabs.remove), "removed from", names(results.list)[[i]]))
    
    # filter
    if(length(metabs.remove > 0)){
      
      filter.list[[i]] <- results.list[[i]][,-metabs.remove]
      
    }else{
      
      filter.list[[i]] <- results.list[[i]]
      
    }
    
    
    if(!is.null(comp.found)){
      
      if(length(metabs.remove > 0)){
        
        comp.list[[i]] <- comp.found[[i]][-metabs.remove,]
        
      }else{
        
        comp.list[[i]] <- comp.found[[i]]
      }  
    }
  }
  
  if(is.null(comp.found)){
    
    names(filter.list) <- names(results.list)
    
    return(filter.list)
    
  }else{
    
    names(comp.list) <- names(comp.found)
    
    return(comp.list)
  }
  
}

# Impute "0" missing values using half-minimum

impute.halfmin <- function(results.list){
  # results.list = list of CF results data frames
  
  mv.results <- list()
  
  for(i in 1:length(results.list)){
    
    mv.results[[i]] <- results.list[[i]]
    
    mv.results[[i]][mv.results[[i]] == 0] <- NA
    
    imp.values <- apply(mv.results[[i]], 2, FUN = function(x) min(x, na.rm = T) / 2)
    
    for(j in 1:ncol(mv.results[[i]])){
      
      if(any(is.na(mv.results[[i]][,j]))){
        
        mv.results[[i]][which(is.na(mv.results[[i]][,j])), j] <- imp.values[j]
        
      }else{next}
      
    }  
    
  }
  
  names(mv.results) <- names(results.list)
  
  return(mv.results)
  
}

# Impute "0" missing values using KNN

impute.knn <- function(results.list){
  # results.list = list of CF results data frames
  
  mv.results <- list()
  
  for(i in 1:length(results.list)){
    
    mv.results[[i]] <- results.list[[i]]
    
    mv.results[[i]][mv.results[[i]] == 0] <- NA
    
    mv.results[[i]] <- pmp::mv_imputation(mv.results[[i]], method="KNN", k=5, rowmax=0.5, colmax=0.5, maxp=NULL, check_df=F)
      
  }
  
  names(mv.results) <- names(results.list)
  
  return(mv.results)
  
}

# Impute "0" values with KNN and "1" values with half-minimum

impute.both <- function(results.filter, results.original){
  # results.filter = results list that has been filtered for high NAs/presence in replicates
  # results.original = original results list containing 1s and 0s for missing values

  
  for(i in 1:length(results.filter)){
    print(i)
    # if value was originally 0, set value to NA and impute using KNN
    
    ## order results.original by column name so matching works (no idea why this works)
    
    results.original[[i]] <- results.original[[i]][,order(colnames(results.original[[i]]))]
    
    for(j in 1:ncol(results.filter[[i]])){
      for(k in 1:ncol(results.original[[i]])){
       
        if(colnames(results.filter[[i]])[j] == colnames(results.original[[i]])[k] & 
           any(results.original[[i]][,k] == 0)){
          print(j)
          results.filter[[i]][which(results.original[[i]][,k] == 0),j] <- NA
        
          } 
      }
    }
  
    results.filter[[i]] <- pmp::mv_imputation(results.filter[[i]], method="KNN", k=5, rowmax=0.5, colmax=0.5, maxp=NULL, check_df=F)
  
    # if value was originally 1, impute using half-minimum
    for(j in 1:ncol(results.filter[[i]])){
      for(k in 1:ncol(results.original[[i]])){
        
        if(colnames(results.filter[[i]])[j] == colnames(results.original[[i]])[k] &
           any(results.original[[i]][,k] == 1)){
          
          missing.index <- which(results.original[[i]][,k] == 1)
          
          results.filter[[i]][missing.index, j] <-  NA
          
          results.filter[[i]][missing.index,j] <- min(results.filter[[i]][,j], na.rm = T) / 2
          
          } 
        }
      }
    
  }
  
  names(results.filter) <- names(results.original)
  
  return(results.filter)
  
}

# Count numbers of annotated bins

count.bins <- function(comp.found){
  # comp.found = list of results files with metabolites found by ClusterFinder
  
  results <- data.frame(matrix(nrow = length(comp.found), ncol = 8))
  colnames(results) <- c("Total Annotated Bins", "POS Annotated Bins", "NEG Annotated Bins",
                         "Database", "Library (MS1/RT)", "LTRS Reference Library", "Molecular Formula", "Unknown")
  row.names(results) <- names(comp.found)
  
  for(i in 1:length(comp.found)){
    
    results[i, c("Total Annotated Bins")] <- nrow(comp.found[[i]])
    results[i, c("POS Annotated Bins")] <- length(which(grepl("pos",row.names(comp.found[[i]]), ignore.case = T)))
    results[i, c("NEG Annotated Bins")] <- length(which(grepl("neg",row.names(comp.found[[i]]), ignore.case = T)))
    results[i, c("Database")] <- length(which(grepl("database",comp.found[[i]]$IDMechanism, ignore.case = T)))
    results[i, c("Library (MS1/RT)")] <- length(which(grepl("library",comp.found[[i]]$IDMechanism, ignore.case = T)))
    results[i, c("LTRS Reference Library")] <- length(which(grepl("LTRS",comp.found[[i]]$IDMechanism, ignore.case = T)))
    results[i, c("Molecular Formula")] <- length(which(grepl("formula",comp.found[[i]]$IDMechanism, ignore.case = T)))
    results[i, c("Unknown")] <- length(which(grepl("unknown",comp.found[[i]]$IDMechanism, ignore.case = T)))
    
  }
  
  return(results)
  
}

# Create "core dataset" of only detected metabolites
## *remove metabolites that have any missing values (represented as values < 1) in the raw data

core.data <- function(raw.results, results.list = NULL, comp.found = NULL){
  # raw.results = list of IROA Clean results from which to find metabolites with original missing values
  # results.list = list of metabolite matrices from each lab with sample names mapped and features (as Bin- name) in columns 
  # comp.found = list of results files with metabolites found by ClusterFinder
  
  raw.filter <- list()
  results.filter <- list()
  
  for(i in 1:length(raw.results)){
    
    feat.nas <- apply(raw.results[[i]], 2, FUN = function(X) length(which(X < 1)))
    feat.remove.index <- which(feat.nas > 0)
    
    print(paste0(names(raw.results)[[i]], 
                 " removed ",
                 length(feat.remove.index),
                 " features with at least 1 missing value out of ",
                 ncol(raw.results[[i]]),
                 " total features"))
    
    # if you supplied a results list, remove features with missing values from the raw data
    
    if(!is.null(results.list)){
      
      if(i == 1){
        
        results.filter[[i]] <- results.list[[i]][,-feat.remove.index]
        
      }else{
        
        results.filter[[length(results.filter) + 1]] <- results.list[[i]][,-feat.remove.index]
        
      }
      
    }else{
      
      # if you didn't supply a results list, filter the raw results only
      
      if(i == 1){
        
        raw.filter[[i]] <- raw.results[[i]][,-feat.remove.index]
        
      }else{
        
        raw.filter[[length(raw.filter) + 1]] <- raw.results[[i]][,-feat.remove.index]
        
      }
      
    }
    
  }
  
  # if you supply a comp.found list, filter those compounds
  
  if(!is.null(comp.found)){
    
    comp.filter <- list()
    
    for(i in 1:length(comp.found)){
      
      if(all.equal(row.names(comp.found[[i]]), colnames(raw.results[[i]]))){
        
        feat.nas <- apply(raw.results[[i]], 2, FUN = function(X) length(which(X < 1)))
        feat.remove.index <- which(feat.nas > 0) 
        
        if(i == 1){
          
          comp.filter[[i]] <- comp.found[[i]][-feat.remove.index,]
          
        }else{
          
          comp.filter[[length(comp.filter) + 1]] <- comp.found[[i]][-feat.remove.index,]
          
        }
        
      }else{
        
        print("feature names don't match between compounds found and results - ensure bins are named BIN- in results and are ordered in both")
        
      }
      
    }
    
  }
  
  # return the necessary list
  
  if(is.null(results.list) & is.null(comp.found)){
    
    names(raw.filter) <- names(raw.results)
    return(raw.filter)
    
  }
  
  if(!is.null(results.list) & is.null(comp.found)){
    
    names(results.filter) <- names(results.list)
    return(results.filter)
    
  }
  
  if(is.null(results.list) & !is.null(comp.found)){
    
    names(comp.filter) <- names(comp.found)
    return(comp.filter)
    
  }
}

# Determine metabolites with high intensity (these might be driving correlations)

high.metabs <- function(results.list, percentile){
  # results.list = list of metabolite results, potentially with classes mapped
  # percentile = threshold of highest metabolite intensities to subset
  
  subset.list <- list()
  
  for(i in 1:length(results.list)){
    
    mean.ints <- colMeans(results.list[[i]])
    
    subset.i <- colnames(results.list[[i]])[which(mean.ints > quantile(mean.ints, percentile))]
    
    if(i == 1){
      subset.list[[i]] <- subset.i
    }else{
      subset.list[[length(subset.list) + 1]] <- subset.i
    }
    
  }
  
  names(subset.list) <- names(results.list)
  return(subset.list)
  
}

# Subset results to exclude metabolites with high intensity (above)

outlier.remove <- function(results.list, high.metabs.result){
  # results.list = list of metabolite results, potentially with classes mapped
  # high.metabs.result = list of high intensity metabolites from high.metabs() function
  
  subset.list <- list()
  
  for(i in 1:length(results.list)){
    
    subset.i <- results.list[[i]][,-which(colnames(results.list[[i]]) %in% high.metabs.result[[i]])]
    
    if(i == 1){
      subset.list[[i]] <- subset.i
    }else{
      subset.list[[length(subset.list) + 1]] <- subset.i
    }
    
  }
  
  names(subset.list) <- names(results.list)
  return(subset.list)
  
}

# Create boxplots of total intensities across samples

intensity.box <- function(results.list, data.type, log2 = F){
  # results.list = list of CF results data frames
  # data.type = character string of type of data (normalized, sc, raw)
  # log2 = logical of whether to log2 transform data or not
  
  # Create total intensity dataframe
  
  total.int <- data.frame()
  
  for(i in 1:length(results.list)){
    
    sample.tics <- as.numeric(rowSums(results.list[[i]]))
    samples <- row.names(results.list[[i]])
    lab <- rep(names(results.list)[[i]], length(sample.tics))
    
    lab.data <- cbind(lab, samples, sample.tics)
    
    total.int <- rbind(total.int, lab.data)
  }
  
  total.int$lab <- as.factor(total.int$lab)
  total.int$sample.tics <- as.numeric(total.int$sample.tics)
  
  # plot
  
  ggplot(total.int, aes(x = lab, y = log2(sample.tics), fill = lab)) + 
    geom_boxplot() +
    theme(text = element_text(size = 20)) +
    #ylim(0, 35) +
    theme_classic() +
    labs(title = paste(data.type, "data"), y = "log2(sum of feature intensities)")
}

# Separate sample meta-data

create.meta <- function(map.results){
  # map.results = list of results with descriptive sample names mapped
  
  # ensure all lists of results have the same row names
  all.labs <- combn(length(map.results), 2)
  for(i in 1:ncol(all.labs)){
    if(all.equal(row.names(map.results[[c(all.labs[1,i])]]), 
                 row.names(map.results[[c(all.labs[2,i])]])) == F){
      print("sample orders don't match - check results that all samples are included and in the same order in all results sheets")
    }
  }
  
  meta <- data.frame(matrix(nrow = nrow(map.results[[1]]), ncol = 3))
  colnames(meta) <- c("group", "concentration", "rep")
  row.names(meta) <- row.names(map.results[[1]])
  
  meta$group <- as.factor(gsub("_.*", "", row.names(meta)))
  meta$concentration <- as.factor(read.table(text = as.character(row.names(meta)), sep = "_")$V2)
  meta$rep <- as.factor(gsub(".*_rep", "", row.names(meta)))
  
  return(meta)
}


# Create boxplots of total intensity, colored by group of choice

intensity.box.group <- function(map.results, meta.data, group.name, take.log = F, data.type, color.vector = NULL, text.size = 1){
  # map.results = list of results with descriptive sample names mapped
  # meta.data = data frame of sample groupings (using create.meta())
  # group.name = desired grouping column name (from meta.data) with which to color boxes
  # take.log = logical indicating whether to take log2 of data or not
  # data.type = character variable of input data type (lipid/metab raw, imputed, normalized) for title
  # color.vector = pre-defined color vector (optional)
  # text.size = size of axis labels
  
  # ensure meta-data and metabolite data have samples in the same order
  
  if(all.equal(row.names(meta.data), row.names(map.results[[1]])) == F){
    print("meta.data and metab.meta can't be coerced into same group order - check sample names in both")
  }
  
  for(i in 1:length(map.results)){
    
    metab.data <- as.data.frame(t(map.results[[i]]))
    
    # order by group levels 
    metab.data <- metab.data[,order(meta.data$group_factor_num)]
    meta.data.plot <- meta.data[order(meta.data$group_factor_num),]
    
    # take log2 of metab data, if provided
    
    if(take.log == F){
      metab.data <- metab.data
    }else{
      metab.data <- log2(metab.data)
    }
    
    # set colors 
    
    if(is.null(color.vector)){
      
      group.values <- meta.data.plot[,c(group.name)]
      
      if(length(unique(group.values)) > 2){
        colors <- brewer.pal(n = length(unique(group.values)), "Accent")
      }else{
        colors <- c("lightgrey", "blue")
      }
      
      names(colors) <- c(unique(group.values))
      
    }else{
      
      group.values <- meta.data.plot[,c(group.name)]
      
      colors <- color.vector
      
    }

    
    # create boxplots COLORS AND AXIS LABELS DON'T MATCH
     
    boxplot(metab.data[,order(meta.data.plot$group_factor_num)], col = colors[group.values], cex.axis = text.size, las = 2, ylim = c(0, 35))
    #legend("top",
           #legend = c(names(colors)),
           #fill = colors,
           #horiz = T, 
           #cex = 0.75)
    #title(paste(names(map.results)[[i]], data.type, "data colored by", group.name))  
    
  }
  
}


# Prep UpSet data

upset.prep <- function(comp.class, column){
  # comp.class = list of results with compounds found and super classes mapped
  # column = character string of name of column to use for matching (usually Name or MF)
  
  # Mark whether metabolites are found in each lab
  ## list unique metabolite names
  
  unlist.matrix <- data.frame() # how to make matrix out of list of lists?
  
  i <- 1
  while(i < length(comp.class) + 1){
    print(i)
    if(i == 1){
      
      unlist.matrix <- comp.class[[i]]
      
    }else{
      
      unlist.matrix <- rbind(unlist.matrix, comp.class[[i]])
      
    }
    
    i <- i +1
    
  }
  
  metabs <- unlist.matrix[,c(column)]
  
  metabs <- data.frame(metab = unique(metabs))
  metabs[is.na(metabs)] <- "unknown"
  
  ## add columns for labs
  for(i in 1:length(comp.class)){
    
    metabs[,c(names(comp.class)[[i]])] <- NA
    
  }
  
  ## mark metabolite presence in lab
  for(k in 1:nrow(metabs)){
    for(j in 1:length(comp.class)){
      
      if(has_element(comp.class[[j]][,c(column)], metabs$metab[k]) == T){
        metabs[k,c(names(comp.class)[[j]])] <- 1
      }
      if(has_element(comp.class[[j]][,c(column)], metabs$metab[k]) == F){
        metabs[k,c(names(comp.class)[[j]])] <- 0
      } 
    }
  }
  
  ## add column for super class
  
  classes <- data.frame()
  
  for(i in 1:length(comp.class)){
    
    lab.class <- data.frame(metab = comp.class[[i]][,c(column)], class = comp.class[[i]]$ramp.class)
    
    classes <- rbind(classes, lab.class)
    
  }
  
  classes <- classes[which(duplicated(classes$metab)==F),]
  classes[is.na(classes)] <- "unknown"
  
  for(i in 1:nrow(metabs)){
    for(j in 1:nrow(classes)){
      
      if(metabs$metab[i] == classes$metab[j]){
        metabs$class[i] <- classes$class[j]
      }
      
    }
  }
  
  return(metabs)
}

# Correlation of duplicate samples

cor.dupes <- function(map.results){
  # map.results = list of results with descriptive sample names mapped
  
  all.dups <- list()
  
  for(i in 1:length(map.results)){
    
    # set sample names
    samples <- gsub("_rep.*", "", row.names(map.results[[i]]))
    
    # create matrix for correlation results
    all.dups[[i]] <- data.frame(matrix(nrow = length(samples), ncol = 2))
    colnames(all.dups[[i]]) <- c("sample", "correlation")
    all.dups[[i]]$sample <- samples
    
    # calculate correlations and add to matrix
    
    for(j in 1:length(samples)){
      
      pair <- t(map.results[[i]][which(grepl(samples[j], row.names(map.results[[i]]))==T),])
      
      cor.pair <- cor.test(pair[,1], pair[,2], method = "pearson")
      
      all.dups[[i]]$correlation[j] <- cor.pair$estimate
      
    }
    
    all.dups[[i]] <- all.dups[[i]][-which(duplicated(all.dups[[i]]$sample)),]
    
  }
  
  names(all.dups) <- names(map.results)
  
  return(all.dups)
}

# Correlation of random pairs

cor.rando <- function(map.results, n){
  # map.results = list of results with descriptive sample names mapped
  # n = number of random pairs
  
  random.dups <- list()

  # create list of actual replicate indexes

  samples <- row.names(map.results[[1]])
  samples.norep <- gsub("_rep.*", "", row.names(map.results[[1]]))
  
  rep.index <- data.frame(matrix(nrow = length(unique(samples.norep)), ncol = 3))
  colnames(rep.index) <- c("rep1", "rep2", "pasted")
  row.names(rep.index) <- unique(samples.norep)
  
  for(i in 1:nrow(rep.index)){
    for(j in 1:length(samples)){
      
      if(grepl(row.names(rep.index)[i], samples[j]) &
         grepl("rep1", samples[j])){
        rep.index[i, c("rep1")] <- which(samples == paste0(row.names(rep.index)[i], "_rep1"))
      }

      if(grepl(row.names(rep.index)[i], samples[j]) &
         grepl("rep2", samples[j])){
        rep.index[i, c("rep2")] <- which(samples == paste0(row.names(rep.index)[i], "_rep2"))
      }
    }
    
    rep.index[i, c("pasted")] <- paste0(rep.index[i,c("rep1")], "|", rep.index[i,c("rep2")])
  }
  
  # create list of all possible pair indexes
  
  pair.index <- as.data.frame(t(combn(nrow(map.results[[1]]), 2)))
  colnames(pair.index) <- c("pair1", "pair2")
  pair.index$pasted <- paste0(pair.index[,c("pair1")], "|", pair.index[,c("pair2")])
  
  # remove pair indexes that reflect actual replicates
  
  fin.index <- pair.index[-which(pair.index$pasted %in% rep.index[,c("pasted")]),]
  
  # sample "n" random pairs 
  
  sample.index <- as.numeric(sample(nrow(fin.index), n, replace = F))
  to.analyze <- fin.index[c(sample.index),]
  
  for(k in 1:length(map.results)){
    
    # create matrix for correlation results
    random.dups[[k]] <- data.frame(matrix(nrow = n, ncol = 2))
    colnames(random.dups[[k]]) <- c("pair", "correlation")
    
    # calculate correlations
    
    for(m in 1:nrow(to.analyze)){
      
      samp1 <- as.numeric(to.analyze[m,c("pair1")])
      samp2 <- as.numeric(to.analyze[m,c("pair2")])
      
      cor.pair <- cor(t(map.results[[k]][c(samp1),]),
                           t(map.results[[k]][c(samp2),]),
                           method = "pearson")
      
      random.dups[[k]]$correlation[m] <- diag(cor.pair)
      random.dups[[k]]$pair[m] <- paste0(row.names(map.results[[k]])[samp1],
                                         "|",
                                         row.names(map.results[[k]])[samp2])
      #plot(t(map.results[[k]][c(samp1),]),
                      #t(map.results[[k]][c(samp2),]))
      #text(t(map.results[[k]][c(samp1),]),
           #t(map.results[[k]][c(samp2),]),
           #labels = colnames(map.results[[k]]))
      
    }
  
  }
  
  names(random.dups) <- names(map.results)
  
  return(random.dups)
}

# Plot histograms of duplicates and random pair correlations

cor.hist <- function(dup.cor, rando.cor, data.type, breaks){
  # dup.cor = results of correlations between duplicate pairs
  # rando.cor = results of correlations between random pairs
  # data.type = character description of type of data input
  # breaks = bin size used for histogram
  
  if(length(dup.cor) != length(rando.cor)){
    message("check your inputs - not the same length")
  }
  
  for(i in 1:length(dup.cor)){
    
    dup.hist <- hist(dup.cor[[i]]$correlation, breaks = breaks, plot = F)
    rand.hist <- hist(rando.cor[[i]]$correlation, breaks = breaks, plot = F)
    
    rand.plot <- plot(rand.hist, xlim = c(0, 1), col= "lightgrey", main = paste(names(dup.cor)[[i]], data.type, xlab = "Correlation"))
    dup.plot <- plot(dup.hist, xlim = c(0, 1), col= "red", add = T, main = "", xlab = "Correlation")

  }
  
}

# Calculate ratios of intensities for samples of the same group but different concentration

ratio.calc <- function(map.results, ratio.map){
  # map.results = list of results with sample names mapped
  # ratio.map = external csv file mapping samples to expected ratios
  
  conc.pairs <- list()
  
  for(i in 1:length(map.results)){
    
    conc.pairs[[i]] <- data.frame(matrix(ncol = 6))
    colnames(conc.pairs[[i]]) <- c("group", "rep", "ratio_label", "expected_ratio", "metab", "real_ratio")
    
    for(k in 1:nrow(ratio.map)){
      
      pair <- map.results[[i]][c(ratio.map$sample_name_1[k], ratio.map$sample_name_2[k]),]
      
      for(j in 1:ncol(map.results[[i]])){
        
        conc.pairs[[i]] <- add_row(group = ratio.map$group[k],
                                        rep = as.factor(ratio.map$rep[k]),
                                        ratio_label = paste(row.names(pair)[1], "/", row.names(pair)[2]),
                                        expected_ratio = as.factor(ratio.map$ratio[k]),
                                        metab = colnames(map.results[[i]])[j],
                                        real_ratio = (pair[1,j] / pair[2,j]), 
                                        conc.pairs[[i]])
        
      }
    }
    
    conc.pairs[[i]] <- conc.pairs[[i]][-which(is.na(conc.pairs[[i]]$group)),]
    conc.pairs[[i]]$expected_ratio <- factor(conc.pairs[[i]]$expected_ratio, c("1.5", "2", "3"))
  }
  
  names(conc.pairs) <- names(map.results)
  
  return(conc.pairs)
}

# Make boxplots of real vs. expected ratios

ratio.box <- function(pair.results, filter.num){
  # pair.results = list of results from ratio.calc()
  # filter.num = cutoff for removing outlier ratios
  
  box.list <- list()
  
  for(i in 1:length(pair.results)){
    
    print(names(pair.results)[[i]])
    
    filtered.conc <- pair.results[[i]][-which(pair.results[[i]]$real_ratio > filter.num),]
    
    print(paste("filtered features with ratio >", filter.num, "=", length(which(pair.results[[i]]$real_ratio > filter.num)),
                "out of total:", nrow(pair.results[[i]])))
    
    boxes <- ggplot(data = filtered.conc,
                    aes(x = ratio_label,
                        y = real_ratio,
                        fill = group,
                        pattern = rep)) +
      geom_boxplot(position = position_dodge()) +
      geom_hline(yintercept = 1.5, linetype = "dashed") +
      geom_hline(yintercept = 2, linetype = "dashed") +
      geom_hline(yintercept = 3, linetype = "dashed") +
      scale_pattern_manual(values = c(one = "stripe", two = "none")) +
      #theme_classic() +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~expected_ratio, scales = "free") +
      ggtitle(paste(names(pair.results)[i], length(which(pair.results[[i]]$real_ratio > filter.num)), "/", 
                    nrow(pair.results[[i]]), "ratios removed for having ratio >", filter.num)) 
    
    if(i == 1){
      box.list[[i]] <- boxes
    }else{
      box.list[[length(box.list) + 1]] <- boxes
    }
    
  }
  
  arranged.boxes <- ggarrange(plotlist = box.list, ncol = 3, nrow = 3)
  return(arranged.boxes)
  
}

# Concatenate data into one matrix, transpose, and include column for lab name

concat.data <- function(map.results, master.map){
  # results.list = list of metabolite matrices from each lab with sample names mapped
  # master.map = list of descriptive sample names
  
  # prep concatenated dataframe
  
  concat.results <- data.frame()
  concat.meta <- c()
  
  # prep objects to check row orders of all matrices in list
  
  results.names <- data.frame(matrix(ncol = length(map.results), nrow = nrow(map.results[[1]])))
  colnames(results.names) <- names(map.results)
  
  perfect.names <- unique(master.map$description)
  perfect.names <- perfect.names[order(perfect.names)]
  
  for(i in 1:length(map.results)){
    
    map.results[[i]] <- map.results[[i]][order(row.names(map.results[[i]])),]
    
    results.names[,i] <- row.names(map.results[[i]]) 
    
  }
  
  # if all of the row orders are the same, add them to the concatenated dataframe
  
  if(any(apply(results.names, 2, function(x) x == perfect.names)) == F){
    print("sample names don't match - check that all samples are included and named properly")
  }else{
    
    for(i in 1:length(map.results)){
      
      results.i <- map.results[[i]]
      meta.i <- rep.int(names(map.results)[[i]], ncol(map.results[[i]]))
      
      if(i == 1){
        
        concat.results <- results.i
        concat.meta <- meta.i
        
      }else{
        
        concat.results <- cbind(concat.results, results.i)
        concat.meta <- c(concat.meta, meta.i)
        
      }
      
    }
    
  }
  
  concat.results <- t(concat.results)
  
  concat.results <- as.data.frame(cbind(concat.meta, concat.results))
  colnames(concat.results)[1] <- "lab"
  
  return(concat.results)
  
}
  

# Interactive PCA of lab clustering

make.pca <- function(metab, meta, class, group.shape.name = NULL,
                     centered = FALSE, scaled = FALSE, do.log2 = FALSE, z.axis = F,
                     data.type = NULL, subset.labels = NULL,
                     point.size = 5, text.size = 10, color.vector = NULL){
  # metab = metabolite data to plot with variables in columns and samples in rows
  # meta = meta data corresponding to metab
  # class = variable name from meta used for coloring scores plot
  # group.shape.name = variable name from meta used to make point shapes
  # data.type = character string describing whether data is normalized, corrected, or raw
  # subset.labels = vector of sample types to subset plots into *MUST match variable names in class column
    ## (e.g., if you only want to look at a few groups from the whole PCA)
  # point.size = size of plot points
  # text.size = size of plot text
  # color.vector = optional color vector to supply--must match the class assigned in "class" argument
  
  # Set data frames for plot
  
  if(do.log2 == T){
    pca_data <- log2(metab)
  }else{
    pca_data <- metab  
  }
  
  pca_data <- prcomp(pca_data, center = centered, scale = scaled)
  exp_var_pca <- round(((pca_data$sdev^2)/sum(pca_data$sdev^2)*100)[1:3],2)
  
  if(z.axis == F){
    pca_plotdata <- data.frame(PC1=pca_data$x[, 1], PC2=pca_data$x[, 2], class=meta[,c(class)])   
  }else{
    pca_plotdata <- data.frame(PC1=pca_data$x[, 1], PC2=pca_data$x[, 2], PC3=pca_data$x[,3], class=meta[,c(class)])
  }
  
  if(!is.null(group.shape.name)){
    pca_plotdata$shape <- as.factor(meta[,c(group.shape.name)])
  }
  
  # Set colors
  
  if(is.null(color.vector)){
    
    if(class(pca_plotdata[,c("class")]) == "factor" |
       class(pca_plotdata[,c("class")]) == "character"){
      
      grafify_kelly_pal <- graf_palettes$kelly
      
      my_colors <- grafify_kelly_pal[1:length(unique(meta[,c(class)]))]
      names(my_colors) <- unique(meta[,c(class)])
      
    }
    
  }else{
    
    my_colors <- color.vector
    
  }
  
  if(is.null(subset.labels)){
    
    subset.labels <- meta[,c(class)]

    }
  
  # Set shapes
  
  if(!is.null(group.shape.name)){
    shapes <- c(1:length(unique(meta[,c(group.shape.name)])))
    names(shapes) <- unique(meta[,c(group.shape.name)])
    pca_plotdata$shape <- NA 
      for(i in 1:nrow(pca_plotdata)){
        for(j in 1:length(shapes)){
         if(grepl(paste0("_", names(shapes)[j], "_"), row.names(pca_plotdata)[i])){
           pca_plotdata$shape[i] <- shapes[j]
         } 
        }
      }
    
    pca_plotdata$shape <- as.factor(pca_plotdata$shape)
    
  }else{
    shapes <- NA
    pca_plotdata$shape <- as.factor(1)
  }
  
  
  
  if(z.axis == FALSE){
    
    # 2D  using ggplot2
    
    pca_plot <- ggplot(pca_plotdata[which(pca_plotdata$class %in% subset.labels),], 
      aes(x=PC1, y=PC2, col=class)) + 
      geom_point(size = point.size, aes(shape = shape, color = class)) +
      #scale_shape_manual(shapes) +
      guides(size = "none") +
      theme_classic() +
      theme(text = element_text(size = text.size)) +
      xlab(paste0("PC1 (", exp_var_pca[1] ," %)")) +
      ylab(paste0("PC2 (", exp_var_pca[2] ," %)")) +
      ggtitle(data.type)
    
    if(class(pca_plotdata[,c("class")]) == "factor" | 
    class(pca_plotdata[,c("class")]) == "character"){
    
        pca_plot_fin <- pca_plot + scale_color_manual(values=my_colors)
    
      }
    
    if(class(pca_plotdata[,c("class")]) == "integer"){
    
        pca_plot_fin <- pca_plot + scale_color_gradient()
    
      }
    
    return(pca_plot_fin)
    
  }else{
    
    # 3D using plotly

      pca_plotly <- plot_ly(pca_plotdata, x = ~PC1, y = ~PC2, z = ~PC3, 
                            color = ~class, mode = "markers", symbol = ~shapes, 
                            colors = group.colors, symbols = shapes,
                            marker = list(point.size))
      pca_plotly <- pca_plotly %>% add_markers()
      pca_plotly <- pca_plotly %>% layout(scene = list(xaxis = list(title = paste0("PC1 (", exp_var_pca[1] ," %)"), font = list(size = text.size)),
                                                       yaxis = list(title = paste0("PC2 (", exp_var_pca[2] ," %)"), font = list(size = text.size)),
                                                       zaxis = list(title = paste0("PC3 (", exp_var_pca[3] ," %)"), font = list(size = text.size))))
    
      return(pca_plotly)
      
      }
    
}

# Evaluate PCA based on regression with study design

eval.pca <- function(metab, meta, class,
                     do.log2 = T, centered = FALSE, scaled = FALSE){
  # metab = metabolite data to plot with variables in columns and samples in rows
  # meta = meta data corresponding to metab
  # class = variable name from meta used for coloring scores plot
  
  
    # Set data frames for plot
    
    if(do.log2 == T){
      pca_data <- log2(metab)
    }else{
      pca_data <- metab  
    }
    
    pca_data <- prcomp(pca_data, center = centered, scale = scaled)
    exp_var_pca <- round(((pca_data$sdev^2)/sum(pca_data$sdev^2)*100)[1:3],2)
    
    pca_plotdata <- data.frame(PC1=pca_data$x[, 1], PC2=pca_data$x[, 2], PC3=pca_data$x[,3], class=meta[,c(class)])
    
  
    # Calculate regression as a function of the class
    
    regression.vals <- data.frame(matrix(nrow = ncol(pca_plotdata) - 1, ncol = 1))
    row.names(regression.vals) <- colnames(pca_plotdata)[1:ncol(pca_plotdata) - 1]
    colnames(regression.vals) <- "P-value"
    
    for(i in 1:(ncol(pca_plotdata) - 1)){
      
      model <- anova(lm(pca_plotdata[,i] ~ pca_plotdata$class))
      
      regression.vals[i,c("P-value")] <- formatC(model$`Pr(>F)`[1], format = "e", digits = 2)
      
    }
  
  return(regression.vals)
  
}
  
# Find loadings of features for PCA

feat.pca <- function(metab, meta, class,
                     centered = FALSE, scaled = FALSE, do.log2 = FALSE, z.axis = F){
  # metab = metabolite data to plot with variables in columns and samples in rows
  # meta = meta data corresponding to metab
  # class = variable name from meta used for coloring scores plot
  # z.axis = logical to include 3 PCs
  
  # Set data frames for plot
  
    if(do.log2 == T){
      pca_data <- log2(metab)
    }else{
      pca_data <- metab  
    }
  
  # Calculate PCs
    pca_data <- prcomp(pca_data, center = centered, scale = scaled)
    exp_var_pca <- round(((pca_data$sdev^2)/sum(pca_data$sdev^2)*100)[1:3],2)
    
    
  # Extract loadings
    if(z.axis == F){
      pca_loadings <- data.frame(PC1=pca_data$rotation[,c("PC1")], PC2=pca_data$rotation[,c("PC2")])   
    }else{
      pca_loadings <- data.frame(PC1=pca_data$rotation[,c("PC1")], PC2=pca_data$rotation[,c("PC2")], PC3=pca_data$rotation[,c("PC3")])
    }
    
    return(pca_loadings)
    
}

# Find scores for PCA

score.pca <- function(metab, meta, class,
                     centered = FALSE, scaled = FALSE, do.log2 = FALSE, z.axis = F){
  # metab = metabolite data to plot with variables in columns and samples in rows
  # meta = meta data corresponding to metab
  # class = variable name from meta used for coloring scores plot
  # z.axis = logical to include 3 PCs
  
  # Set data frames for plot
  
  if(do.log2 == T){
    pca_data <- log2(metab)
  }else{
    pca_data <- metab  
  }
  
  # Calculate PCs
  pca_data <- prcomp(pca_data, center = centered, scale = scaled)
  exp_var_pca <- round(((pca_data$sdev^2)/sum(pca_data$sdev^2)*100)[1:3],2)
  
  
  # Extract loadings
  if(z.axis == F){
    pca_scores <- data.frame(PC1=pca_data$x[,c("PC1")], PC2=pca_data$x[,c("PC2")])   
  }else{
    pca_scores <- data.frame(PC1=pca_data$x[,c("PC1")], PC2=pca_data$x[,c("PC2")], PC3=pca_data$x[,c("PC3")])
  }
  
  return(pca_scores)
  
}

# Apply MSTUS

mstus <- function(data, low.bound = 0.2, high.bound = 0.95){
  # data = metabolite abundance data with samples in rows and metabolites in columns
  # low.bound = lower quantile for MSTUS as decimal
  # high.bound = upper quantile for MSTUS as decimal
  
  # calculate low and high quantiles for each sample and find metabolites outside bounds
  outliers <- function(x) {
    
    lohi <- t(apply(x,1,quantile,c(low.bound, high.bound)))
    
  }
  
  lohi <- outliers(data)
  
  # mark metabolites outside bounds
  lobinfilt <- hibinfilt <- matrix(0,nrow=nrow(data), ncol=ncol(data))
  
  for ( i in 1:nrow(data)) {
    
    lobinfilt[ i, which( data[ i , ] < lohi[ i, 1 ] ) ] = 1
    hibinfilt[ i, which( data[ i , ] > lohi[ i, 2 ] ) ] = 1
    
  }
  
  sumlo <- colSums(lobinfilt)
  sumhi <- colSums(hibinfilt)  
  
  # find indexes of features with intensities less than the lower bound in 80%+ of samples or greater than the higher bound in < 10% of samples
  bads <- c( which( sumlo > ( 0.8 * nrow(data) ) ), which( sumhi > 0 & sumhi < ( 0.1 * nrow(data) ) ) )
  length(bads) 
  length(unique(bads)) 
  
  # calculate scale factor without features outside bounds
  scalefactors <- rowSums(data[,-bads])
  
  # normalize according to scale factors
  normmetab <- data / scalefactors
  
  # add column and row names
  colnames(normmetab) <- colnames(data)
  row.names(normmetab) <- row.names(data)
  
  return(normmetab)

}
  
  
  
  
  
  
  
  
  
  
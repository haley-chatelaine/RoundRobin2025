# Create functions
print("create plotting functions")

## i. m/z ppm ranges *use mode-specific list of ionized m/z's (typically M+H or M-H) 

mz.range <- function(ion.list, ppm, rt = NULL){
  # ion.list = matrix of metabolites and theoretical ionization m/z for each *mode-specific*
  # ppm = desired ppm tolerance for m/z range
  # rt = desired RT tolerance
  
  ## Extract m/z ranges
  ion.list$lower.mz <- NA
  ion.list$upper.mz <- NA
  
  for(i in 1:nrow(ion.list)){
    
    ## calculate lower m/z bound
    ion.list$lower.mz[i] <- (-(as.numeric(ion.list$mz[i])) * (ppm/1000000)) + ion.list$mz[i]
    
    ## calculate upper m/z bound
    ion.list$upper.mz[i] <- (ion.list$mz[i] * (ppm/1000000)) + ion.list$mz[i]
    
  }
  
  ## Extract RT ranges, if available
  
  ion.list$lower.rt <- NA
  ion.list$upper.rt <- NA
  
  if(!is.null(rt)){
    
    for(i in 1:nrow(ion.list)){
      
      ## calculate lower rt bound
      ion.list$lower.rt[i] <- as.numeric(ion.list$rt[i]) - 0.5*rt
      
      if(ion.list$lower.rt[i] < 0){ ion.list$lower.rt[i] <- 0 }
      
      ## calculate upper m/z bound
      ion.list$upper.rt[i] <- as.numeric(ion.list$rt[i]) + 0.5*rt
      
    }
    
  }else{
    
    ion.list$lower.rt[i] <- NA
    
    ## calculate upper m/z bound
    ion.list$upper.rt[i] <- NA
    
  }
  
  ## retain only metabolite name and ranges
  part.list <- subset(ion.list, select = c(Metabolite, lower.mz, upper.mz, lower.rt, upper.rt))
  
  ## make into a list of lists
  fin.list <- lapply(split(part.list, part.list$Metabolite), as.list)
  
  return(fin.list)
}

## ii. random m/z ppm ranges 

mz.random.range <- function(data, size, ppm){
  # data = raw data (OnDiskMSnExp) object
  # size = number of random m/z's to generate
  # ppm = desired ppm tolerance for m/z range
  
  ion.list <- data.frame(matrix(nrow = size, ncol = 4))
  colnames(ion.list) <- c("random.number", "mz", "lower.mz", "upper.mz")
  
  # Randomly select m/z's of features in data of size supplied
  
  set.seed(3)
  data.index <- sample(1:length(data@featureData@data[["basePeakMZ"]]), size)
  
  for(i in 1:nrow(ion.list)){
    
    ion.list$mz[i] <- data@featureData@data[["basePeakMZ"]][[(data.index[i])]]
    
  }
  
  # Create lower/upper bounds, based on ppm supplied
  
  ion.list$lower.mz <- NA
  ion.list$upper.mz <- NA
  
  for(i in 1:nrow(ion.list)){
    
    ## calculate lower m/z bound
    ion.list$lower.mz[i] <- (-(as.numeric(ion.list$mz[i])) * (ppm/1000000)) + ion.list$mz[i]
    
    ## calculate upper m/z bound
    ion.list$upper.mz[i] <- (ion.list$mz[i] * (ppm/1000000)) + ion.list$mz[i]
    
  }
  
  # Name metabolites based on m/z
  
  ion.list$random.number <- paste0("feat.", seq(1, size, 1), "_", round(ion.list$mz, 2))
  
  # Add empty RT column to not confuse the EICs
  
  ion.list$lower.rt <- NA
  ion.list$upper.rt <- NA
  
  ## retain only metabolite name and ranges
  part.list <- subset(ion.list, select = c(random.number, mz, lower.mz, upper.mz, lower.rt, upper.rt))
  
  
  ## make into a list of lists
  
  fin.list <- lapply(split(part.list, part.list$random.number), as.list)
  
  return(fin.list)
}

## iii. plot EICs

plot.eics <- function(data, ion.list, group_colors = NULL, step.name = NULL, file.suffix){
  # data = input MSnExp or XCMSnExp
  # ion.list = mass ranges derived from mz.range
  # group_colors = vector of colors for groups, if desired
  # step.name = name of processing step for inclusion on EIC title
  # file.suffix = file name suffix for saved pdf (can be set in a config file)
  
  ## create empty list
  all.chr <- list()
  
  for(i in 1:length(ion.list)){
    
    print(i)
    
    if(is.na(ion.list[[i]]$lower.rt) == TRUE){
      
      ## create chromatograms for each metabolite of interest
      chr.i <- chromatogram(data, 
                            mz = c(ion.list[[i]]$lower.mz, ion.list[[i]]$upper.mz))
      
    }
    
    if(is.na(ion.list[[i]]$lower.rt) == FALSE){
      
      chr.i <- chromatogram(data, 
                            mz = c(ion.list[[i]]$lower.mz, ion.list[[i]]$upper.mz), 
                            rt = c(ion.list[[i]]$lower.rt, ion.list[[i]]$upper.rt))
      
    }
    
    ## create list of total chromatogram objects
    
    ### IF it's the first file, add it to the list
    if(i == 1){
      all.chr[[i]] <- chr.i
    }else{
      
      ### ELSE append file to list
      
      all.chr[[length(all.chr) + 1]] <- chr.i 
    }
  }
  
  ## set metabolite names
  names(all.chr) <- names(ion.list)
  
  ## plot EICs of chromatograms in list
  fun.plot <- function(x, y) plot(x, col = group_colors[x$group], main = paste(y, step.name, file.suffix), cex.main = 0.5)
  
  mapply(fun.plot, all.chr, names(all.chr))
  
}
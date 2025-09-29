# These functions are for qualitative comparisons of CF and XCMS results

# i. Function for calculating overlapping ranges
    ## logic: if range a + range b > the min and max of either values

overlaps <- function(a.min, a.max, b.min, b.max){
  # a.min = minimum value of list a
  # a.max = maximum value of list a
  # b.min = minimum value of list b
  # b.max = maximum value of list b  
  
  max(as.numeric(a.max), as.numeric(b.max)) - min(as.numeric(a.min), as.numeric(b.min)) < (as.numeric(a.max) - as.numeric(a.min)) + (as.numeric(b.max) - as.numeric(b.min))
  
}

# ii. Create list of metabolites found on both platforms

id_compare <- function(cf_metabs, metab_meta, ppm, rt){
  # cf_metabs = identified compounds sheet from ClusterFinder results
  # metab_meta = metabolite meta-data (mz, rt, mode) from open-source preprocess
  # ppm = desired ppm range
  # rt = desired rt range
  
  matched.ids <- data.frame(matrix(ncol = 12))
  colnames(matched.ids) <- c("name", "cf_mz", "cf_rt", "xcms_mz_med", "xcms_mz_max", "xcms_mz_min", "xcms_rt_med", "xcms_rt_max", "xcms_rt_min", "cf_mode", "xcms_mode", "C12_num_c")
  
      # create matrix of anticipated values from CF, based on ppm and rt windows
      cf_ranges <- data.frame(matrix(ncol = 6, nrow = nrow(cf_metabs)))
      colnames(cf_ranges) <- c("name", "mz_min", "mz_max", "rt_min", "rt_max", "cf_mode")
      
      ## ppm functions
      ppm_max <- function (x) as.numeric(x) + as.numeric(x) * (ppm/1000000)
      ppm_min <- function (x) as.numeric(x) - as.numeric(x) * (ppm/1000000)
      
      ## rt functions
      rt_max <- function (x) x + 0.5*rt
      rt_min <- function (x) x - 0.5*rt
      
      ## calculate and populate CF matrix
      cf_ranges$name <- cf_metabs$Name
      cf_ranges$mz_min <- as.numeric(lapply(cf_metabs$BinC12mz, ppm_min))
      cf_ranges$mz_max <- as.numeric(lapply(cf_metabs$BinC12mz, ppm_max))
      cf_ranges$rt_min <- as.numeric(lapply(cf_metabs$BinRTsec, rt_min))
      cf_ranges$rt_max <- as.numeric(lapply(cf_metabs$BinRTsec, rt_max))
      cf_ranges$cf_mode <- cf_metabs$Mode
      
      # create columns for ppm and rt ranges around XCMS medians
      metab_meta$min.mz <- as.numeric(lapply(metab_meta$med.mz, ppm_min))
      metab_meta$max.mz <- as.numeric(lapply(metab_meta$med.mz, ppm_max))
      metab_meta$min.rt <- as.numeric(lapply(metab_meta$med.rt, rt_min))
      metab_meta$max.rt <- as.numeric(lapply(metab_meta$med.rt, rt_max))
  
  # determine if metabolites from XCMS fall within ranges of CF (for C12 annotations)
  
  for(i in 1:nrow(cf_ranges)){
    for(j in 1:nrow(metab_meta)){
      
      if(
        # m/z
        overlaps(cf_ranges$mz_min[i], cf_ranges$mz_max[i], metab_meta$min.mz[j], metab_meta$max.mz[j]) == TRUE
        &&
        # rt
        overlaps(cf_ranges$rt_min[i], cf_ranges$rt_max[i], metab_meta$min.rt[j], metab_meta$max.rt[j]) == TRUE
        &&
        # mode
        cf_ranges$cf_mode[i] == metab_meta$mode[j])
        
      {matched.ids <- add_row(matched.ids, name = cf_metabs$Name[i], 
                              cf_mz = cf_metabs$BinC12mz[i],  
                              cf_rt = cf_metabs$BinRTsec[i], 
                              xcms_mz_med = metab_meta$med.mz[j],
                              xcms_mz_max = metab_meta$max.mz[j],
                              xcms_mz_min = metab_meta$min.mz[j],
                              xcms_rt_med = metab_meta$med.rt[j],
                              xcms_rt_max = metab_meta$max.rt[j],
                              xcms_rt_min = metab_meta$min.rt[j],
                              cf_mode = cf_metabs$Mode[i], 
                              xcms_mode = metab_meta$mode[j],
                              C12_num_c = cf_metabs$NumC[i])}else{next}
      
    }
  }
  
  matched.ids <- matched.ids[-c(which(is.na(matched.ids$name)==T)),]
  
  # determine if C13 cognate is likely in XCMS data
  
  for(i in 1:nrow(metab_meta)){
    for(j in 1:nrow(matched.ids)){
      
      # calculate C13 mass
      
      c13.j <- matched.ids$cf_mz[j] + matched.ids$C12_num_c[j]
      
      # create ppm window around C13 mass
      
      c13.j.mzmin <- ppm_min(c13.j)
      c13.j.mzmax <- ppm_max(c13.j)
      
      # create RT window
      
      c13.j.rtmin <- rt_min(matched.ids$cf_rt[j])
      c13.j.rtmax <- rt_max(matched.ids$cf_rt[j])
      
      # search metab_meta for C13 ppm and rt window and record name as "name | C13 pair"
      
      if(
        # m/z
        overlaps(c13.j.mzmin, c13.j.mzmax, metab_meta$min.mz[i], metab_meta$max.mz[i]) == TRUE
        &&
        # rt
        overlaps(c13.j.rtmin, c13.j.rtmax, metab_meta$min.rt[i], metab_meta$max.rt[i]) == TRUE)
        
      {matched.ids <- add_row(matched.ids, name = paste0(matched.ids$name[j], " | C13 cognate"), 
                              cf_mz = matched.ids$cf_mz[j],  
                              cf_rt = matched.ids$cf_rt[j], 
                              xcms_mz_med = metab_meta$med.mz[i],
                              xcms_mz_min = metab_meta$min.mz[i], 
                              xcms_mz_max = metab_meta$max.mz[i], 
                              xcms_rt_med = metab_meta$med.rt[i], 
                              xcms_rt_min = metab_meta$min.rt[i], 
                              xcms_rt_max = metab_meta$max.rt[i], 
                              cf_mode = matched.ids$cf_mode[j], 
                              xcms_mode = metab_meta$mode[i],
                              C12_num_c = matched.ids$C12_num_c[j])}else{next}
      
    }
    
  }
  
  # find duplicate CF features and only retain the XCMS feature that is closest in m/z and rt 
  
    # make a list of CF features with duplicates
  
    
  
  return(matched.ids)
  
}

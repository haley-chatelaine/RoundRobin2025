# This script contains functions for conducting RF analyses based on SAMPLE GROUP in each lab

# remove metabolites with low CV

cv.filter <- function(results.list, cv.cutoff){
  # results.list = list of metabolite abundances from each lab
  # cv.cutoff = cutoff, given as decimal, for cv quartile
  
  
  filter.list <- list()
  
  for(i in 1:length(results.list)){
    
    sds <- apply(results.list[[i]], 2, sd)
    means <- colMeans(results.list[[i]])
    cvs <- sds/means
    
    cv.percentile <- quantile(cvs, cv.cutoff)
    
    filter.list[[i]] <- results.list[[i]][,which(cvs > cv.percentile)]
    
  }
  
  names(filter.list) <- names(results.list)
  
  return(filter.list)
  
}

# remove metabolites with m/z > number

mz.filter <- function(results.list, comp.found, mz.cutoff){
  # results.list = list of metabolite abundances from each lab
  # comp.found = list of metabolite meta-data from each lab
  # mz.cutoff = m/z below which to keep
  
  filter.list <- list()
  
  for(i in 1:length(results.list)){
    
    filter.list[[i]] <- results.list[[i]][,which(comp.found[[i]]$BinC12mz < mz.cutoff)]
    
  }
  
  names(filter.list) <- names(results.list)
  
  return(filter.list)
  
}

# Split training and testing data and add group column from meta.data

split.data <- function(results.list, meta.data, split = 0.75){
  # results.list = list of metabolite abundances from each lab
  # meta.data = meta.data of sample groups
  # split = fraction to split to train and test
  
  split.list <- list()
  
  for(i in 1:length(results.list)){
    
    # combine meta data and results
    #input.data <- cbind(group = meta.data$group, results.list[[i]])
    input.data <- results.list[[i]]
    
    # Step 1: Get row numbers for the training data
    input.data.indexes <- data.frame(index = 1:nrow(input.data), group = meta.data$group)
    
        ## Ensure 1 sample from each group is included in training
        trainRowNumbers <- input.data.indexes %>% 
                            group_by(group) %>% 
                            sample_n(size = 1) %>% 
                            ungroup() 
        
        remainingRowNumbers <- input.data.indexes[-which(input.data.indexes$index %in% trainRowNumbers$index),]
        
        ## Randomly select the remaining needed training samples from the remaining samples
        trainRowNumbers <- rbind(trainRowNumbers,
                                 input.data.indexes[sample(remainingRowNumbers$index, size = floor(split * nrow(input.data)) - 7),])
    
    # Step 2: Create the training  dataset
    trainData <- input.data[trainRowNumbers$index,]
    trainMeta <- meta.data[trainRowNumbers$index,]
    
    # Step 3: Create the test dataset
    testData <- input.data[-trainRowNumbers$index,]
    testMeta <- meta.data[-trainRowNumbers$index,]
    
    # add to final list  
    split.list[[i]] <- list(trainData = trainData, trainMeta = trainMeta, testData = testData, testMeta = testMeta)
    
  }

  return(split.list)
  
}

# Format metabolite lists so columns don't contain special characters or start with numbers

format.data <- function(results.list){
  
  for(i in 1:length(results.list)){
    
    # replace special characters ("-" and "()") in column names with "."
    colnames(results.list[[i]]) <- gsub("[[:punct:]]", "\\.", colnames(results.list[[i]]))
    
    # replace spaces with "."
    colnames(results.list[[i]]) <- gsub(" ", "\\.", colnames(results.list[[i]]))
    
    # add "x" before column names that start with digits or numbers
    for(j in 1:ncol(results.list[[i]])){
      
      if(grepl("^[[:digit:]]+", colnames(results.list[[i]])[j])){colnames(results.list[[i]])[j] <- paste0("x", colnames(results.list[[i]])[j])}
      if(grepl("^[[:punct:]]+", colnames(results.list[[i]])[j])){colnames(results.list[[i]])[j] <- paste0("x", colnames(results.list[[i]])[j])}
      
    }
    
  }
  
  return(results.list)
  
}

# Run random forest using training/testing split 

tt.rf <- function(results.list, meta.data, variable = "group", split = 0.75, return.model = F, log = T){
  # results.list = list of metabolites
  # meta.data = describe samples and groups
  # variable = desired variable to predict
  # split = fraction to use as training data
  # return.model = logical indicating whether you want the whole model in the results **needed for feature importance
    # otherwise, returns data frame with accuracy
  # log = logical indicating whether to log transform the data
  
  # Split and format data
  
    ## format
    
    results.list <- format.data(results.list)

    ## make a list of matrices
    
    results.list <- lapply(results.list, as.matrix)    
    
    ## log transform, if indicated
    
    if(log == T){
      
      results.list <- lapply(results.list, log)
      
    }else{
      
      results.list <- results.list
      
    }
    
    ## split
    split.list <- split.data(results.list, meta.data, split)
  
  # Run models
    
    ## set up results
      ### confusion matrix results
        rf.list <- list()
        
    ## Loop through modeling
    for(i in 1:length(split.list)){
      #print(i)
      
      # get the training and testing data
        trainData <- split.list[[i]]$trainData
        trainMeta <- split.list[[i]]$trainMeta
        testData <- split.list[[i]]$testData
        testMeta <- split.list[[i]]$testMeta
        
        #print(trainData$group)
        
      # run the model
        rf.model <- randomForest(x = trainData, y = trainMeta[,c(variable)], 
                                 ntree=1000, replace = F)
    
      # Predict
      predicted <- predict(rf.model, testData)
      
      # Compute the confusion matrix
      #testData$group <- as.factor(testData$group)
      confusion.matrix <- table(observed = testMeta[,c(variable)], predicted = predicted)
      
      caret.matrix <- caret::confusionMatrix(predicted, testMeta[,c(variable)])
      
        ## Calculate true and false positive/negative rates from confusion matrix
        
          ### if any groups weren't observed in the testing set (all 0s), remove from the confusion matrix
          zeros <- apply(confusion.matrix, 2, function(x) length(which(x == 0)))
          
          if(any(zeros == 7)){
            confusion.matrix <- confusion.matrix[-which(zeros == 7), -which(zeros == 7)]
          }
          
        
        true_positives <- diag(confusion.matrix)  
        true_positives_sum <- sum(true_positives)
        
        false_positives <- colSums(confusion.matrix) - true_positives
        false_positives_sum <- sum(false_positives)
        
        false_negatives <- rowSums(confusion.matrix) - true_positives
        false_negatives_sum <- sum(false_negatives)
        
        true_negatives <- sum(confusion.matrix) - true_positives - false_positives - false_negatives
        true_negatives_sum <- sum(true_negatives)
        
        specificity <- true_negatives_sum / (false_positives_sum + true_negatives_sum)
        
        sensitivity <- true_positives_sum / (true_positives_sum + false_negatives_sum)
      
        ## Add to results
          if(return.model == T){
            
            rf.list[[i]] <- list(model = rf.model, 
                                 confusion.matrix = confusion.matrix, 
                                 accuracy = caret.matrix[["overall"]][["Accuracy"]],
                                 sensitivity = sensitivity, 
                                 specificity = specificity)
          }else{
            
            rf.list[[i]] <- list(accuracy = caret.matrix[["overall"]][["Accuracy"]],
                                 sensitivity = sensitivity, 
                                 specificity = specificity)
      }

  }
  
  names(rf.list) <- names(results.list)
    
  return(rf.list)
} 

# Run models using LOO validation

loo.rf <- function(results.list, meta.data, include.reps = F, variable = "group", log = T){
  # results.list = metabolite abundance data
  # meta.data = sample descriptions
  # include.reps = logical indicating whether to include replicate samples in the model
  # variable = meta data variable to use for predictions
  # log = logical indicating whether or not to log transform data
  
  # prepare data
  
    ## format
    
    results.list <- format.data(results.list)
    
    ## make a list of matrices
    
    results.list <- lapply(results.list, as.matrix)
    
    ## log transform, if indicated
    
    if(log == T){
      
      results.list <- lapply(results.list, log)
      
    }else{
      
      results.list <- results.list
      
    }
  
    ## adjust for including replicates or not
    
    new.results.list <- list()
    
    if(include.reps == F){
      
      for(i in 1:length(results.list)){
        
        new.results.list[[i]] <- results.list[[i]][which(grepl("rep1", row.names(results.list[[i]]))),]
        new.results.list[[i]] <- new.results.list[[i]][-which(grepl("abc", row.names(new.results.list[[i]]))),]
        
        new.meta.data <- meta.data[which(meta.data$rep == 1),]
        new.meta.data <- new.meta.data[-which(new.meta.data$group == "abc"),]
        new.meta.data$group <- factor(new.meta.data$group, levels = c("a", "b", "c", "ab", "bc", "ac"))
      }
      
    }else{
      
      new.results.list <- results.list
      new.meta.data <- meta.data
      
    }
  
    names(new.results.list) <- names(results.list)
    
    
  # run models for each lab, leaving one sample out as testing
 
    ## prepare results for accuracy table 
      accuracy <- data.frame(matrix(ncol = length(new.results.list) + 1, nrow = nrow(new.results.list[[1]])))
      colnames(accuracy) <- c("truth", names(new.results.list))
      row.names(accuracy) <- paste0("test_", row.names(new.results.list[[1]]))
      accuracy$truth <- new.meta.data[,c(variable)] 
    
    ## prepare results for feature importance
      feat.imp <- list()
      length(feat.imp) <- length(new.results.list)
      names(feat.imp) <- names(new.results.list)
      
  # add meta data to results.list
  #fin.results.list <- lapply(new.results.list, function(x) cbind("variable" = new.meta.data[,c(variable)], x))  
      
  for(i in 1:length(new.results.list)){
    
    # prep feature importance results
    feat.imp[[i]] <- data.frame(matrix(nrow = ncol(new.results.list[[i]]), ncol = nrow(new.results.list[[i]])))
    colnames(feat.imp[[i]]) <- paste0("test_", row.names(new.results.list[[i]]))
    row.names(feat.imp[[i]]) <- colnames(new.results.list[[i]])
    
    for(j in 1:nrow(new.results.list[[i]])){
      
      # split data to leave one sample out
      trainData <- new.results.list[[i]][-j,]
      trainMeta <- new.meta.data[-j,]
      
      testData <- new.results.list[[i]][j, , drop = F]
      testMeta <- new.meta.data[j,]
      
      # run the model
      rf.model <- randomForest(x = trainData, y = trainMeta[,c(variable)], 
                               ntree=1000, replace = F)
      
      # Predict
      predicted <- predict(rf.model, testData)
      
      # add to results
      accuracy[j, i + 1] <- as.character(predicted)
        
      # extract variable importance
      feat.imp[[i]][,j] <- rf.model$importance
        
    }
    
  }
  
  # report accuracy
  
  truth <- accuracy$truth
  labs <- subset(accuracy, select = -c(truth))
  
  lab.accuracy <- apply(labs, 2, function(x) mean(x == truth))
  
  for(i in 1:length(lab.accuracy)){
    
    print(paste("Average Accuracy for", colnames(labs)[i], " = ", round(lab.accuracy[i], 2)))
    
  }
  
  return(list(accuracy = accuracy, feature.importance = feat.imp))
  
}
    
# Find the most important features across labs and models
## calculate average importance for each feature across all models for a given lab

feat.imp <- function(features.list, report.type = "mean", rank.cutoff = NULL){
  # features.list = list of feature importances for each lab across multiple model iterations
  # report.type = type of report to give; 
    ## "mean" = mean gini score for feature across models
    ## "rank" = how often a feature was in the top X% of gini scores across models
  # rank.cutoff = percentage cutoff to be determined important if using rank report type
  
  # create empty results table for all metabolites 
  ## (to be filled in with the number of times they were determined "important" across models for each lab)
  i = 1
  total.metabs <- row.names(features.list[[i]])
  
 for(i in 2:length(features.list)) {
    
    total.metabs <- c(total.metabs, row.names(features.list[[i]]))
    
 }
  total.metabs <- unique(total.metabs)
  
  results <- data.frame(matrix(nrow = length(total.metabs), ncol = length(features.list)))
  colnames(results) <- names(features.list)
  row.names(results) <- total.metabs
  
    # fill results with mean gini scores
    
    if(report.type == "mean"){
      for(i in 1:length(features.list)){
        for(j in 1:nrow(results)){
          for(k in 1:nrow(features.list[[i]])){
            
            if(row.names(features.list[[i]])[k] == row.names(results)[j]){
              
              results[j, i] <- round(rowMeans(features.list[[i]])[k], 2)
              
            }  
          }
        }
      }
    }

  # fill results with frequency of gini scores in top X%
  
  if(report.type == "rank"){
    for(i in 1:length(features.list)){
      
      gini.cutoffs <- apply(features.list[[i]], 2, function(x) quantile(x, rank.cutoff))
      
      for(j in 1:nrow(results)){
        for(k in 1:nrow(features.list[[i]])){
            
            if(row.names(features.list[[i]])[k] == row.names(results)[j]){
              
              # find how many elements are greater than the gini.cutoffs for each feature
              
              compare <- features.list[[i]][k,] > gini.cutoffs
              
              results[j, i] <- sum(compare)
                      
                  }  
                }
              }
            }
          }
          
  return(results)
  
}





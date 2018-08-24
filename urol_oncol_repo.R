####
###
## 
# Custom functions -- by Damiano Fantini -- 2018/Aug/20
#
# ** these functions were developed by Damiano Fantini, Ph.D.
# ** as part of the bioinformatic analyses included in the original
# ** paper: [Fantini et al, Urologic Oncology, 2018]
#
# ** This code can be used, or modified and used only if the 
# ** original author (ie, Damiano Fantini) is informed 
# ** (email: damiano.fantini at gmail.com), and acknowledged
# ** in the final output of the analysis/research.

# ** Typically, this means either:
# ** 1) acknowledging the author of this code 
#       by including his name (Dr. Damiano Fantini) 
#       and a link to the GitHub or data-pulse.com repos
#       (www.github.com/dami82/)
#       in the acknowledgments or "Data availability" sections of the paper
#       AND citing the work [Fantini et al, Urologic Oncology, 2018] 
# ** OR:
# ** 2) including Dr. Damiano Fantini as author of the paper
#

####
###
##
# Warning: the following libraries may be expected to be available and/or loaded on the system

# library(TCGAretriever)
# library(mutSignatures)
# ibrary(gplots)
# library(ggplot2)
# library(caret)
# library(DMwR)
# library(rpart)
# library(pROC)
# library(qvalue)
# library(survival)
# library(survminer)

####
###
##
# Custom functions // code starts here

matrix_to_df <- function(x) {
  cls <- colnames(x)
  rws <- rownames(x)
  out <- do.call(rbind, lapply(1:ncol(x), function(ci){
    do.call(rbind, lapply(1:nrow(x), function(ri){
      data.frame(A=rws[ri], 
                 B=cls[ci],
                 C=x[ri,ci],
                 stringsAsFactors = FALSE)
    }))
  }))
  return(out)
}

prep_tcga_muts <- function(x, colFeat_name, rowFeat_name, valuesFeat_name) {
  all_vals <- unique(x[, valuesFeat_name])
  encoded_vals <- data.frame(id = LETTERS[1:length(all_vals)],
                             vals = all_vals,
                             stringsAsFactors = FALSE)
  
  rowVals <- unique(x[,rowFeat_name])
  colVals <- unique(x[,colFeat_name])
  
  message("", appendLF = TRUE)
  message("Prepare ", appendLF = FALSE)
  
  flag_1 <- unique(as.integer(seq(1, length(colVals), length.out = 20)))
  
  myMat <- sapply(1:length(colVals), function(j1) {
    ci <- colVals[j1]
    tmp <- x[x[, colFeat_name] == ci, ]
    
    if(j1 %in% flag_1)
      message(".", appendLF = FALSE)
    
    sapply(rowVals, function(ri) {
      tmp2 <- tmp[tmp[, rowFeat_name] == ri, ]
      if (nrow(tmp2) < 1) {
        outVal <- NA
      } else {
        tmpOut <- unique(tmp2[, valuesFeat_name])
        outVal <- sapply(tmpOut, function(zz) {
          encoded_vals[encoded_vals[, "vals"] == zz, "id"]
        })
        outVal <- paste(sort(outVal), collapse = "")
      }
      outVal
    })
  })
  
  myDF <- as.data.frame(myMat, stringsAsFactors = FALSE)
  colnames(myDF) <- colVals
  
  message("", appendLF = TRUE)
  message("Reformat ", appendLF = FALSE)
  
  flag_2 <- unique(as.integer(seq(1, nrow(myDF), length.out = 20)))
  myDF2 <- do.call(rbind, lapply(1:nrow(myDF), function(j2) {
    rn <- rownames(myDF)[j2]
    tmp <- myDF[rn,]
    
    if (j2 %in% flag_2)
      message(".", appendLF = FALSE)
    
    do.call(rbind, lapply(colnames(myDF), function(cn) {
      data.frame(id1 = rn, id2 = cn,
                 value = as.character(tmp[cn]),
                 stringsAsFactors = FALSE)
    }))
  }))
  colnames(myDF2)[1:2] <- c(rowFeat_name, colFeat_name)
  
  message(" Done!", appendLF = TRUE)
  
  return(list(encoded_vals = encoded_vals,
              mutGene_tab = myDF,
              mutGene_df = myDF2))
}

initialize_par <- function() {
  tryCatch(dev.off(), error = function(e) NULL) 
  plot(1) 
  oldPAR <- par(no.readonly = T) 
  tryCatch(dev.off(), error = function(e) NULL)
  return(oldPAR)
}

custScale <- function(x, outMax = 1.5, outMin = 0.3) {
  mMin <- min(x, na.rm = T)
  mMax <- max(x, na.rm = T)
  
  #Transpose to linear range 0-1
  x1 <- (x-mMin) / (mMax - mMin)
  
  #Transpose to linear range outMin - outMax
  x2 <- (x1 * (outMax - outMin)) + outMin
  
  return(x2)
}

rescale_surv <- function(surv_data, 
                         time_col = "OS_MONTHS", 
                         status_col = "OS_STATUS", 
                         event_lev = "DECEASED", 
                         tm_factor = 1, 
                         threshold = 60) 
{
  # Initial checks
  tmp <- surv_data
  if(!time_col %in% names(tmp) | !status_col %in% names(tmp)) {
    message("column(s) not found")
    return(NULL)
  } 
  if (is.factor(tmp[,status_col])){
    if (!event_lev[1] %in% levels(tmp[,status_col])){
      message("wrong level")
      return(NULL)
    } else if (!event_lev[1] %in% unique(tmp[,status_col])) {
      message("wrong level")
      return(NULL)
    }
  }
  
  # Thresholding
  tmp[,status_col] <- as.character(tmp[,status_col])
  tmp <- tmp[!is.na(tmp[,status_col]),]
  two_levs <- unique(tmp[,status_col])
  if (length(two_levs) != 2) {
    message("Wrong num of levs")
    return(NULL)
  }
  non_event_lev <- unique(tmp[,status_col])[unique(tmp[,status_col]) != event_lev]
  
  get_idx <- tmp[,time_col] > threshold
  get_idx[is.na(get_idx)] <- FALSE
  
  tmp[get_idx, status_col] <- non_event_lev
  tmp[get_idx, time_col] <- threshold
  
  # Scaling
  tmp[,time_col] <- tmp[,time_col] * tm_factor
  
  # Returning
  return(tmp)
  
}

compare_numVar_by_status <- function(df, col_Var, col_Status, ...) {
  my_f <- df[,col_Status]
  my_n <- as.numeric(df[,col_Var])
  my_list <- split(my_n, f = my_f)
  
  # Compute p-val
  all_levs <- names(my_list)
  my_test <- list()
  for (i in 1:(length(all_levs)-1)) {
    for (j in (i+1):length(all_levs)) {
      
      t_01 <- all_levs[i]
      t_02 <- all_levs[j]
      t_tmp <- paste(t_01, "__vs__", t_02, sep = "")
      
      tt_tmp <- t.test(my_list[[t_01]], my_list[[t_02]])
      my_test[[t_tmp]] <- tt_tmp$p.value
    }
  }
  
  # Plot
  boxplot(my_list, ...)
  
  return(my_test)
}

chisquare_Var_by_status <- function(df, col_Var, col_Status, cellnote = TRUE, ...) {
  my_f <- df[,col_Status]
  my_n <- df[,col_Var]
  my_tab <- table(my_n, my_f)
  
  my_mat <- as.matrix(my_tab)
  
  # Compute p-val
  all_levs <- dimnames(my_tab)$my_f
  my_test <- list()
  for (i in 1:(length(all_levs)-1)) {
    for (j in (i+1):length(all_levs)) {
      
      t_01 <- all_levs[i]
      t_02 <- all_levs[j]
      t_tmp <- paste(t_01, "__vs__", t_02, sep = "")
      
      tt_tmp <- suppressWarnings(chisq.test(my_mat[,c(i,j)]))
      
      my_test[[t_tmp]] <- tt_tmp$p.value
    }
  }
  
  # Plot
  if (cellnote) {
    heatmap.2(my_mat, cellnote = my_mat, ...)
  } else {
    heatmap.2(my_mat, cellnote = matrix("", ncol = ncol(my_mat), nrow = nrow(my_mat)), ...)    
  }
  return(my_test)
}


med.zsco <- function(num) {
  med <- median(num)
  mad <- sum(abs(num-med)) / length(num)
  if (mad == 0) {
    rep(0, length(num))
  } else {
    as.numeric((num-med) / (1.486*mad) )
  }
}


calc.dist <- function(data, scale = TRUE) {
  
  # store names for later
  my.names <- colnames(data)
  
  # if scaling is set to 1, let's convert all rows
  data <- do.call(rbind, lapply(1:nrow(data), (function(i){
    med.zsco(data[i,])
  })))
  
  out <- dist(t(data))
  
  #return
  return(out)
}

get.mostVar <- function(data, n = 1000) {
  
  # compute coefficient of variation
  coeff.var <- sapply(1:nrow(data), (function(i){
    nn <- as.numeric(data[i,])
    sd(nn, na.rm = T) / mean(nn, na.rm = T)
  }))
  
  # reorder based on coeff.var
  out <- data[order(coeff.var, decreasing = T),]
  
  #
  if (!is.null(n) &  n < nrow(data))
    data <- data[(1:n),]
  
  return(out)
}


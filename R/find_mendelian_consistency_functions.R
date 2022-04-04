# Whoeps, 20th Jul 2021

#' @export
extract_used_samples <- function(csvlink){
  # This is for processing of one of David's typical outputs. #
  # Load david's genotypes

  dgt = read.table(david_list, header=T, sep=',', stringsAsFactors = F)

  # define entries to keep
  cols_to_keep = c('genoT')

  # filter
  dgtf = dgt %>%   select(matches(paste(cols_to_keep, collapse="|")))

  # rename columns
  colnames(dgtf) = str_replace(colnames(dgtf),'genoT_','')

  return(colnames(dgtf))

}

#' @export
cut_cm_to_samples <- function(cm_f, used_samples_uniq_f, n_overjump, include_children){

  # Kick children out if not wanted (hehe)
  if (include_children){
    children_ids = c('HG00514','HG00733','NA19240')
  } else {
    children_ids = NULL
  }
  # Change callmatrix samplenames

  # Go in reverse order throu sample columns
  for (i in length(colnames(cm_f)):n_overjump){

    # Rename HG002 if you encounter it
    if (colnames(cm_f)[i]=='HG002'){
      colnames(cm_f)[i] = 'NA24385'
      next
    } else if (colnames(cm_f)[i] %in% children_ids){
      next
    }

    # Use the samplenames from David's table to adjust ours.
    newlabel = used_samples_uniq_f[grepl(substr(colnames(cm_f)[i],3,8),used_samples_uniq_f)]

    # If the sample wasn't found in David's table, we drop it.
    if (length(newlabel)>0){
      colnames(cm_f)[i] = newlabel
    } else {
      print(paste0('Dropped sample ',colnames(cm)[i]))
      cm_f[,colnames(cm_f)[i]] = NULL
    }
    # Reset newlabel
    newlabel = NULL
  }

  return(cm_f)
}


#' Make a vector to check inheritance plausibilities
#' @return child_expect, a vector describing excpected child genotypes given the parents.
#' @author Wolfram Hoeps
#' @export
make_child_expect_vector <- function(){
  # ok who cares let's do it the hard way.
  #parents = paste(p1, p2)
  child_expect <- vector(mode="list")
  child_expect[['0|0 0|0']] = c('0|0')
  child_expect[['0|0 0|1']] = c('0|0','0|1')
  child_expect[['0|0 1|1']] = c('0|1')
  child_expect[['0|1 0|0']] = c('0|0','0|1')
  child_expect[['0|1 0|1']] = c('0|0','0|1','1|1')
  child_expect[['0|1 1|1']] = c('0|1','1|1')
  child_expect[['1|1 0|0']] = c('0|1')
  child_expect[['1|1 0|1']] = c('0|1','1|1')
  child_expect[['1|1 1|1']] = c('1|1')
  return (child_expect)
}

#' Test a gt for validity
#' @export
test_mendel <- function(ce, gt_parent1, gt_parent2, gt_child){

  # Confidence in gt. If lowconf is present, this is 1.
  conf_parent1 = grep('lowconf', gt_parent1)
  conf_parent2 = grep('lowconf', gt_parent2)
  conf_child =   grep('lowconf', gt_child)

  gt_parent1 = substr(gt_parent1,1,3)
  gt_parent2 = substr(gt_parent2,1,3)
  gt_child = substr(gt_child, 1,3)
  gt_parent1 = gsub('1\\|0','0\\|1', gt_parent1)
  gt_parent2 = gsub('1\\|0','0\\|1', gt_parent2)
  gt_child   = gsub('1\\|0','0\\|1', gt_child)

  valid_gts = c('0|0','0|1','1|1')
  c1 = gt_parent1 %in% valid_gts
  c2 = gt_parent2 %in% valid_gts
  c3 = gt_child %in% valid_gts

  if (c1){
    if (c2){
      if (c3){
        valid = gt_child %in% ce[[paste(gt_parent1, gt_parent2)]]

        lowconf_exists = sum(conf_parent1, conf_parent2, conf_child)
        if (lowconf_exists>0){
          valid = paste0(valid, "_lowconf")
        }
        return(valid)
      }
    }
  }
  return('CPX')
}

#' @export
add_mendelfails <- function(cm_f){
  ### This should definitely be rewritten to a less hardcoded version ###

  ce = make_child_expect_vector()
  callmatrix = cm_f
  callmatrix$mendel1 = 'UNK'
  callmatrix$mendel2 = 'UNK'
  callmatrix$mendel3 = 'UNK'

  for (row in 1:nrow(callmatrix)){
    #callmatrix[row,]$mendel =
    callmatrix[row,]$mendel1 = (test_mendel(ce, callmatrix[row,]$HG00512, callmatrix[row,]$HG00513,callmatrix[row,]$HG00514 ))
    callmatrix[row,]$mendel2 = (test_mendel(ce, callmatrix[row,]$HG00731, callmatrix[row,]$HG00732,callmatrix[row,]$HG00733 ))
    callmatrix[row,]$mendel3 = (test_mendel(ce, callmatrix[row,]$NA19238, callmatrix[row,]$NA19239,callmatrix[row,]$NA19240 ))

    #callmatrix[row,]$mendel4 = as.logical(test_mendel(ce, callmatrix[row,]$GM19650, callmatrix[row,]$HG00864,callmatrix[row,]$HG03371 ))

  }

  callmatrix$mendelfails_highconf = rowSums(callmatrix[,c('mendel1','mendel2','mendel3')] == 'FALSE')
  callmatrix$mendelfails_lowconf = rowSums(callmatrix[,c('mendel1','mendel2','mendel3')] == 'FALSE_lowconf')

  return(callmatrix)
}


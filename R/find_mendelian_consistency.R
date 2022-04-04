#' Test callset for mendelian consistency
#'
#' A wrapper function to test an inversion callset for mendelian consistency.
#' Warning: The trios (t1 = [HG00512, HG00513, HG00514], t2 = [HG00731, HG00732, HG00733],
#' t3 = [NA19238, NA19239, NA19240]) are hard-coded in this release, which was specifically
#' designed for the inversion callset published in Porubsky et al., 2021.
#'
#' @param infile_link link to an inversion callset file with n=14+n_samples columns (format: see data/testsample.csv)
#' @param outfile_link output file to be written: the same table, but with columns specifying mendelian consistency.
#' @author Wolfram Hoeps
#' @export
run_main <- function(infile_link, outfile_link=F){

  # Load libraries
  library(stringr)
  library(dplyr)
  library(pheatmap)
  library(matrixStats)
  library(reshape2)

  gt_mat = read.table(infile_link, sep='\t', stringsAsFactors=F, header=1)#, row.names=NULL)

  # Obscure renaming related to row.names. Seems to be a new problem. Did something
  # change with a newer R version in that regard?
  #colnames(gt_mat) <- colnames(gt_mat)[2:ncol(gt_mat)]
  #gt_mat <- gt_mat[ , - ncol(gt_mat)]
  print(head(gt_mat))
  gt_mat_mendel = add_mendelfails(gt_mat)

  # Mendelfail: At least one trio fails mendelian consistency
  gt_mat_mendel$Mendelfail_highconf = gt_mat_mendel$mendelfails_highconf > 0
  gt_mat_mendel$Mendelfail_lowconf = gt_mat_mendel$mendelfails_lowconf > 0

  gt_mat_mendel$CPX_calls = rowSums(gt_mat_mendel[,c('mendel1', 'mendel2', 'mendel3')] == 'CPX')

  if (outfile_link != F){
    write.table(gt_mat_mendel, outfile_link, sep='\t', col.names = T, row.names = F, quote = F)
  } else {
    return(gt_mat_mendel)
  }
}



# runs only when script is run by itself
if (sys.nframe() == 0){

    library(optparse)
    source('find_mendelian_consistency_functions.R')
    # INPUT
    option_list = list(
    make_option(c("-i", "--infile"), type="character", default=NULL,
                help="Genotype file the be mendeltested", metavar="character"),
    make_option(c("-o", "--outfile"), type="character", default="./outputcorr/",
                help="Outputfile for that genotype file", metavar="character")
    )
    opt <- parse_args(OptionParser(option_list=option_list))

    infile = opt$infile
    outfile = opt$outfile

    run_main(infile, outfile)

}
# } else {
#    # Hardcode path
#     infile = "/Users/hoeps/PhD/projects/huminvs/analyses_paper/mendel_finalcallsets/data/variants_freeze4inv_sv_inv_hg38_processed_arbigent_filtered_manualDotplot_filtered_PAVgenAdded_withInvCategs.tsv"
#     outfile = "/Users/hoeps/PhD/projects/huminvs/analyses_paper/mendel_finalcallsets/processed/hg38.tsv"
# }




# Mendelian consistency tests.

## Overview of this workflow

This is a repo collecting scripts and a wrapper for testing an inversion callset, such as reported in Porubsky et al., 2021 (https://doi.org/10.1101/2021.12.20.472354), for mendelian consistency.

Specifically, all inversion calls of all sample 'trios' (three trios in our case) are examined, noting whether the haplotype of the offspring can be explained by mendelian inheritance by the two parents.

## Running the script


### Input Data

An example input file is provided in data/

### Instructions

Run the wrapper script using:
 `cd R`
 `Rscript find_mendelian_consistency.R -i ../data/testsample.tsv -o ../testsample_out.tsv`

### Output

A tsv file with additional columns:
 column   description
 mendel1  mendelian consistency of trio 1 [T/F]
 mendel2  mendelian consistency of trio 2 [T/F]
 mendel3  mendelian consistency of trio 3 [T/F]
 mendelfails_highconf number of high-quality non-mendelian calls
 mendelfails_lowconf  number of low-quality non-mendelian calls
 Mendelfail_highconf  mendelfails_highconf >= 1? [T/F]
 Mendelfail_lowconf mendelfails_lowconf >= 1? [T/F]
 CPX_calls  number of genotypes with a complex call.


## References

This code is associated to our publication on human inversion recurrence:
Porubsky et al., 2021 (https://doi.org/10.1101/2021.12.20.472354)

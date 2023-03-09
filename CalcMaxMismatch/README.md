The Fortran program in `CalcMismatch.f90` is a simplification of the function `CalcMaxMismatch()` 
in the `sequoia` R package. It uses the mean (or median) MAF rather than the MAF at every SNP,
which in most cases does not affect the calculated quantiles of the number of mismatches. 

The quantile used by `sequoia()` as threshold was 0.999^(1/nrow(GenoM)) up to version 2.4, and is 
0.9999^(1/nrow(GenoM))) from version 2.5 onwards. 

The file `MismatchSpecs.txt` shows the format of the input specification; 
comments are preceded by an exclamation mark. 

The file `Mismatch_distr.txt` shows an example of the output.  
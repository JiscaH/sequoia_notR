## General
This Fortran program is equivalent to the R package, except that the input & output does not go via R but via text files and command-line arguments. 
Very large datasets can be incredibly slow or even impossible to load into R; this non-R version works for up to 50K - 100K genotyped individuals. 


## Compilation
A Fortran compiler is needed, such as e.g. gfortran. In Windows you also need a linux emulator, e.g. cygwin (although other approaches may be possible).
```
gfortran -O3 Sequoia_SA.f90 -o sequoia
```

(NOTE: previous versions relied on `-std=95 -fall-intrinsics`, these flags are no longer
necessary and will in fact cause compilation to fail)


## File formats
Template files can be generated using the R package `sequoia`:
```
library(sequoia)
data(SimGeno_example, LH_HSg5)
SeqOUT <- sequoia(SimGeno_example, LH_HSg5, Err = 0.005, Module = "pre")
writeSeq(SeqList = SeqOUT, GenoM = SimGeno_example, 
         folder = "SequoiaFileTemplates")
```
Genotype data can be formatted from standard PLINK format as follows:
```
plink --bfile FileNameIN --recode A --out FileNameOUT
mv FileNameOUT.raw tmp.raw
cat tmp.raw | tr -s ' ' | cut -d ' ' -f2,7- > FileNameOUT.raw
sed -i '1d' FileNameOUT.raw 
sed -i 's/NA/-9/g' FileNameOUT.raw
rm tmp.raw
```
which

- recodes to 1 column per SNP, coded as 0/1/2 copies of minor allele, missing = NA
- drops columns 2--6 (family ID, sex, phenotype, etc)
- drops header row
- replaces missing value code NA by -9


## Running
Parameter values are stored in `SequoiaSpecs.txt`. Many of these can be overruled on the command line, as described in the PDF manual. 
```
./sequoia --geno geno.txt --par --ped --verbose
```

All command line options are shown with
```
./sequoia --help
```


## Version
Its version will typically be functionally identical to the latest R package beta version. To obtain an older version, please send an email. 

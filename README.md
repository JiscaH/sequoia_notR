## General
This Fortran program is equivalent to the R package, except that the input & output does not go via R but via text files and command-line arguments. 


## Compilation
A Fortran compiler is needed, such as e.g. gfortran. In Windows you also need a linux emulator, e.g. cygwin (although other approaches may be possible).
```
gfortran -std=f95 -fall-intrinsics -O3 Sequoia_SA.f90 -o sequoia
```


## File formats
Template files can be generated using the R package `sequoia`:
```
library(sequoia)
data(SimGeno_example, LH_HSg5)
SeqOUT <- sequoia(SimGeno_example, LH_HSg5, Err = 0.005, Module = "pre")
writeSeq(SeqList = SeqOUT, GenoM = SimGeno_example, 
         folder = "SequoiaTemplatesFolder")
```

## Running
Parameter values are stored in `SequoiaSpecs.txt`. Many of these can be overruled on the command line, as described in the PDF manual. 
```
./sequoia --geno GriffinGenotypes --par --ped --verbose
```


## Version
Its version will typically be functionally identical to the latest R package beta version. To obtain an older version, please send an email. 
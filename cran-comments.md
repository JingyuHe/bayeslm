## Resubmission of an archived package

The `bayeslm` package was archived on 10/2/2020. The last CRAN check results are available at 
https://cran-archive.r-project.org/web/checks/2020/2020-10-02_check_results_bayeslm.html

To summarize, there were four check notes: 

```
Check: R code for possible problems 
Result: NOTE 
    Possibly missing '()' after 'return' in the following function:
     'predict.bayeslm.fit' 
```

```
Check: line endings in C/C++/Fortran sources/headers 
Result: NOTE 
    Found the following sources/headers not terminated with a newline:
     inst/include/bayeslm.h
    Some compilers warn on such files. 
```

```
Check: R code for possible problems 
Result: NOTE 
    Possibly missing '()' after 'return' in the following function:
     'predict.bayeslm.fit' 
```

```
Check: for GNU extensions in Makefiles 
Result: NOTE 
    GNU make is a SystemRequirements.
```

```
Check: installed package size 
Result: NOTE 
     installed size is 11.2Mb
     sub-directories of 1Mb or more:
     libs 11.0Mb 
```

The first two of these have been corrected. The second two are addressed below.

## R CMD check results

0 errors | 0 warnings | 2 notes

```
* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements
```

This package uses `RcppParallel` and, as a result, we cannot avoid producing 
this note when running `R CMD check`.

```
* checking installed package size ... NOTE
  installed size is 20.2Mb
  sub-directories of 1Mb or more:
    libs  19.9Mb
```

This note only occurs on linux distributions and, to the best of our understanding, 
is an artifact of the way that R packages are built with debug information on linux.

# rdadapt

The goal of rdadapt is to implment the procedures given in Kwon and Kwon (2019a)
and Kwon and Kwon (2019b).

## Installation

You can install the released version of rdadapt from github with:

``` r
devtools::install_github("https://github.com/soonwookwon/rdadapt")
```

## Some conventions

Many of the functions given in the package take vectors (many cases of length
2, but not necessarily) `C` and `gamma` as arguments. The convention is that the
elements that come first *always* correpond to smaller parameter spaces than
those come later.

Test line

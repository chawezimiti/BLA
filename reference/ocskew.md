# Octile Skewness

A function to calculate an estimate of the octile skewness from a set of
data.

## Usage

``` r
ocskew(x)
```

## Arguments

- x:

  A vector of numeric values.

## Value

The octile skewness.

## Author

Richard Murray Lark \<murray.lark@nottingham.ac.uk\>

## Examples

``` r
x<-evapotranspiration$`ET(mm)`
ocskew(x)
#>     87.5% 
#> 0.1520965 
```

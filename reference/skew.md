# Skewness

Function to calculate an estimate of the coefficient of skewness from a
set of data.

## Usage

``` r
skew(x)
```

## Arguments

- x:

  A vector of numeric values.

## Value

The coefficient of skewness.

## Author

Richard Murray Lark \<murray.lark@nottingham.ac.uk\>

## Examples

``` r
x<-evapotranspiration$`ET(mm)`
skew(x)
#> [1] 0.5380429
```

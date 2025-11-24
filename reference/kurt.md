# kurtosis

A function to calculate an estimate of the coefficient of kurtosis from
a set of data.

## Usage

``` r
kurt(x)
```

## Arguments

- x:

  A vector of numeric values.

## Value

The reduced coefficient of kurtosis.

## Author

Richard Murray Lark \<murray.lark@nottingham.ac.uk\>

## Examples

``` r
x<-evapotranspiration$`ET(mm)`
kurt(x)
#> [1] 0.2089355
```

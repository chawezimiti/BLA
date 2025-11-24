# Summary statistics

A function to calculate summary statistics of a set of data.

## Usage

``` r
summastat(x, sigf, varname, plot = TRUE)
```

## Arguments

- x:

  A vector of numeric values.

- sigf:

  The number of significant figures to report (optional).

- varname:

  The name of the variable (optional), character so in quotes e.g. "Clay
  content". If not used then the variable is called x on plots.

- plot:

  If `TRUE`, a plot is part of the output. If `FALSE`, plot is not part
  of output (default is `TRUE`).

## Value

A matrix containing the mean value, median value, first and third
quartiles, sample variance, sample standard deviation, coefficient of
skewness, octile skewness, coefficient of kurtosis and the number of
probable outliers in a data set. A histogram with a boxplot over it and
QQ plot of the variable x if `plot=TRUE`.

## Author

Richard Murray Lark \<murray.lark@nottingham.ac.uk\>

## Examples

``` r
x<-evapotranspiration$`ET(mm)`
summastat(x,2)

#>      Mean Median Quartile.1 Quartile.3 Variance SD Skewness Octile skewness
#> [1,]  290    280        230        340     7000 84     0.54            0.15
#>      Kurtosis No. outliers
#> [1,]     0.21            0

```

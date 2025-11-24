# Hessian matrix

This is a set of functions that support the censored bivariate normal
model.

## Usage

``` r
seHessian(a, hessian = FALSE, silent = FALSE)
```

## Arguments

- a:

  The hessian matrix.

- hessian:

  If \`True\`, hessian matrix is used.

- silent:

  Condition of matrix.

## Value

The standard error from hessian matrix

## Examples

``` r
hessian<-matrix(c(37.45965, 83.0686,83.06863,188.92427),2,2)

seHessian(hessian) # calculates the standard error
#> [1] 1.0341811 0.4605052
```

# Determination of the most limiting factor to biological response

This function determines the most limiting factor based on von Liebig
law of the minimum given results of the predicted boundary line values
for the different factors of interest. Boundary lines for various
factors are fitted and the factor that predicts the minimum response for
a particular point is considered as the most limiting factor (Casanova
et al. 1995).

## Usage

``` r
limfactor(...)
```

## Arguments

- ...:

  vectors with predicted values from the boundary line models for each
  factor being evaluated.

## Value

A dataframe consisting of the most limiting factor and the minimum
predicted response

## Author

Chawezi Miti \<chawezi.miti@nottingham.ac.uk\>

## Examples

``` r
N<-rnorm(10,50,5)#assuming these are predicted responses using the fitted BL for N,P,K
K<-rnorm(10,50,4)
P<-rnorm(10,50,6)

limfactor(N,K,P)
#> [[1]]
#>          Rs Lim_factor
#> 1  40.12164          K
#> 2  45.05137          P
#> 3  39.34444          P
#> 4  41.56472          P
#> 5  43.68471          N
#> 6  48.82520          N
#> 7  44.86946          K
#> 8  52.25433          K
#> 9  45.88153          K
#> 10 51.37455          N
#> 
#> [[2]]
#> [1] 61.92102
#> 
```

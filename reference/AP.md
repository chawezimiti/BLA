# Calculate perimeter and area of data scatter

This function determines the perimeter and area around the boundary
points

## Usage

``` r
AP(points)
```

## Arguments

- points:

  A dataframe with two columns.

## Value

The area and perimeter covered bivariate data points, and the selected
boundary points

## Examples

``` r
x<-data.frame(x=evapotranspiration$`ET(mm)`,y=evapotranspiration$`yield(t/ha)`)
AP(x)
#> $Perimeter
#> [1] 33616.54
#> 
#> $Area
#> [1] 256.6961
#> 
```

# Predict boundary response

This function predicts the most efficient response at a level of factor,
`x`, given the parameters of the fitted boundary line.

## Usage

``` r
predictBL(object, x)
```

## Arguments

- object:

  An output in form of a list from the boundary line fitting using the
  [`blqr()`](https://chawezimiti.github.io/BLA/reference/blqr.md),
  [`blbin()`](https://chawezimiti.github.io/BLA/reference/blbin.md),
  [`bolides()`](https://chawezimiti.github.io/BLA/reference/BOLIDES.md)
  or [`cbvn()`](https://chawezimiti.github.io/BLA/reference/cbvn.md)
  functions.

- x:

  A numeric vector of values for the factor with which response is to be
  predicted.

## Value

A vector predicted value of response.

## Author

Chawezi Miti \<chawezi.miti@nottingham.ac.uk\>

## Examples

``` r
x<-evapotranspiration$`ET(mm)`
y<-evapotranspiration$`yield(t/ha)`
z<-bolides(x,y, start = c(0.5,0.02), model= "blm", xmax = 350)


Results<-predictBL(z,x)

head(Results) # prediction for first 6 lines
#> [1]  8.392391  9.833827 10.515234  9.021382  8.025480  9.021382

```

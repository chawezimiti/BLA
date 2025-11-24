# Starting values for optimization functions

This functions helps to determine initial values for a selected boundary
line model when using the functions
[`blbin()`](https://chawezimiti.github.io/BLA/reference/blbin.md),
[`blqr()`](https://chawezimiti.github.io/BLA/reference/blqr.md),
[`bolides()`](https://chawezimiti.github.io/BLA/reference/BOLIDES.md),
[`cbvn()`](https://chawezimiti.github.io/BLA/reference/cbvn.md) and
[`ble_profile()`](https://chawezimiti.github.io/BLA/reference/ble_profile.md)
to determine model parameters.

## Usage

``` r
startValues(model = "explore", p = NULL, digits = 2, ...)
```

## Arguments

- model:

  Selects the functional form of the boundary line. It includes `"blm"`
  for linear model, `"lp"` for linear plateau model, `"mit"` for the
  Mitscherlich model, `"schmidt"` for the Schmidt model, `"logistic"`
  for logistic model, `"logisticND"` for logistic model proposed by
  Nelder (1961), `"inv-logistic"` for the inverse logistic model,
  `"double-logistic"` for the double logistic model, `"qd"` for
  quadratic model, `"trapezium"` for the trapezium model and `"explore"`
  for function use exploration. The default is `"explore"`.

- p:

  The number of selected points used to obtain start values for the
  logistic mitcherlich and schmidt models. It is `NULL` for other
  models.

- digits:

  Number of decimal points for logistic type models (default is 2).

- ...:

  Additional graphical parameters. Applies to the logistic, mitcherlich
  and schmidt models to control the text on the plot.

## Value

A list containing the parameters of the suggested model.

## Details

This function uses the
[`locator()`](https://rdrr.io/r/graphics/locator.html) function. Once
the model is selected, the points that make up the boundary points are
selected using mouse click on the plots.

## References

Fekedulegn, D., Mac Siurtain, M.P., & Colbert, J.J. 1999. Parameter
estimation of nonlinear growth models in forestry. Silva Fennica 33(4):
327–336.

Lark, R. M., & Milne, A. E. (2016). Boundary line analysis of the effect
of water filled pore space on nitrous oxide emission from cores of
arable soil. European Journal of Soil Science, 67 , 148-159.

## Author

Chawezi Miti \<chawezi.miti@nottingham.ac.uk\>

## Examples

``` r
startValues(model="explore")
#> y = f(P₁,P₂|x)
#> 
#> $Param1
#> P₁ 
#>  0 
#> 
#> $Param2
#> P₂ 
#>  0 
#> 

```

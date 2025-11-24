# Boundary line model determination using quantile regression

This function fits a boundary model to the upper bounds of a scatter
plot of `x` and `y` by estimating the conditional quantile (0-1) of the
response variable, `y`, across values of the predictor variables, `x`.
This is achieved using optimization procedure and hence requires some
starting guess parameters of a proposed model.

## Usage

``` r
blqr(x,y,model, equation=NULL,start,tau=0.95,optim.method="Nelder-Mead",
     xmin=min(bound$x),xmax=max(bound$x), plot=TRUE,line_col="red",lwd=1,
     line_smooth=1000,...)
```

## Arguments

- x:

  A numeric vector of values for the independent variable.

- y:

  A numeric vector of values for the response variable.

- model:

  Selects the functional form of the boundary line. It includes `"blm"`
  for linear model, `"lp"` for linear plateau model, `"mit"` for the
  Mitscherlich model, `"schmidt"` for the Schmidt model, `"logistic"`
  for logistic model, `"logisticND"` for logistic model proposed by
  Nelder (1961), `"inv-logistic"` for the inverse logistic model,
  `"double-logistic"` for the double logistic model, `"qd"` for
  quadratic model and the `"trapezium"` for the trapezium model.For
  custom models, set `model = "other"`.

- equation:

  A custom model function writen in the form of an R function. Applies
  only when argument `model="other"`, else it is `NULL`.

- start:

  A numeric vector of initial starting values for optimization in
  fitting the boundary model. Its length and arrangement depend on the
  suggested model:

  - For the `"blm"` model, it is a vector of length 2 arranged as
    intercept and slope.

  - For the `"lp"` model, it is a vector of length 3 arranged as
    intercept, slope and maximum response.

  - For the `"logistic"` and `"inv-logistic"` models, it is a vector of
    length 3 arranged as the scaling parameter, shape parameter and
    maximum response.

  - For the `"logisticND"` model proposed by Nelder (1961), it is a
    vector of length 3 arranged as the scaling parameter, shape
    parameter and maximum response.

  - For the `"double-logistic"` model, it is a vector of length 6
    arranged as the scaling parameter one, shape parameter one, maximum
    response, maximum response, scaling parameter two and shape
    parameter two.

  - For the `"qd"` model, it is a vector of length 3 arranged as
    constant, linear coefficient and quadratic coefficient.

  - For the `"trapezium"` model, it is a vector of length 3 arranged as
    intercept one, slope one, maximum response, intercept two and slope
    two.

  - For the `"mit"` model, it is a vector of length 3 arranged as the
    intercept, shape parameter and the maximum response.

  - For the `"schmidt"` model, it is a vector of length 3 arranged as
    scaling parameter, shape parameter (x-value at maximum response )
    and maximum response.

- tau:

  The quantile value (0- 1) that represents the boundary
  (`default is tau = 0.95`).

- optim.method:

  Describes the method used to optimize the model as in the
  [`optim()`](https://rdrr.io/r/stats/optim.html) function. The methods
  include `"Nelder-Mead"`, `"BFGS"`, `"CG"`, `"L-BFGS-B"`, `"SANN"` and
  `"Brent"`.

- xmin:

  Numeric value that describes the minimum `x` value to which the
  boundary line is to be fitted (default is `min(x)`).

- xmax:

  A numeric value that describes the maximum `x` value to which the
  boundary line is to be fitted (default is `max(x)`). `xmin` and `xmax`
  determine the subset of the data set used to fit boundary model.

- plot:

  If `TRUE`, a plot is part of the output. If `FALSE`, plot is not part
  of output (default is `TRUE`).

- line_col:

  Selects the color of the boundary line.

- lwd:

  Determines the thickness of the boundary line on the plot (default is
  1).

- line_smooth:

  Parameter that describes the smoothness of the boundary line. (default
  is 1000). The higher the value, the smoother the line.

- ...:

  Additional graphical parameters.

## Value

A list of length 5 consisting of the fitted model, equation form,
parameters of the boundary line, the weighted residue sum square.
Additionally, a graphical representation of the boundary line on the
scatter plot is produced.

## Details

Some inbuilt models are available for the `blqr()` function. The suggest
model forms are as follows:

1.  Linear model (`"blm"`) \$\$y=\beta_1 + \beta_2x\$\$ where
    \\\beta_1\\ is the intercept and \\\beta_2\\ is the slope.

2.  Linear plateau model (`"lp"`) \$\$y= {\rm min}(\beta_1+\beta_2x,
    \beta_0)\$\$ where \\\beta_1\\ is the intercept , \\\beta_2\\ is the
    slope and \\\beta_0\\ is the maximum response.

3.  The logistic (`"logistic"`) and inverse logistic (`"inv-logistic"`)
    models \$\$ y= \frac{\beta_0}{1+e^{\beta_2(\beta_1-x)}}\$\$ \$\$ y=
    \beta_0 - \frac{\beta_0}{1+e^{\beta_2(\beta_1-x)}}\$\$ where
    \\\beta_1\\ is a scaling parameter , \\\beta_2\\ is a shape
    parameter and \\\beta_0\\ is the maximum response.

4.  Logistic model (`"logisticND"`) (Nelder (1961)) \$\$ y=
    \frac{\beta_0}{1+(\beta_1 \times e^{-\beta_2x})}\$\$ where
    \\\beta_1\\ is a scaling parameter, \\\beta_2\\ is a shape parameter
    and \\\beta_0\\ is the maximum response.

5.  Double logistic model (`"double-logistic"`) \$\$ y=
    \frac{\beta\_{0,1}}{1+e^{\beta_2(\beta_1-x)}} -
    \frac{\beta\_{0,2}}{1+e^{\beta_4(\beta_3-x)}}\$\$ where \\\beta_1\\
    is a scaling parameter one, \\\beta_2\\ is a shape parameter one,
    \\\beta\_{0,1}\\ and \\\beta\_{0,2}\\ are the maximum response ,
    \\\beta_3\\ is a scaling parameter two and \\\beta_4\\ is a shape
    parameter two.

6.  Quadratic model (`"qd"`) \$\$y=\beta_1 + \beta_2x + \beta_3x^2\$\$
    where \\\beta_1\\ is a constant, \\\beta_2\\ is a linear coefficient
    and \\\beta_3\\ is the quadratic coefficient.

7.  Trapezium model (`"trapezium"`) \$\$y={\rm min}(\beta_1+\beta_2x,
    \beta_0, \beta_3 + \beta_4x)\$\$ where \\\beta_1\\ is the intercept
    one, \\\beta_2\\ is the slope one, \\\beta_0\\ is the maximum
    response, \\\beta_3\\ is the intercept two and \\\beta_3\\ is the
    slope two.

8.  Mitscherlich model (`"mit"`) \$\$y= \beta_0 - \beta_1\*\beta_2^x\$\$
    where \\\beta_1\\ is the intercept, \\\beta_2\\ is a shape parameter
    and \\\beta_0\\ is the maximum response.

9.  Schmidt model (`"schmidt"`) \$\$y= \beta_0 +
    \beta_1(x-\beta_2)^2\$\$ where \\\beta_1\\ is ascaling parameter,
    \\\beta_2\\ is a shape parameter (x-value at maximum response ) and
    \\\beta_0\\ is the maximum response .

10. Custom model ("other") This option allows you to create your own
    model form using the function `function()`. The custom model should
    be assigned to the argument `equation`. Note that the parameters for
    the custom model should be `a` and `b` for a two parameter model;
    `a`, `b` and `c` for a three parameter model; `a`, `b`, `c` and `d`
    for a four parameter model and so on.

The function
[`blbin()`](https://chawezimiti.github.io/BLA/reference/blbin.md)
utilities the optimization procedure of the
[`optim()`](https://rdrr.io/r/stats/optim.html) function to determine
the model parameters. There is a tendency for optimization algorithms to
settle at a local optimum. To remove the risk of settling for local
optimum parameters, it is advised that the function is run using several
starting values and the results with the smallest error (weighted
residue sum square) can be taken as a representation of the global
optimum.

The common errors encountered due to poor start values

1.  function cannot be evaluated at initial parameters

2.  initial value in 'vmmin' is not finite

## References

Cade, B. S., & Noon, B. R. (2003). A gentle introduction to quantile
regression for ecologists. Frontiers in Ecology and the Environment,
1(8), 412-420.

Nelder, J.A. 1961. The fitting of a generalization of the logistic
curve. Biometrics 17: 89–110.

Phillips, B.F. & Campbell, N.A. 1968. A new method of fitting the von
Bertelanffy growth curve using data on the whelk. Dicathais, Growth 32:
317–329.

Schmidt, U., Thöni, H., & Kaupenjohann, M. (2000). Using a boundary line
approach to analyze N2O flux data from agricultural soils. Nutrient
Cycling in Agroecosystems, 57, 119-129.

## Author

Chawezi Miti \<chawezi.miti@nottingham.ac.uk\>

## Examples

``` r
x<-log(SoilP$P)
y<-SoilP$yield
start<-c(4,3,13.6)

blqr(x,y, start=start,model = "lp", tau=0.99,
      xlab=expression("ET mm ha"^-1),
      ylab=expression("Wheat yield/ ton ha"^-1),
      pch=16, col="grey")

#> $Model
#> [1] "lp"
#> 
#> $Equation
#> [1] y = min (β₁ + β₂x, β₀)
#> 
#> $Parameters
#>     Estimate
#> β₁  6.523497
#> β₂  2.605818
#> β₀ 13.345196
#> 
#> $RSS
#> [1] 253.8857
#> 
```

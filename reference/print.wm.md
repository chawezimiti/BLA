# Print class wm

This is an S3 print function that Prints only the first 4 elements of wm
class objects.

## Usage

``` r
# S3 method for class 'wm'
print(x, ...)
```

## Arguments

- x:

  Print object.

- ...:

  Other parameters associated with the
  [`print()`](https://rdrr.io/r/base/print.html) function.

## Value

A object containing only the first four items.

## Examples

``` r
numbers<- 1:10
class(numbers)<-"wm"
numbers
#> [1] 1 2 3 4
```

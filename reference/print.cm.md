# Print class cm

This is an S3 print function that Prints only the first 4 elements of cm
class objects.

## Usage

``` r
# S3 method for class 'cm'
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
class(numbers)<-"cm"
numbers
#> [1] 1 2 3 4
```

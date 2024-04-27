test_that("Output should always be of class 'wm'", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- bolides(x,y, theta = c(0.5,0.02), model= "blm", xmax = 350, plot=F)  # Replace "my_function" with the name of your function

  # Check if the output is of class "wm"
  expect_true(inherits(result, "wm"),
              info = "Output should be of class 'wm'")
})

# Check if the output is a list
test_that("Output should always be a list", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- bolides(x,y, theta = c(0.5,0.02), model= "blm", xmax = 350, plot=F)  # Replace "my_function" with the name of your function

  expect_true(is.list(result),
              info = "Output should be a list")
})


# Check if the output contains Model
test_that("Output should always contains Model", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- bolides(x,y, theta = c(0.5,0.02), model= "blm", xmax = 350, plot=F)  # Replace "my_function" with the name of your function

  expect_true("Model" %in% names(result),
              info = "Output should contain 'Model' element")
})


# Check if the output contains Parameters
test_that("Output should always contains Parameters", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- bolides(x,y, theta = c(0.5,0.02), model= "blm", xmax = 350, plot=F)  # Replace "my_function" with the name of your function

  expect_true("Parameters" %in% names(result),
              info = "Output should contain 'Parameters' element")
})


# Check if the output contains Equation
test_that("Output should always contains Equation", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- bolides(x,y, theta = c(0.5,0.02), model= "blm", xmax = 350, plot=F)  # Replace "my_function" with the name of your function

  expect_true("Equation" %in% names(result),
              info = "Output should contain 'Equation' element")
})

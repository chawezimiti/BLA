
# Check if the output is of class "cm"
test_that("Output should always be of class 'cm'", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- blbin(x,y, bins= c(100,350,25) ,start= c(0.5,0.02), model = "blm", xmax = 250, plot=F)  # Replace "my_function" with the name of your function

  expect_true(inherits(result, "cm"),
              info = "Output should be of class 'cm'")
})


# Check if the output is a list
test_that("Output should always be a list", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- blbin(x,y, bins= c(100,350,25) ,start= c(0.5,0.02), model = "blm", xmax = 250, plot=F)  # Replace "my_function" with the name of your function

  expect_true(is.list(result),
              info = "Output should be a list")
})

# Check if the output contains Model
test_that("Output should always contains Model", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- blbin(x,y, bins= c(100,350,25) ,start= c(0.5,0.02), model = "blm", xmax = 250, plot=F)  # Replace "my_function" with the name of your function

  expect_true("Model" %in% names(result),
              info = "Output should contain 'Model' element")
})


# Check if the output contains Parameters
test_that("Output should always contains Parameters", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- blbin(x,y, bins= c(100,350,25) ,start= c(0.5,0.02), model = "blm", xmax = 250, plot=F)  # Replace "my_function" with the name of your function

  expect_true("Parameters" %in% names(result),
              info = "Output should contain 'Parameters' element")
})

# Check if the output contains Equation
test_that("Output should always contains Equation", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- blbin(x,y, bins= c(100,350,25) ,start= c(0.5,0.02), model = "blm", xmax = 250, plot=F)  # Replace "my_function" with the name of your function

  expect_true("Equation" %in% names(result),
              info = "Output should contain 'Equation' element")
})


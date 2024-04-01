
# Check if the output is of class "cm"
test_that("Output should always be of class 'cm'", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- blqr(x,y,theta = c(0.5,0.02), model = "blm", tau = 0.95, plot=F)

  expect_true(inherits(result, "cm"),
              info = "Output should be of class 'cm'")
})

# Check if the output is a list
test_that("Output should always be a list", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- blqr(x,y,theta = c(0.5,0.02), model = "blm", tau = 0.95, plot=F)

  expect_true(is.list(result),
              info = "Output should be a list")
})

# Check if the output contains Model
test_that("Output should always contain Model", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- blqr(x,y,theta = c(0.5,0.02), model = "blm", tau = 0.95, plot=F)

  expect_true("Model" %in% names(result),
              info = "Output should contain 'Model' element")
})


# Check if the output contains Parameters
test_that("Output should always contain Parameters", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- blqr(x,y,theta = c(0.5,0.02), model = "blm", tau = 0.95, plot=F)

  expect_true("Parameters" %in% names(result),
              info = "Output should contain 'Parameters' element")
})

# Check if the output contains Equation
test_that("Output should always contain Equation", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  result <- blqr(x,y,theta = c(0.5,0.02), model = "blm", tau = 0.95, plot=F)

  expect_true("Equation" %in% names(result),
              info = "Output should contain 'Equation' element")
})

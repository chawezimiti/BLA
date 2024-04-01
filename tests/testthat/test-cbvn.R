# Check if the output is of class "cm"
test_that("Output should always be of class 'cm'", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  vals<-data.frame(x,y)
  theta<-c(0.5,0.02,mean(x),mean(y),sd(x),sd(y),cor(x,y))

  result <- cbvn(vals,theta = theta, sigh = 0.4,model= "blm", plot=F)

  expect_true(inherits(result, "cm"),
              info = "Output should be of class 'cm'")
})

# Check if the output is a list
test_that("Output should always be a list", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  vals<-data.frame(x,y)
  theta<-c(0.5,0.02,mean(x),mean(y),sd(x),sd(y),cor(x,y))

  result <- cbvn(vals,theta = theta, sigh = 0.4,model= "blm", plot=F)

  expect_true(is.list(result),
              info = "Output should be a list")
})

# Check if the output contains Model
test_that("Output should always contains Model", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  vals<-data.frame(x,y)
  theta<-c(0.5,0.02,mean(x),mean(y),sd(x),sd(y),cor(x,y))

  result <- cbvn(vals,theta = theta, sigh = 0.4,model= "blm", plot=F)

  expect_true("Model" %in% names(result),
              info = "Output should contain 'Model' element")
})

# Check if the output contains Parameters
test_that("Output should always be a list", {
  x<-evapotranspiration$`ET(mm)`
  y<-evapotranspiration$`yield(t/ha)`
  vals<-data.frame(x,y)
  theta<-c(0.5,0.02,mean(x),mean(y),sd(x),sd(y),cor(x,y))

  result <- cbvn(vals,theta = theta, sigh = 0.4,model= "blm", plot=F)
  expect_true("Parameters" %in% names(result),
              info = "Output should contain 'Parameters' element")
})


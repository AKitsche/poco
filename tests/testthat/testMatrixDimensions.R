library(testthat)
library(poco)
context("Contrast matrices test")

test_that("check the dimension of the numertor and denominator contrast matrices",{
  
  #example data from binomial data
  data(MetoCRXL2)
  library(MCPAN)
  MCPAN_Est <- binomest(Success ~ RegionTreat,
                        data=MetoCRXL2, 
                        success="1", 
                        method="Wald")
  base <- PowConBinom(p=MCPAN_Est$estp, 
              n=MCPAN_Est$n, 
              n.sub = 12, 
              TreatMat = "Tukey", 
              SubMat = "GrandMean",
              rhs = 0.5, 
              alternative = "less", 
              alpha = 0.05,
              type="anypair")
  
  
  
  expect_equal(ncol(base$CMat), length(MCPAN_Est$estp))
  expect_equal(ncol(base$DMat), length(MCPAN_Est$estp))
  expect_equal(dim(base$CMat), dim(base$DMat))
})

test_that("check the constrains for the numerator and denominator matrices",{
  
  #example data from binomial data
  data(MetoCRXL2)
  library(MCPAN)
  MCPAN_Est <- binomest(Success ~ RegionTreat,
                        data=MetoCRXL2, 
                        success="1", 
                        method="Wald")
  base <- PowConBinom(p=MCPAN_Est$estp, 
                      n=MCPAN_Est$n, 
                      n.sub = 12, 
                      TreatMat = "Tukey", 
                      SubMat = "GrandMean",
                      rhs = 0.5, 
                      alternative = "less", 
                      alpha = 0.05,
                      type="anypair")

    
  expect_equal(.rowSums(base$DMat[base$DMat<0], m=nrow(base$DMat),n=length(base$DMat[base$DMat<0])/nrow(base$DMat)), rep(-1,times=nrow(base$DMat)))
  expect_equal(.rowSums(base$CMat[base$CMat<0], m=nrow(base$CMat),n=length(base$CMat[base$CMat<0])/nrow(base$CMat)), rep(-1,times=nrow(base$DMat)))
  expect_equal(.rowSums(base$DMat[base$DMat>0], m=nrow(base$DMat),n=length(base$DMat[base$DMat>0])/nrow(base$DMat)), rep(1,times=nrow(base$DMat)))
  expect_equal(.rowSums(base$CMat[base$CMat>0], m=nrow(base$CMat),n=length(base$CMat[base$CMat>0])/nrow(base$CMat)), rep(1,times=nrow(base$DMat)))
  expect_equal(rowSums(base$CMat), rowSums(base$DMat))
})
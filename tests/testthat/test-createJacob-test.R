context("check if the jacobian matrix contains NA")
load("fluxing_webs/fluxweb/data/groups.level.RData")
test_that('all values defined',{

  mat.fluxes = fluxing(groups.level$mat, 
                       groups.level$biomasses, 
                       0.71*groups.level$bodymasses, 
                       groups.level$efficiencies)
  growth.rates = rep(NA, dim(groups.level$mat)[1])
  growth.rates[colSums(groups.level$mat) == 0] = 0.5
  
  jacob = create.jacob(mat.fluxes, groups.level$biomasses, 0.71*groups.level$bodymasses, groups.level$efficiencies, growth.rates)
  expect_equal(0, sum(is.na(jacob)))

}) 

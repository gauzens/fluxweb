context("check if the jacobian matrix contains NA")
load("groups.level.RData")
test_that('all values defined',{

  mat.fluxes = fluxing(groups.level$mat, 
                       groups.level$biomasses, 
                       0.71*groups.level$bodymasses, 
                       groups.level$efficiencies)
  met.types = rep("animal", nrow(mat.fluxes))
  met.types[groups.level$efficiencies == 0.545] = "plant"
  met.types[1] = "detritus" # not true, only for testing purpose
  jacob = create.jacob(mat.fluxes, groups.level$biomasses, groups.level$efficiencies, met.types)
  expect_equal(0, sum(is.na(jacob)))

}) 

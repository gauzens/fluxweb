#' sensitivity analysis
#'
#' Assess how sensitive are results from argument function to variability of input parameter through coefficient of variation.
#'
#' @param fun.name Function to analyse.
#' @param ... Arguments to be passed to \code{fun.name}. Argument names must exactly match those of fun.name.
#' @param param.name Parameter from \code{...} on wich variation is applied.
#' @param var Define the interval of incertainty for the uniform law around \code{x} as \code{[x - x*var, x + x*var]}.
#' @param n Number of replicates.
#' @param full.output logical, if \code{TRUE} all of \code{n} estimations of \code{fun.name} are returned. Only their mean otherwise.
#'
#' @return mean coefficient of variation in comparison to non randomised inputs among all the replicates
#'
#' @details
#'
#' At each replicate, a coefficient of variation is computed (relatively to results obtained form \code{fun.name} without random variation).
#' if \code{full.output} is \code{FALSE} (default) an object of the same type as the one produced by \code{fun.name} is returned, containing all of variation coefficients.
#' If \code{full.output} is \code{TRUE}, a list of size \code{n} with of objects containing variation coeficient is returned.
#'
#' Argument for \code{...} should be passed with their names.
#'
#' 
#' @examples
#'

#' # first compute species per unit biomass metabolic rates using the metabolic theory:
#'losses = 0.1 * species.level$bodymasses^(-0.25)
#'
#'
#'res = sensitivity(fluxing, "mat", 0.1, 5, full.output = TRUE, mat = species.level$mat, biomasses = species.level$biomasses, losses = losses, efficiencies = species.level$efficiencies)
#'res = sensitivity(fluxing, "efficiencies", 0.01, 50, mat = species.level$mat, biomasses = species.level$biomasses, losses = losses, efficiencies = species.level$efficiencies)
#'
#' # growth rates of basal sppecies
#'growth.rate = rep(0.5, length(species.level$biomasses[colSums(species.level$mat) == 0]))
#'
#'val.mat = fluxing(species.level$mat, species.level$biomasses, losses, species.level$efficiencies)
#'#sensitivity(stability.value, "efficiencies", 0.01, 50, val.mat = val.mat, biomasses = species.level$biomasses, losses = losses, efficiencies = species.level$efficiencies, growth.rate = growth.rate)
#'
#'
#'cvs = c()
#'for (var in seq(0, 0.6, 0.05)){
#'  cvs = c(cvs, sensitivity(stability.value, "mat", var, 50, val.mat = val.mat, biomasses = species.level$biomasses, losses = losses, efficiencies = species.level$efficiencies, growth.rate = growth.rate))
#'}
#'
#'plot(abs(cvs) ~ seq(0, 0.6, 0.05))
#'
#' @export
#'


sensitivity = function(fun.name, param.name, var, n, full.output = FALSE, ...){
  is.scalar = function(x){
    return(is.vector(x) && length(x) == 1)
  }
  # get arguments for fun.name
  args = list(...)
  # store result from the function without any randomness
  res.init = do.call(fun.name, args)
  # creates vector of results
  if (full.output){
    if (is.scalar(res.init)){
      res.output = rep(NA, n)
    }else{
      res.output = list()
    }
  }
  vec.cv = rep(NA, n)
  simple.res = 0
  for (i in 1:n){

    ## first, modification of the focal parameter
    if (is.matrix(args[[param.name]])){
      param = args[[param.name]]
      # fill deviation matrix with uniform distribution in [-var, var]
      deviation = matrix(stats::runif(dim(param)[1]*dim(param)[1], min = 1 - var, max = 1 + var), dim(param)[1], dim(param)[1])
      # output = initial value * deviation
      args[[param.name]] = args[[param.name]] * (deviation)
    }

    if (is.vector(args[[param.name]])){ # scalars are also vectors...
      param = args[[param.name]]
      # fill deviation vector with uniform distribution in [-var, var]
      deviation = stats::runif(length(param), min = 1 - var, max = 1 + var)
      # output = initial value * deviation
      args[[param.name]] = args[[param.name]] * (deviation)
      # print(args[[param.name]])
    }

    ## then call to the function
    res = do.call(fun.name, args)
    neutral.res = res.init
    neutral.res[neutral.res == 0] = NA
    # then, for this draw, departure is the coefficient of variation (base = results without random variation)
    # for matrix or vector, I take the mean value of all cv
    res.cv = abs((res - res.init)/neutral.res)


    if (full.output){
      if (is.scalar(res.init)){
        res.output[i] = res.cv
      }else{
        res.output[[i]] = res.cv
      }
    }
    #  sensitivity value is then the means of the coefficient of variation
    simple.res = simple.res + res.cv/n

  }
  if (full.output){
    return(res.output)
  }else{
    return(simple.res)
  }

}

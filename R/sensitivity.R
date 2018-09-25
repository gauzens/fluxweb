#' sensitivity analysis
#'
#' Assesses how sensitive the results from argument function are to variability of input parameter through coefficient of variation.
#'
#' @param fun.name Function to analyse.
#' @param ... Arguments to be passed to \code{fun.name}. Argument names must exactly match those of fun.name.
#' @param param.name Parameter from \code{...} on which variation is applied.
#' @param var Define the interval of uncertainty for the uniform law around \code{x} as \code{[x - x*var, x + x*var]}.
#' @param n Number of replicates.
#' @param full.output Logical, if \code{TRUE} all of \code{n} estimations of \code{fun.name} are returned. Only their mean otherwise.
#'
#' @return a list of two elements of the same type as \code{param.name}: 
#' first element contains the mean coefficient of variation in comparison to non randomised inputs among all the replicates, 
#' second element contains the standard deviation of these coefficient of variation
#'
#' @details
#'
#' At each replicate, a coefficient of variation is computed (relative to results obtained form \code{fun.name} without random variation).
#' if \code{full.output} is \code{FALSE} (default) a list of two objects of the same type as the one produced by \code{fun.name} is returned, 
#' first element contains the mean coefficient of variation in comparison to non randomised inputs among all the replicates, 
#' second element contains the standard deviation of these coefficients of variation
#' If \code{full.output} is \code{TRUE}, a list of size \code{n} with of objects containing the coefficients of variation  is returned.
#'
#' Argument for \code{...} should be passed with their names.
#' 
#' @examples
#'

#' # first compute species per unit biomass metabolic rates using the metabolic theory:
#'losses = 0.1 * species.level$bodymasses^(-0.25)
#'
#'
#'res = sensitivity(fluxing, "mat", 0.1, 5, full.output = TRUE, 
#'                  mat = species.level$mat, 
#'                  biomasses = species.level$biomasses, 
#'                  losses = losses, 
#'                  efficiencies = species.level$efficiencies)
#'res = sensitivity(fluxing, "efficiencies", 0.01, 50, 
#'                  mat = species.level$mat, 
#'                  biomasses = species.level$biomasses, 
#'                  losses = losses, 
#'                  efficiencies = species.level$efficiencies)
#'
#' # growth rates of basal species
#' growth.rates = rep(NA, dim(species.level$mat)[1])
#' growth.rates[colSums(species.level$mat) == 0] = 0.5
#'
#'val.mat = fluxing(species.level$mat, species.level$biomasses, losses, species.level$efficiencies)
#'
#'
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
  
  # two scalars containing the mean coefficient of variation and the sd
  # mean and sd are calculated from the distribution of n replicates.
  simple.res.mean = 0
  simple.res.sd = 0
  # store the initial value of the parameter to change
  param = args[[param.name]]
  for (i in 1:n){
    ## first, modification of the focal parameter
    if (is.matrix(args[[param.name]])){
      # param = args[[param.name]]
      # fill deviation matrix with uniform distribution in [-var, var]
      deviation = matrix(stats::runif(dim(param)[1]*dim(param)[1], min = 1 - var, max = 1 + var), dim(param)[1], dim(param)[1])
      # output = initial value * deviation
      args[[param.name]] = param * (deviation)
    }

    if (is.vector(args[[param.name]])){ # scalars are also vectors...
      # param = args[[param.name]]
      # fill deviation vector with uniform distribution in [-var, var]
      deviation = stats::runif(length(param), min = 1 - var, max = 1 + var)
      # print(deviation)
      # output = initial value * deviation
      args[[param.name]] = param * (deviation)
      # print(rbind(param, deviation, args[[param.name]]))
    }

    ## then call to the function
    res = do.call(fun.name, args)
    neutral.res = res.init
    neutral.res[neutral.res == 0] = NA
    # then, for this draw, departure is the coefficient of variation (base = results without random variation)
    # for matrix or vector, I take the mean value of all cv
    res.cv = (res - res.init)/neutral.res


    if (full.output){
      if (is.scalar(res.init)){
        res.output[i] = res.cv
      }else{
        res.output[[i]] = res.cv
      }
    }
    #  sensitivity value is then the means of the coefficient of variation
    simple.res.mean = simple.res.mean + res.cv/n
    simple.res.sd = simple.res.sd + res.cv*res.cv
  }
  simple.res.sd = sqrt(simple.res.sd/n - simple.res.mean * simple.res.mean)
  if (full.output){
    return(res.output)
  }else{
    return(list(simple.res.mean, simple.res.sd))
  }

}

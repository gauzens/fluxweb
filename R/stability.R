#' Estimates network stability
#'
#' Computes resilience of the system through Jacobian matrix eigenvalues.
#'
#'
#' @param val.mat A matrix describing fluxes between species (usually a result of \code{\link[fluxweb]{fluxing}} function).
#' @param biomasses A vector of species biomasses.
#' @param losses A vector or an array of species energy losses (excluding predation).
#' @param efficiencies A vector or an array of conversion efficiencies of species in the adjacency matrix. These values describe the proportion of consumed energy that is converted to biomass of the consumer.
#' @param growth.rate A vector defining growth rate of basal species.
#' @param bioms.prefs Logical, if \code{TRUE} (default) preferences are scaled according to species biomasses.
#' @param bioms.losses Logical, if \code{TRUE} (default) losses are scaled with biomass.
#' @param ef.level Set to \code{"prey"} if efficiencies are defined by prey, \code{"pred"} if they are a property of the predator.
#' @param full.output Logical, if \code{TRUE} function return supplementary informations.
#'
#' @return Maximum eigenvalue of the Jacobian matrix of a Lotka Voltera like system of equations. If full.output, Jacobian eigenvalues and eigenvectors are returned.
#'
#'
#'@details
#'
#'
#For more information about the set of the underlying system of differential equations and mathematical derivation, please read Gauzens et al. 2017 SI2. <https://www.biorxiv.org/content/early/2017/12/06/229450>
#'
#'
#'\itemize{
#'
#'\item{\code{losses}:} Express species energetic losses not related to consumption. Usually metabolic or death rates.
#'When an array is provided, losses associated to each species correspond to line sums.
#'
#'\item{\code{efficiencies}:} Determines how efficient species are to convert energy (see \code{ef.level} for more details).
#'Providing an array will assume values depending on both prey and predator identity.
#'
#'\item{\code{growth.rate}:} Growth rates of basal species defined. Length of the vector should be equal to the number of species. 
#'expects positive numeric values for index corresponding to basal species, NA otherwise
#'
#'\item{\code{bioms.pref}:} If \code{TRUE}, preferences \eqn{w_{ij}} of predator j on prey i are scaled according to species biomass using the following formula:
#'\deqn{
#'w_{i,j} = \frac{mat[i,j] * biomasses[i]}{\sum_k mat[i,k]* biomasses[k]}
#'}
#'\item{\code{bioms.losses}:} If \code{TRUE}, function will assume that losses are defined per biomass unit.
#'Thus, total losses will be thereafter multiplied by biomass values for each species.
#'
#'\item{\code{ef.level}:} If \code{"prey"} (resp \code{"pred"}), the total amount of energy that can be metabolized from a trophic link
#'will be determined by prey (resp pred) identity. \code{"link.specific"} assumes that efficiencies are defined for each trophic interaction
#'and implies \code{efficiencies} parameter to be a matrix
#'
#'\item{\code{full.output}:} If \code{TRUE}, function result is a list of eigenvalues and eigenvectors of the Jacobian matrix.
#'}
#'
#' @examples
# # first compute species per unit biomass metabolic rates using the metabolic theory:
#' losses = 0.15 * groups.level$bodymasses^(-0.25)
#'
#' # growth rates of basal sppecies
#' growth.rates = rep(NA, dim(groups.level$mat)[1])
#' growth.rates[colSums(groups.level$mat) == 0] = 0.5
#'
#' val.mat = fluxing(groups.level$mat, 
#'                   groups.level$biomasses, 
#'                   losses, 
#'                   groups.level$efficiencies, 
#'                   bioms.pref = TRUE, 
#'                   ef.level = "pred")
#'                   
#' stability.value(val.mat, 
#'                 groups.level$biomasses, 
#'                 losses, 
#'                 groups.level$efficiencies, 
#'                 growth.rates, 
#'                 ef.level = "pred")
#'
#' @author Benoit Gauzens, \email{benoit.gauzens@gmail.com}
#'
#' @export


stability.value = function(val.mat,
                           biomasses,
                           losses,
                           efficiencies,
                           growth.rate,
                           bioms.prefs = TRUE,
                           bioms.losses = TRUE,
                           ef.level = "prey",
                           full.output = FALSE){

  # mat
  if (! is.numeric(val.mat)){
    stop("'val.mat' must be numeric")
  }
  if (dim(val.mat)[1] != dim(val.mat)[2]){
    stop("val.mat should be a square matrix")
  }

  # biomasses
  if (! is.vector(biomasses)){
    stop("biomasses should be a vector")
  } else{
    if (length(biomasses) != dim(val.mat)[1]){
      stop("length of biomasses vector should equal to dimensions of val.mat")
    }
  }
  if (! is.numeric(biomasses)){
    stop("'biomasses' must be numeric")
  } else if (any(biomasses <= 0)){
    stop("'biomasses' must be all >0")
  }

  # losses
  if (! is.numeric(losses)){
    stop("'losses' should be numeric")
  } else if (any(losses < 0)){
    stop("'losses' contain negative value(s)")
  }

  # efficiencies
  if (! is.numeric(efficiencies)){
    stop("'efficiencies' must be numeric")
    if (!(any(efficiencies < 0) || any(efficiencies > 1))) {
      stop("'efficiencies' must be all in interval [0,1]")
    }

    if (is.vector(efficiencies)){
      if (length(efficiencies) != dim(val.mat)[1]){
        stop("'efficiencies' vector length sould be equal to number of species (dimension of val.mat)")
      }
    } else if (dim(efficiencies != dim(val.mat))){
      stop("'efficiencies' matrix dimension different from 'val.mat'")
    }
  }

  if (!(ef.level %in% c('prey', 'pred', 'link.specific'))){
    stop("ef.level should be set to 'pred', 'prey' or 'link.specific'")
  }
  if (ef.level == 'prey' && is.matrix(efficiencies)){ # if user did not change the ef.level = "prey" optional argument but provide a matrix for efficiencies:
    warning("'ef.level' is set to 'prey' and expect a vector of efficiencies but get a matrix instead.\n ef.level was then set to 'link.specific'")
    ef.level = 'link.specific'
  }

  # growth.rate
  if (! is.numeric(growth.rate)){
    stop("growth.rate should be numeric")
  } else if (any(growth.rate[colSums(val.mat) == 0] < 0)){
    stop("all growth.rate values are expected to be positive (or 0)")
  }
  if (!is.vector(growth.rate)){
    stop("growth.rate should be a vector")
  } else if (length(growth.rate) != dim(val.mat)[1]){
    stop("length of growth.rate vector should be equal to the number of species")
  }
  if (sum(!is.na(growth.rate[colSums(val.mat)!=0])) > 0){
    warning('some growth rates are defined for non basal species. Values will be ignored')
  }
  if (sum(is.na(growth.rate[which(colSums(val.mat)==0)])) > 0){
    stop('growth rates of some basal species are not defined')
  }
  
  ### define loss vector as the sum of species losses:
  if (! is.vector(losses)){
    losses = rowSums(losses)
  }
  if (bioms.losses == T){
    losses = losses*biomasses
  }

  nb_s = dim(val.mat)[1]
  nb_b = sum(colSums(val.mat) == 0)
  jacob = matrix(0, nb_s, nb_s)
  jacob = create.jacob(val.mat, biomasses, losses, efficiencies, growth.rate, bioms.prefs = TRUE, bioms.losses = TRUE, ef.level = "prey")

  # return(eigen(jacob))
  if (full.output){
    # should also return the Jacobian
    return(eigen(jacob))
  } else{
    return(max(Re(eigen(jacob)$values)))
  }
}


#' making network stability
#'
#' Find the smallest scalar multiplying a variable from losses insuring system stability
#'
#' @param val.mat A matrix describing fluxes between species (usually a result of \code{\link[fluxweb]{fluxing}} function).
#' @param biomasses A vector of species biomasses.
#' @param losses A vector or an array of species energy losses (excluding predation).
#' @param efficiencies A vector or an array of conversion efficiencies of species in the adjacency matrix. These values describe the proportion of consumed energy that is converted to biomass of the consumer.
#' @param growth.rate A vector defining growth rate of basal species.
#' @param losses.scale Defines a Column from \code{losses} on which scalar multiplication will be tested. (default \code{NULL} if the value is independent of losses).
#' @param bioms.prefs Logical, if \code{TRUE} (default) preferences are scaled accordingly to species biomasses.
#' @param bioms.losses Logical, if \code{TRUE} (default) losses are scaled with biomass.
#' @param ef.level Set to \code{"prey"} if efficiencies are defined by prey, \code{"pred"} if they are a property of the predator.
#' @param interval Search interval for returned value.
#' @param ... Optional parameters for function \code{\link[stats]{uniroot}}
#'
#'
#' @return
#' A list from \code{\link[stats]{uniroot}} function.
#'
#' @details
#' The function assumes a monotonous increase of stability with multiplication by a scalar value. Solution is estimated from the \code{\link[stats]{uniroot}} function, and stability using the \code{\link[fluxweb]{fluxing}} function
#' Thus, accordingly to \code{\link[stats]{uniroot}} solving criteria, if stability values at the two extremum parts of the interval are of same sign, an error is raised.
#'
#' Behavior of the multiplicative term depends on the type of losses:
#' \itemize{
#' \item{\code{losses.scale = NULL} and \code{is.vector(losses)}:} multiplication will be applied to the \code{losses} vector.
#' \item{\code{losses.scale = NULL} and \code{is.matrix(losses)}:} multiplication will be independent of any columns from \code{losses}.
#' \item{\code{losses.scale = FALSE} :} value used for multiplication always independent of losses.
#' \item{other values:} should refer to an element of losses.
#' }
#'
#' @seealso \code{\link[stats]{uniroot}} for root estimate and \code{\link[fluxweb]{stability.value}} for assessing system stability.
#'
#'
#' @examples
#'
# # first compute species per unit biomass metabolic rates using the metabolic theory:
#' losses = 0.15 * groups.level$bodymasses^(-0.25)
#'
#' # growth rates of basal sppecies
#' growth.rates = rep(NA, dim(groups.level$mat)[1])
#' growth.rates[colSums(groups.level$mat) == 0] = 0.5
#'
#' val.mat = fluxing(groups.level$mat, 
#'                   groups.level$biomasses, 
#'                   losses, 
#'                   groups.level$efficiencies, 
#'                   bioms.pref = TRUE, 
#'                   ef.level = "pred")
#' make.stability(val.mat, 
#'                groups.level$biomasses, 
#'                losses, 
#'                groups.level$efficiencies, 
#'                growth.rates, 
#'                ef.level = "pred")
#'
#' @export
#'
#'

make.stability = function(val.mat,
                          biomasses,
                          losses,
                          efficiencies,
                          growth.rate,
                          losses.scale = NULL,
                          bioms.prefs = TRUE,
                          bioms.losses = TRUE,
                          ef.level = "prey",
                          interval = c(1e-12,1),
                          ...){

  stability.value.wrapper = function(x, unvariant, col, val.mat, biomasses, efficiencies, growth.rate, bioms.prefs, bioms.losses, ef.level){
    tot.losses = x*col + unvariant
    return(stability.value(val.mat, biomasses, tot.losses, efficiencies, growth.rate, bioms.prefs, bioms.losses, ef.level))
  }

  # mat
  if (! is.numeric(val.mat)){
    stop("'val.mat' must be numeric")
  }
  if (dim(val.mat)[1] != dim(val.mat)[2]){
    stop("val.mat should be a square matrix")
  }
  nb_s = dim(val.mat)[1]
  # biomasses
  if (! is.vector(biomasses)){
    stop("biomasses should be a vector")
  } else{
    if (length(biomasses) != dim(val.mat)[1]){
      stop("length of biomasses vector should equal to dimensions of val.mat")
    }
  }
  if (! is.numeric(biomasses)){
    stop("'biomasses' must be numeric")
  } else if (any(biomasses <= 0)){
    stop("'biomasses' must be all >0")
  }

  # losses
  if (! is.numeric(losses)){
    stop("'losses' should be numeric")
  } else if (any(losses < 0)){
    stop("'losses' contain negative value(s)")
  }

  # efficiencies
  if (! is.numeric(efficiencies)){
    stop("'efficiencies' must be numeric")
    if (!(any(efficiencies < 0) || any(efficiencies > 1))) {
      stop("'efficiencies' must be all in interval [0,1]")
    }

    if (is.vector(efficiencies)){
      if (length(efficiencies) != dim(val.mat)[1]){
        stop("'efficiencies' vector length sould be equal to number of species (dimension of val.mat)")
      }
    } else if (dim(efficiencies != dim(val.mat))){
      stop("'efficiencies' matrix dimension different from 'val.mat'")
    }
  }

  #ef.level
  if ((ef.level != "prey") && (ef.level != "pred")){
    stop("ef.level should be either 'pred' or 'prey'")
  }

  # growth.rate
  if (! is.numeric(growth.rate)){
    stop("growth.rate should be numeric")
  } else if (any(growth.rate[colSums(val.mat) == 0] < 0)){
    stop("all growth.rate values are expected to be positive (or 0)")
  }
  if (!is.vector(growth.rate)){
    stop("growth.rate should be a vector")
  } else if (length(growth.rate) != dim(val.mat)[1]){
    stop("length of growth.rate vector should be equal to the number of species")
  }
  if (sum(!is.na(growth.rate[colSums(val.mat)!=0])) > 0){
    warning('some growth rates are defined for non basal species. Values will be ignored')
  }
  if (sum(is.na(growth.rate[which(colSums(val.mat)==0)])) > 0){
    stop('growth rates of some basal species are not defined')
  }
  
  # compute losses
  # here, operation is:
  # final.losses = colSums(invariant + value*col), where value is the value to estimate

  if (missing(losses.scale)){
    if (is.vector(losses)){
      unvariant = rep(0, nb_s)
      col = losses
    }
    if (is.matrix(losses)){
      #generate a warning as behaviour will be the same as losses.scale = false?
      unvariant = colSums(losses)
      col = rep[1, nb_s]
    }
  } else{
    if (is.numeric(losses.scale) && !is.integer(losses.scale)){
      stop("losses.scale should be a either character or non integer numeric")
    }
    if (is.vector(losses) && losses.scale == FALSE){
      unvariant = losses
      col = rep(1, nb_s)
    }
    if (is.matrix(losses) && losses.scale == FALSE){
      unvariant = colSums(losses)
      col = rep(1, nb_s)
    }
    if (is.matrix(losses) && is.integer(losses.scale)){
      if ((losses.scale <1 ) || (losses.scale > dim(val.mat[1]))){
        stop("losses.scale should be between 1 and species number")
      }
      col = losses[,losses.scale]
      unvariant = colSums(losses[, -losses.scale])
    }
    if (is.matrix(losses) && is.character(losses.scale)){
      if (!(losses.scale %in% names(losses))){
        stop("losses.scale not a column name of losses")
      }
      col = losses[,names(losses.scale) == losses.scale]
      unvariant = colSums(losses[, colnames(losses) != losses.scale])
    }
  }

  # get stability values to at the interval boundaries
  stab.upper = stability.value.wrapper(max(interval), unvariant, col, val.mat, biomasses, efficiencies, growth.rate, bioms.prefs, bioms.losses, ef.level)
  stab.lower =  stability.value.wrapper(min(interval), unvariant, col, val.mat, biomasses, efficiencies, growth.rate, bioms.prefs, bioms.losses, ef.level)

  if (stab.lower < 0 && stab.upper <0){
    stop("unable to assess minimal stability value: system always stable in specified interval.")
  }
   if (stab.lower > 0 && stab.upper > 0){
    stop("unable to assess minimal stability value: system never stable in specified interval.")
  }

  min.stab.value = stats::uniroot(stability.value.wrapper, interval, ..., unvariant, col, val.mat, biomasses, efficiencies, growth.rate, bioms.prefs, bioms.losses, ef.level,
                           lower = min(interval), upper = max(interval), f.lower = stab.lower, f.upper = stab.upper)
  return(min.stab.value)

}


#' Internal function to compute Jacobian
#'
#'@keywords internal
#'
#'
#'

create.jacob = function(val.mat, biomasses, losses, efficiencies, growth.rate, bioms.prefs = TRUE, bioms.losses = TRUE, ef.level = "prey"){
  nb_s = dim(val.mat)[1]
  nb_b = sum(colSums(val.mat) == 0)
  jacob = matrix(0, nb_s, nb_s)
  basal.index = which(colSums(val.mat) == 0) # store the position of basal species in arguments
  # double loop here... maybe possible to vectorise
    for (i in 1:nb_s){
      for (j in 1: nb_s){
        if (i == j){
          if (i %in% basal.index){ # then species i is basal
            jacob[i,j] = growth.rate[i] - sum(val.mat[i,])/biomasses[i]
          } else{
            if (ef.level == "pred") {
            jacob[i,j] = -losses[i] + val.mat[i,j] * (efficiencies[i] - 1) / biomasses[i] +
              efficiencies[i] * sum(val.mat[,i])/biomasses[i] - sum(val.mat[i,]) / biomasses[i]
            }
            if (ef.level == "prey" && is.vector(efficiencies)){
              jacob[i,j] = -losses[i] + val.mat[i,j] * (efficiencies[i] - 1) / biomasses[i] +
                sum(val.mat[,i] * efficiencies) / biomasses[i] - sum(val.mat[i,]) / biomasses[i]
            }
            if (ef.level == "link.specific"){
              jacob[i,j] = -losses[i] + val.mat[i,j] * (efficiencies[i,i] - 1) / biomasses[i] +
                sum(val.mat[,i] * efficiencies[,i]) / biomasses[i] - sum(val.mat[i,]) / biomasses[i]
            }
          }
        } else{
          if (ef.level == "pred") {jacob[i,j] = (efficiencies[i]*val.mat[j,i] - val.mat[i,j])/biomasses[j]}
          if (ef.level == "prey" && is.vector(efficiencies)) {jacob[i,j] = (efficiencies[j]*val.mat[j,i] - val.mat[i,j])/biomasses[j]}
          # should change the following line to consider a flag for efficiencies = "link.level" instead of checking if it's a matri or not
          if (ef.level == "link.specific") {jacob[i,j] = (efficiencies[j,i]*val.mat[j,i] - val.mat[i,j])/biomasses[j]}
        }
      }
    }
  return(jacob)
}


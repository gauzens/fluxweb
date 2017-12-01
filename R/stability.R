#' estimate network stability
#'
#' compute resiliance of the system through jacobian eigenvalues
#'
#'
#' @param val.mat A matrix describing fluxes between species (usually a result of \code{\link[fluxweb]{fluxing}} function).
#' @param biomasses A vector of species biomasses.
#' @param losses A vector or an array of species energy losses (excluding predation).
#' @param efficiencies A vector or an array of conversion efficiencies of species in the adjacency matrix. These values describe the proportion of consumed energy that is converted to biomass of the consumer.
#' @param growth.rate A vector defining growth rate of basal species.
#' @param bioms.prefs Logical, if \code{TRUE} (default) preferences are scaled accordingly to species biomasses.
#' @param bioms.losses Logical, if \code{TRUE} (default) losses are scaled with biomass.
#' @param ef.level Set to \code{"prey"} if efficiences are defined by prey, \code{"pred"} if they are a property of the predator.
#' @param full.output Logical, if \code{TRUE} function return supplementary informations
#'
#' @return maximum eigenvalue of the jacobian matrix of a lotka Voltera like system of equation If \code{full.output}, Jacobian eigenvalues and eigenvectors are returned.
#'
#'
#'@details
#'
#'
#For more information about the set of the underlying system of differential equations and mathematical derivation, please read REF TO SI of the paper
#'
#'
#'\itemize{
#'
#'\item{\code{losses}:} express species energetic losses not related to predation. Usually metabolic or death rates.
#'When an array is provided, losses associated to each species correspond to line sums.
#'
#'\item{\code{efficiencies}:} Determines how efficient species are to convert energy (see \code{ef.level} for more details).
#'Providing an array will assume values depending on both prey and predator identity
#'
#'\item{\code{growth.rate}:} Growth rates of basal species defined in \code{growth.rate} should appear in the same order as in other arguments.
#'For example the second value specified in \code{growth.rate} should set the groth rate of the second basal species found in \code{biomasses}. IS THAT CLEAR...
#'
#'\item{\code{bioms.pref}:} If \code{TRUE}, preferences \eqn{w_{ij}} of predator j on prey i are scaled accordingly to species biomass unsing the following formula:
#'\deqn{
#'w_{i,j} = \frac{mat[i,j] * biomasses[i]}{\sum_k mat[i,k]* biomasses[k]}
#'}
#'\item{\code{bioms.losses}:} If \code{TRUE}, function will assume that losses are defined per biomass unit.
#'Thus, total losses will be thereafter multiplied by biomass values for each species.
#'
#'\item{\code{ef.level}:} if \code{"prey"} (resp \code{"pred"}), the total amount of energy that can be metabolised from a trophic link
#'will be determined by prey (resp pred) identity. \code{"link.specific"} assumes that efficiencies are defined for each trophic interaction
#'and ask \code{efficiencies} parameter to be a matrix
#'
#'\item{\code{full.output}:} If \code{TRUE}, function result is a list of eigenvalues and eigenvectors of the jacobian matrix. (return also the jacobian?)
#'
#'}
#'
#'
#'
#' @examples
# # first compute species per unit biomass metabolic rates using the metabolic theory:
#' losses = 0.1 * species.level$bodymasses^(-0.25)
#'
#' # growth rates of basal sppecies
#' growth.rates = rep(0.5, length(species.level$biomasses[colSums(species.level$mat) == 0]))
#'
#' val.mat = fluxing(species.level$mat, species.level$biomasses, losses, species.level$efficiencies, bioms.pref = TRUE, ef.level = "pred")
#' stability.value(val.mat, species.level$biomasses, losses, species.level$efficiencies, growth.rates, ef.level = "pred")
#'
#' @author Benoit gauzens, \email{benoit.gauzens@gmail.com}
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

  # efficiences
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
  } else if (any(growth.rate < 0)){
    stop("all growth.rate values are expected to be positive (or 0)")
  }
  if (!is.vector(growth.rate)){
    stop("growth.rate should be a vector")
  } else if (length(growth.rate) != sum(colSums(val.mat) == 0)){
    stop("length of growth.rate vector should be equal to the number of primary producer (species with no prey)")
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
    # should also return the jacobian
    return(eigen(jacob))
  } else{
    return(max(Re(eigen(jacob)$values)))
  }
}


#' making network stability
#'
#' Find the smallest multiplicator of a variable from losses insuring system stability - !Version under devlopment!
#'
#' @param val.mat A matrix describing fluxes between species (usually a result of \code{\link[fluxweb]{fluxing}} function).
#' @param biomasses A vector of species biomasses.
#' @param losses A vector or an array of species energy losses (excluding predation).
#' @param efficiencies A vector or an array of conversion efficiencies of species in the adjacency matrix. These values describe the proportion of consumed energy that is converted to biomass of the consumer.
#' @param growth.rate A vector defining growth rate of basal species.
#' @param losses.scale Defines a Column from \code{losses} mulitplicator value will apply to. (default \code{NULL} if multiplicator independant of losses).
#' @param bioms.prefs Logical, if \code{TRUE} (default) preferences are scaled accordingly to species biomasses.
#' @param bioms.losses Logical, if \code{TRUE} (default) losses are scaled with biomass.
#' @param ef.level Set to \code{"prey"} if efficiences are defined by prey, \code{"pred"} if they are a property of the predator.
#' @param interval Search interval for returned value.
#' @param ... Optional parameters for function \code{\link[stats]{uniroot}}
#'
#'
#' @return
#' A list from \code{\link[stats]{uniroot}} function.
#'
#' @details
#' The function assume a monotonous increase of stability with multiplicator value. Solution is estimated from the \code{\link[stats]{uniroot}} function, and stability using the \code{\link[fluxweb]{fluxing}} function
#' Thus, accordingly to \code{\link[stats]{uniroot}} solving criteria, if stability values at the two extremu parts of the interval are of same sign, an error is raised.
#'
#' Behavior of the multiplicative term depends on the type of losses:
#' \itemize{
#' \item{\code{losses.scale = NULL} and \code{is.vector(losses)}:} multiplicator will be applied to the \code{losses} vector.
#' \item{\code{losses.scale = NULL} and \code{is.matrix(losses)}:} multiplicator will be independant of any columns from \code{losses}.
#' \item{\code{losses.scale = FALSE} :} multiplicator always independant of losses.
#' \item{other values:} should refer to an element of losses.
#'
#' This is an ugly part, but I don't see how to do better right now. Any ideas?
#'
#' }
#'
#' @seealso \code{\link[stats]{uniroot}} for root estimate and \code{\link[fluxweb]{stability.value}} for assessing system stability.
#'
#'
#' @examples
#'
# # first compute species per unit biomass metabolic rates using the metabolic theory:
#' losses = 0.1 * species.level$bodymasses^(-0.25)
#'
#' # growth rates of basal sppecies
#' growth.rates = rep(0.5, length(species.level$biomasses[colSums(species.level$mat) == 0]))
#'
#' val.mat = fluxing(species.level$mat, species.level$biomasses, losses, species.level$efficiencies, bioms.pref = TRUE, ef.level = "pred")
#'make.stability(val.mat, species.level$biomasses, losses, species.level$efficiencies, growth.rates, ef.level = "pred")
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

  # efficiences
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
  } else if (any(growth.rate < 0)){
    stop("all growth.rate values are expected to be positive (or 0)")
  }
  if (!is.vector(growth.rate)){
    stop("growth.rate should be a vector")
  } else if (length(growth.rate) != sum(colSums(val.mat) == 0)){
    stop("length of growth.rate vector should be equal to the number of primary producer (species with no prey)")
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
    stop("unable to assess minimal stability value: system never stable: system never stable in specified interval.")
  }

  min.stab.value = stats::uniroot(stability.value.wrapper, interval, ..., unvariant, col, val.mat, biomasses, efficiencies, growth.rate, bioms.prefs, bioms.losses, ef.level,
                           lower = min(interval), upper = max(interval), f.lower = stab.lower, f.upper = stab.upper)
  return(min.stab.value)

}


#' Internal function to compute jacobian
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
    for (i in 1:nb_s){
      for (j in 1: nb_s){
        if (i == j){
          if (i %in% basal.index){ # then species i is basal
            jacob[i,j] = growth.rate[which(basal.index == i)] - sum(val.mat[i,])/biomasses[i]
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


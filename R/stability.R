#' Estimates network stability
#'
#' Computes resilience of the system through Jacobian matrix eigenvalues.
#'
#'
#' @param val.mat A matrix describing fluxes between species (usually a result of \code{\link[fluxweb]{fluxing}} function).
#' @param biomasses A vector of species biomasses.
#' @param efficiencies A vector or an array of conversion efficiencies of species in the adjacency matrix. These values describe the proportion of consumed energy that is converted to biomass of the consumer.
#' @param metabolic.types A vector containing information on species type (\code{"detritus"}, \code{"plant"} or \code{"animal"})
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
#'
#'\item{\code{efficiencies}:} Determines how efficient species are to convert energy (see \code{ef.level} for more details).
#'Providing an array will assume values depending on both prey and predator identity.
#'
#'
#'\item{\code{ef.level}:} If \code{"prey"} (resp \code{"pred"}), the total amount of energy that can be metabolized from a trophic link
#'will be determined by prey (resp pred) identity. \code{"link.specific"} assumes that efficiencies are defined for each trophic interaction
#'and implies \code{efficiencies} parameter to be a matrix
#'
#'\item{\code{full.output}:} If \code{TRUE}, function result is a list of eigenvalues and eigenvectors of the Jacobian matrix.
#'}
#'
#' @examples
# # first generate some species per unit biomass metabolic rates:
#' losses = 0.15 * groups.level$bodymasses^(-0.25)
#'# define metbolic types:
#' met.types = rep('animal', length(losses))
#' met.types[groups.level$efficiencies == 0.545] = 'plant'
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
#'                 groups.level$efficiencies,
#'                 metabolic.types = met.types,
#'                 ef.level = "pred")
#'
#' @author Benoit Gauzens, \email{benoit.gauzens@gmail.com}
#'
#' @export


stability.value = function(val.mat,
                           biomasses,
                           efficiencies,
                           metabolic.types,
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

  nb_s = dim(val.mat)[1]
  nb_b = sum(colSums(val.mat) == 0)
  jacob = matrix(0, nb_s, nb_s)
  jacob = create.jacob(val.mat, biomasses, efficiencies, metabolic.types, ef.level = "prey")

  # return(eigen(jacob))
  if (full.output){
    # should also return the Jacobian
    return(eigen(jacob))
  } else{
    return(max(Re(eigen(jacob)$values)))
  }
}


#' Compute the Jacobian matrix
#' 
#' @param val.mat A matrix describing fluxes between species (usually a result of \code{\link[fluxweb]{fluxing}} function).
#' @param biomasses A vector of species biomasses.
#' @param efficiencies A vector or an array of conversion efficiencies of species in the adjacency matrix.
#' @param metabolic.types A vector containing information on species type (\code{"detritus"}, \code{"plant"} or \code{"animal"})
#' @param ef.level Set to \code{"prey"} if efficiencies are defined by prey, \code{"pred"} if they are a property of the predator, \code{link.specific} if they depend on both consumer and resource identity.
#'
#' @return A matrix representing the Jacobian matrix of the dynamical system associated with the fluxes from the \code{\link[fluxweb]{fluxing}} function.
#' 
#' 
#' @examples
#' # First compute species per unit biomass metabolic rates:
#' losses = 0.15 * groups.level$bodymasses^(-0.25)
#'
#'
#' val.mat = fluxing(groups.level$mat, 
#'                   groups.level$biomasses, 
#'                   losses, 
#'                   groups.level$efficiencies, 
#'                   bioms.pref = TRUE, 
#'                   ef.level = "prey")
#'                   
#' # define metabolic types                   
#' 
#' met.types = rep("animal", nrow(val.mat))
#' met.types[groups.level$efficiencies == 0.545] = "plant"
#' met.types[groups.level$efficiencies == 0.158] = "detritus"
#' 
#' create.jacob(val.mat, 
#'                groups.level$biomasses, 
#'                groups.level$efficiencies, 
#'                met.types,
#'                ef.level = "prey")
#'                
#' 
#' @author Benoit Gauzens, \email{benoit.gauzens@gmail.com}
#' @export
#'                
#'                

create.jacob = function(val.mat, biomasses, efficiencies, metabolic.types, ef.level = "prey"){
  nb_s = dim(val.mat)[1]
  nb_b = sum(colSums(val.mat) == 0)
  jacob = matrix(0, nb_s, nb_s)

  plants = metabolic.types == "plant"
  detritus = metabolic.types == "detritus"
  
  if (ef.level == "pred"){
    jacob = t(val.mat) * efficiencies - val.mat # ei * Fji - Fij  
    jacob = sweep(jacob, MARGIN = 2, biomasses, '/')
    diag(jacob) = (diag(val.mat) / biomasses) * (efficiencies - 1)
  }
  if (ef.level == "prey"){
    jacob = sweep(t(val.mat), MARGIN = 2, efficiencies, '*') - val.mat # ej * Fji - Fij 
    jacob = sweep(jacob, MARGIN = 2, biomasses, '/')
    diag(jacob) = (diag(val.mat) / biomasses) * (efficiencies - 1)
  }
  if (ef.level == "link.specific"){
    jacob = t(val.mat) * efficiencies - val.mat
    jacob = sweep(jacob, MARGIN = 2, biomasses, '/')
    diag(jacob) = (diag(val.mat) / biomasses) * (diag(efficiencies) - 1)
  }
  diag(jacob)[plants] = 0
  diag(jacob)[detritus] = rowSums(val.mat[detritus,, drop = FALSE]) / biomasses[detritus]
  
  return(jacob)
}


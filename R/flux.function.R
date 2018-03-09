#' generate fluxes
#'
#' Creates a valuated graph adjacency matrix from its binary version.
#'
#' @param mat Network adjacency matrix describing interactions among species. Interactions can be either binary or weighted.
#' @param losses A vector or an array of species energy losses (excluding consumption).
#' @param biomasses Vector of species biomasses.
#' @param efficiencies A vector or an array of conversion efficiencies of species in the adjacency matrix. These values describe the proportion of consumed energy that is converted to biomass of the consumer.
#' @param bioms.prefs Logical - if \code{TRUE}, consumer preferences are scaled according to species biomasses.
#' @param bioms.losses Logical - if \code{TRUE}, losses are scaled with species biomasses.
#' @param ef.level Set to \code{"prey"} if efficiencies are defined by prey, \code{"pred"} if they are a property of the predator.
#'
#'
#' @return Returns an adjacency matrix where entries are the computed energy fluxes between consumer species and their respective resources.
#'
#'
#'@details
#'This function computes fluxes in food webs based on an equilibrium hypothesis: for each species, sum of ingoing fluxes (gains from predation) balances the sum of outgoing fluxes.
#'Outgoing fluxes are defined by consumption and the \code{losses} argument. Usually \code{losses} relate to species metabolic rates and/or natural death rates. For each species \code{i}, sum of ingoing fluxes \code{F_i} is computed as:
#'\deqn{
#'F_{i} = \frac{1}{e_i} (L_i + \sum_j W_{ij}F_j) \quad if \quad \code{ef.level == "pred"}
#'}
#'\deqn{
#'F_{i} = \frac{L_i + \sum_j W_{ij}F_j}{\sum_j W_{ji}e_j} \quad if \quad \code{ef.level == "pred"}
#'}
#'\code{W} set the matrix of preferences estimated from \code{mat}, according to \code{bioms.prefs}. \code{L} is the vector depicting sum of losses
#'(scaled or not by biomasses, accordingly to \code{bioms.losses}) and \code{e} is the vector of species efficiencies.
#'
#'
#'\itemize{
#'\item{\code{mat}:} Either a binary or a valuated matrix can be used. A non zero value for mat[i,j] means that species i is consumed by species j.
#'Matrix entries would assess predator preferences on its prey, thus providing a binary matrix assumes no preferences.
#'
#'\item{\code{losses}:} Express species energetic losses not related to consumption. Usually metabolic or death rates.
#'When an array is provided, losses associated to each species correspond to line sums.
#'
#'\item{\code{efficiencies}:} Determines how efficient species are to convert energy (see \code{ef.level} for more details).
#'Providing an array will assume values depending on both prey and predator identity.
#'
#'\item{\code{bioms.pref}:} If \code{TRUE}, preferences \eqn{W_{ij}} of predator j on prey i are scaled accordingly to species biomass using the following formula:
#'\deqn{
#'W_{i,j} = \frac{mat[i,j] * biomasses[i]}{\sum_k mat[i,k]* biomasses[k]}
#'}
#'If \code{FALSE}, a normalisation on column values is performed.
#'
#'\item{\code{bioms.losses}:} Set to true, function will assume that losses are defined per biomass unit.
#'Thus, total losses will be thereafter multiplied by biomass values for each species.
#'
#'\item{\code{ef.level}:} If \code{"prey"} (resp \code{"pred"}), the total amount of energy that can be metabolized from a trophic link
#'will be determined by prey (resp predator) identity. \code{"link.specific"} assumes that efficiencies are defined for each trophic interaction
#'and implies \code{efficiencies} parameter to be a matrix.
#'
#'}
#'
#'
#' @examples
#' # first compute species per unit biomass metabolic rates using the metabolic theory:
#'losses = 0.1 * species.level$bodymasses^(-0.25)
#'
#'# call of the function:
#'fluxing(species.level$mat, 
#'        species.level$biomasses, 
#'        losses, 
#'        species.level$efficiencies, 
#'        bioms.pref = TRUE, 
#'        ef.level = "prey")
#'
#' @export
#'
#'
#' @author Benoit gauzens, \email{benoit.gauzens@gmail.com}


fluxing = function(mat, biomasses = NULL, losses, efficiencies, bioms.prefs = TRUE, bioms.losses = TRUE, ef.level = "prey"){

  # mat
  if (! is.numeric(mat)){
    stop("'mat' must be numeric")
  }
  if (dim(mat)[1] != dim(mat)[2]){
    stop("mat should be a square matrix")
  }

  # biomasses
  if (!is.null(biomasses)){
    if (! is.vector(biomasses)){
      stop("biomasses should be a vector")
    } else{
      if (length(biomasses) != dim(mat)[1]){
        stop("length of biomasses vector should equal to dimensions of mat")
      }
    }
    if (! is.numeric(biomasses)){
      stop("'biomasses' must be numeric")
    } else if (any(biomasses < 0)){
      stop("'biomasses' must be all >=0")
    }
  } else if (bioms.prefs){
    stop("bioms.prefs set to TRUE but no biomasses provided")
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
      if (length(efficiencies) != dim(mat)[1]){
        stop("'efficiencies' vector length sould be equal to number of species (dimension of mat)")
      }
    } else if (dim(efficiencies != dim(mat))){
      stop("'efficiencies' matrix dimension different from 'mat'")
    }
  }

  if (!(ef.level %in% c('prey', 'pred', 'link.specific'))){
    stop("ef.level should be set to 'pred', 'prey' or 'link.specific'")
  }
  if (ef.level == 'prey' && is.matrix(efficiencies)){ # if user did not change the ef.level = "prey" optional argument but provide a matrix for efficiencies:
    warning("'ef.level' is set to 'prey' and expect a vector of efficiencies but get a matrix instead.\n ef.level was then set to 'link.specific'")
    ef.level = 'link.specific'
  }

  ### first arrange mat: apply the biomass scaling of preferences if needed
  ### columns should sum to 1 for predators, 0 to preys
  column.sum = colSums(mat)
  # in the following, the as.matrix() is needed becaue sometimes mat[, column.sum > 0] is only one column and, thanks to R, is automatically casted to a vector
  if (bioms.prefs){
    # apply 'functional response' of preferencs
    # mat[, column.sum > 0] = apply(mat[, column.sum > 0], 2, function(vec, bioms) vec*biomasses/sum(vec*biomasses), biomasses) #! in the function I should use bioms instead biomasses
    mat[, column.sum > 0] = apply(as.matrix(mat[, column.sum > 0]), 2, function(vec) vec*biomasses/sum(vec*biomasses)) #! in the function biomasses is already defined more globaly, so no need of another parameter

      } else { # here optimise with else if not all element of colsums are equal to either 1 or 0...
    # sum of entries have to sum to one for each predator (normalisaton of preferences)
    colomn.sum = colSums(mat)
    mat[, colomn.sum>0] = sweep(apply(as.matrix(mat[, colomn.sum>0])), 2, colomn.sum[colomn.sum>0], "/")
  }

  ### define loss vector as the sum of species losses:
  if (! is.vector(losses)){ # this is for allowing user to input a loss matrix (different kinds of physiological loss in the same parameter)
    losses = rowSums(losses)
  }
  if (bioms.losses == T){
    losses = losses*biomasses
  }

  ### then solving the system
  # warning here: even if efficiencies are defined at the predator level I need a vector of legth = to number of species (with some arbitrary values for basal species)
  if (ef.level == "pred"){
    F = solve(diag(efficiencies) - mat) %*% losses
  }
  if (ef.level == "prey"){
    vec.in = as.vector(t(mat) %*% efficiencies)
    vec.1p = rep(0, dim(mat)[1])
    vec.1p[colSums(mat) == 0] = 1
    F = solve(diag(vec.in + vec.1p) - mat) %*% losses
  }
  if (ef.level == "link.specific"){
    U = mat * efficiencies
    vec.one = rep(1, dim(efficiencies)[1])
    vec.1p = rep(0, dim(mat)[1])
    vec.1p[colSums(mat) == 0] = 1
    vec.in = as.vector(t(U)%*%vec.one + vec.1p)
    F = solve(diag(vec.in) - mat) %*% losses
  }
  
  if (any(F < 0)){
    stop("model chosen is unable to determine fluxes accoringly to data")
  }
  ### set individual fluxes (each element of ith line from mat.norm is multiplied bu the ith element of F)
  flux.mat = sweep(mat, 2, F, "*")
  return(flux.mat)
}



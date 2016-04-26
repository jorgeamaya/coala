#' @importFrom R6 R6Class
variation_par_class <- R6Class("variation_par", inherit = parameter_class,
  private = list(
    base_par = list(),
    func = "variation",
    add_parameter = function(parameter) {
      if (is.numeric(parameter) && length(parameter) == 1) {
        expr <- parameter
      } else if (is.character(parameter) && length(parameter) == 1) {
        expr <- parameter
      } else if (is.par(parameter)) {
        private$base_par[[length(private$base_par) + 1]] <- parameter
        expr <- parameter$get_expression()
      } else {
        stop("Unexpected type of parameter")
      }
      expr
    }),
  public = list(
    initialize = function(parameter, variance, n_groups) {
      expr_mean <- private$add_parameter(parameter)
      expr_var <- private$add_parameter(variance)
      private$expr <- parse(text = paste0(private$func, "(", expr_mean,
                                          ", ", expr_var,
                                          ", ", n_groups, ")"))
    },
    get_base_par = function() private$base_par
  )
)

#' @importFrom stats rgamma
variation <- function(mean, variance, n_groups = NULL) {
  gamma_shape <-  mean ^ 2 / variance
  gamma_rate <- mean / variance
  if (is.null(n_groups)) return(rgamma(1, gamma_shape, gamma_rate))

  group_median_probs <- (1:n_groups * 2 - 1) / (2 * n_groups)
  group_medians <- qgamma(group_median_probs, gamma_shape, gamma_rate)
  group <- sample.int(n_groups, 1)
  return(group_medians[group] * mean * n_groups / sum(group_medians))


  #probs <- c(group - 1, group) / n_groups
  #probs <- mean(probs)
  probs[probs == 1.0] <- .9999
  interval <- qgamma(probs, gamma_shape, gamma_rate)
  1 / (n_groups * diff(interval))
  #interval
}

is.par_variation <- function(object) inherits(object, "variation_par")


#' Variable Parameters
#'
#' This function can be used to let the values of a parameter vary between
#' the different loci. When used, the values for the enclosed parameter
#' will follow a gamma distribution with mean of the parameters original
#' value, and the variance specified as argument \code{variance}. This requires
#' that the original value is positive. When using this, the simulators
#' are called separately for each locus, which can dramatically increase the
#' time needed to simulate models with many loci.
#'
#' @param par A parameter whichs value will be made variable between the loci.
#' @param variance The variance of the gamma distribution, which the values used
#'   for simulation will follow.
#' @export
#' @seealso For parameters that are identical for all loci: \code{\link{parameter}}
#' @examples
#' model <- coal_model(5, 10) +
#'   feat_mutation(par_variation(par_const(5), 10)) +
#'   sumstat_nucleotide_div()
#' simulate(model)
par_variation <- function(par, variance, n_groups = NULL) {
  variation_par_class$new(par, variance, n_groups)
}

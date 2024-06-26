#' Compute Conditional Probability of Observed Mediator Given True Mediator, for Every Subject
#'
#'  Compute the conditional probability of observing mediator \eqn{M^* \in \{1, 2 \}} given
#'  the latent true mediator \eqn{M \in \{1, 2 \}} as
#'  \eqn{\frac{\text{exp}\{\gamma_{kj0} + \gamma_{kjZ} Z_i\}}{1 + \text{exp}\{\gamma_{kj0} + \gamma_{kjZ} Z_i\}}}
#'  for each of the \eqn{i = 1, \dots,} \code{n} subjects.
#'
#' @param gamma_matrix A numeric matrix of estimated regression parameters for the
#'   observation mechanism, \code{M* | M} (observed mediator, given the true mediator)
#'   ~ \code{Z} (misclassification predictor matrix). Rows of the matrix
#'   correspond to parameters for the \code{M* = 1} observed mediator, with the
#'   dimensions of \code{z_matrix}. Columns of the matrix correspond to the true
#'   mediator categories \eqn{j = 1, \dots,} \code{n_cat}. The matrix should be
#'   obtained by \code{COMMA_EM}, \code{COMMA_PVW}, or \code{COMMA_OLS}.
#' @param z_matrix A numeric matrix of covariates in the observation mechanism.
#'   \code{z_matrix} should not contain an intercept.
#'
#' @return \code{misclassification_prob} returns a dataframe containing four columns.
#'   The first column, \code{Subject}, represents the subject ID, from \eqn{1} to \code{n},
#'   where \code{n} is the sample size, or equivalently, the number of rows in \code{z_matrix}.
#'   The second column, \code{M}, represents a true, latent mediator category \eqn{M \in \{1, 2 \}}.
#'   The third column, \code{Mstar}, represents an observed outcome category \eqn{M^* \in \{1, 2 \}}.
#'   The last column, \code{Probability}, is the value of the equation
#'   \eqn{\frac{\text{exp}\{\gamma_{kj0} + \gamma_{kjZ} Z_i\}}{1 + \text{exp}\{\gamma_{kj0} + \gamma_{kjZ} Z_i\}}}
#'   computed for each subject, observed mediator category, and true, latent mediator category.
#'
#' @include pistar_compute.R
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' sample_size <- 1000
#' cov1 <- rnorm(sample_size)
#' cov2 <- rnorm(sample_size, 1, 2)
#' z_matrix <- matrix(c(cov1, cov2), nrow = sample_size, byrow = FALSE)
#' estimated_gammas <- matrix(c(1, -1, .5, .2, -.6, 1.5), ncol = 2)
#' P_Ystar_Y <- misclassification_prob(estimated_gammas, z_matrix)
#' head(P_Ystar_Y)
misclassification_prob <- function(gamma_matrix,
                                   z_matrix){

  n_cat = 2
  sample_size = nrow(z_matrix)

  if (is.data.frame(z_matrix))
    z_matrix <- as.matrix(z_matrix)
  if (!is.numeric(z_matrix))
    stop("'z_matrix' should be a numeric matrix.")

  if (is.vector(z_matrix))
    z_matrix <- as.matrix(z_matrix)
  if (!is.matrix(z_matrix))
    stop("'z_matrix' should be a matrix or data.frame.")

  Z = matrix(c(rep(1, sample_size), c(z_matrix)),
             byrow = FALSE, nrow = sample_size)

  subject = rep(1:sample_size, n_cat * n_cat)
  Y_categories = rep(1:n_cat, each = sample_size * n_cat)
  Ystar_categories = rep(c(1:n_cat, 1:n_cat), each = sample_size)
  pistar_matrix = pistar_compute(gamma_matrix, Z, sample_size, n_cat)
  pistar_df = data.frame(Subject = subject,
                         Y = Y_categories,
                         Ystar = Ystar_categories,
                         Probability = c(pistar_matrix))

  return(pistar_df)
}

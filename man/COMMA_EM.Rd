% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COMMA_EM.R
\name{COMMA_EM}
\alias{COMMA_EM}
\title{EM Algorithm Estimation of the Binary Mediator Misclassification Model}
\usage{
COMMA_EM(
  Mstar,
  outcome,
  outcome_distribution,
  interaction_indicator,
  x_matrix,
  z_matrix,
  c_matrix,
  beta_start,
  gamma_start,
  theta_start,
  sigma_start = NULL,
  tolerance = 1e-07,
  max_em_iterations = 1500,
  em_method = "squarem"
)
}
\arguments{
\item{Mstar}{A numeric vector of indicator variables (1, 2) for the observed
mediator \code{M*}. There should be no \code{NA} terms. The reference category is 2.}

\item{outcome}{A vector containing the outcome variables of interest. There
should be no \code{NA} terms.}

\item{outcome_distribution}{A character string specifying the distribution of
the outcome variable. Options are \code{"Bernoulli"}, \code{"Normal"}, or
\code{"Poisson"}.}

\item{interaction_indicator}{A logical value indicating if an interaction between
\code{x} and \code{m} should be used to generate the outcome variable, \code{y}.}

\item{x_matrix}{A numeric matrix of predictors in the true mediator and outcome mechanisms.
\code{x_matrix} should not contain an intercept and no values should be \code{NA}.}

\item{z_matrix}{A numeric matrix of covariates in the observation mechanism.
\code{z_matrix} should not contain an intercept and no values should be \code{NA}.}

\item{c_matrix}{A numeric matrix of covariates in the true mediator and outcome mechanisms.
\code{c_matrix} should not contain an intercept and no values should be \code{NA}.}

\item{beta_start}{A numeric vector or column matrix of starting values for the \eqn{\beta}
parameters in the true mediator mechanism. The number of elements in \code{beta_start}
should be equal to the number of columns of \code{x_matrix} and \code{c_matrix} plus 1.
Starting values should be provided in the following order: intercept, slope
coefficient for the \code{x_matrix} term, slope coefficient for first column
of the \code{c_matrix}, ..., slope coefficient for the final column of the \code{c_matrix}.}

\item{gamma_start}{A numeric vector or matrix of starting values for the \eqn{\gamma}
parameters in the observation mechanism. In matrix form, the \code{gamma_start} matrix rows
correspond to parameters for the \code{M* = 1}
observed mediator, with the dimensions of \code{z_matrix} plus 1, and the
gamma parameter matrix columns correspond to the true mediator categories
\eqn{M \in \{1, 2\}}. A numeric vector for \code{gamma_start} is
obtained by concatenating the gamma matrix, i.e. \code{gamma_start <- c(gamma_matrix)}.
Starting values should be provided in the following order within each column:
intercept, slope coefficient for first column of the \code{z_matrix}, ...,
slope coefficient for the final column of the \code{z_matrix}.}

\item{theta_start}{A numeric vector or column matrix of starting values for the \eqn{\theta}
parameters in the outcome mechanism. The number of elements in \code{theta_start}
should be equal to the number of columns of \code{x_matrix} and \code{c_matrix} plus 2
(if \code{interaction_indicator} is \code{FALSE}) or 3 (if
\code{interaction_indicator} is \code{TRUE}). Starting values should be
provided in the following order: intercept, slope coefficient for the \code{x_matrix} term,
slope coefficient for the mediator \code{m} term,
slope coefficient for first column of the \code{c_matrix}, ...,
slope coefficient for the final column of the \code{c_matrix},
and, optionally, slope coefficient for \code{xm}).}

\item{sigma_start}{A numeric value specifying the starting value for the
standard deviation. This value is only required if \code{outcome_distribution}
is \code{"Normal"}. Otherwise, this value is set to \code{NULL}.}

\item{tolerance}{A numeric value specifying when to stop estimation, based on
the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.}

\item{max_em_iterations}{A numeric value specifying when to stop estimation, based on
the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.}

\item{em_method}{A character string specifying which EM algorithm will be applied.
Options are \code{"em"}, \code{"squarem"}, or \code{"pem"}. The default and
recommended option is \code{"squarem"}.}
}
\value{
\code{COMMA_EM} returns a data frame containing four columns. The first
column, \code{Parameter}, represents a unique parameter value for each row.
The next column contains the parameter \code{Estimates}, followed by the standard
error estimates, \code{SE}. The final column, \code{Convergence}, reports
whether or not the algorithm converged for a given parameter estimate.
}
\description{
Jointly estimate \eqn{\beta}, \eqn{\gamma}, and \eqn{\theta} parameters from
the true mediator, observed mediator, and outcome mechanisms, respectively,
in a binary mediator misclassification model.
}
\examples{
set.seed(20240709)
sample_size <- 2000

n_cat <- 2 # Number of categories in the binary mediator

# Data generation settings
x_mu <- 0
x_sigma <- 1
z_shape <- 1
c_shape <- 1

# True parameter values (gamma terms set the misclassification rate)
true_beta <- matrix(c(1, -2, .5), ncol = 1)
true_gamma <- matrix(c(1, 1, -.5, -1.5), nrow = 2, byrow = FALSE)
true_theta <- matrix(c(1, 1.5, -2, -.2), ncol = 1)

example_data <- COMMA_data(sample_size, x_mu, x_sigma, z_shape, c_shape,
                           interaction_indicator = FALSE,
                           outcome_distribution = "Bernoulli",
                           true_beta, true_gamma, true_theta)
                           
beta_start <- matrix(rep(1, 3), ncol = 1)
gamma_start <- matrix(rep(1, 4), nrow = 2, ncol = 2)
theta_start <- matrix(rep(1, 4), ncol = 1)

Mstar = example_data[["obs_mediator"]]
outcome = example_data[["outcome"]]
x_matrix = example_data[["x"]]
z_matrix = example_data[["z"]]
c_matrix = example_data[["c"]]
                           
EM_results <- COMMA_EM(Mstar, outcome, "Bernoulli", FALSE,
                       x_matrix, z_matrix, c_matrix,
                       beta_start, gamma_start, theta_start)

EM_results

}

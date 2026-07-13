#' infometrics: Information-Theoretic Methods for Econometric Estimation
#'
#' Implements Information-Theoretic (IT) estimators for econometric models
#' following the unified framework of Golan (2008). Provides Generalized
#' Maximum Entropy (GME) and Generalized Cross-Entropy (GCE) estimators for
#' linear regression (\code{\link{linreg}}), instrumental variables
#' (\code{\link{linreg_iv}}), panel data (\code{\link{panel_gce}}), multinomial
#' response (\code{\link{multinomial_gce}}, \code{\link{mixed_gce}}), matrix
#' balancing (\code{\link{matrix_ce}}, \code{\link{matrix_gce}}), Markov
#' transition matrices (\code{\link{markov_ce}}, \code{\link{markov_gce}}), and
#' pure/noisy inverse problems (\code{\link{inverse_ce}},
#' \code{\link{inverse_noise}}), all via the concentrated (dual) formulation.
#'
#' @references
#' Golan, A. (2008). Information and Entropy Econometrics — A Review and
#' Synthesis. \emph{Foundations and Trends in Econometrics}, \strong{2}(1-2),
#' 1-145.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats optim pchisq quantile sd fitted terms model.frame model.matrix model.response .getXlevels
NULL

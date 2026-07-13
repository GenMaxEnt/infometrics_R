#' infometrics: Information-Theoretic Methods for Econometric Estimation
#'
#' Implements Information-Theoretic (IT) estimators for econometric models
#' following the unified framework of Golan (2008). Provides entropy measures
#' and the ME, CE, GME, and GCE estimators in both primal and dual form,
#' plus multinomial extensions via \code{\link{me_mnl}} and
#' \code{\link{gme_mnl}}.
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

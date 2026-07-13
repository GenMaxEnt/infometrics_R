# ============================================================
# gce_table.R — side-by-side comparison table for gce_clogit fits
#
# No external dependencies. Produces a formatted text table
# analogous to stargazer but built entirely from the fields
# already computed by gce_clogit() and its S3 methods.
# ============================================================


#' Side-by-side comparison table for \code{gce_clogit} models
#'
#' @description
#' Prints a formatted text table comparing coefficients and model diagnostics
#' across one or more \code{\link{gce_clogit}} fits.  Analogous to
#' \code{stargazer} output but built without external dependencies, directly
#' from the fields produced by \code{gce_clogit()}.
#'
#' Coefficients missing from a particular model are shown as blank cells so
#' the table aligns correctly when models have different regressors (e.g.,
#' income-only vs.\ price+catch vs.\ full mixed model).
#'
#' @param ... One or more objects of class \code{"infometrics_clogit"} returned
#'   by \code{\link{gce_clogit}}.
#' @param column.labels Character vector of column headers, one per model.
#'   Defaults to \code{"(1)"}, \code{"(2)"}, \ldots
#' @param dep.var.label Label for the dependent variable row.
#'   Default: \code{"Dependent variable"}.
#' @param title Optional character string printed as a centred title above the
#'   table.
#' @param digits Integer. Decimal places for coefficients and standard errors.
#'   Default \code{4}.
#' @param stars Logical. Append significance stars to estimates
#'   (\code{*} p<0.10, \code{**} p<0.05, \code{***} p<0.01).
#'   Default \code{TRUE}.
#' @param out Optional file path. When supplied, the table is written to that
#'   file (one line per row) in addition to being printed on the console.
#'
#' @return Invisibly returns the character vector of table lines so callers can
#'   post-process or write to a file.
#'
#' @details
#' The coefficient block shows, for each parameter name found in any model:
#' \itemize{
#'   \item Row 1: point estimate with significance stars (p-values from the
#'         standard normal, using the analytic Hessian SEs from
#'         \code{gce_clogit}).
#'   \item Row 2: standard error in parentheses.
#'   \item Blank cells for models that do not include that parameter.
#' }
#'
#' The diagnostics block reports (Golan 2008, Sections 6--7):
#' \describe{
#'   \item{\bold{Log-likelihood}}{Value of \eqn{\ell = -L_\text{GCE}} at
#'         convergence.}
#'   \item{\bold{McFadden R²}}{\eqn{1 - \ell_{CL}/\ell_0} where
#'         \eqn{\ell_0 = N\log(1/J)}.}
#'   \item{\bold{S(p\_hat)}}{Normalised signal entropy \eqn{\in [0,1]};
#'         1 = uniform (no information), 0 = degenerate.}
#'   \item{\bold{S(w\_hat)}}{Normalised error entropy \eqn{\in [0,1]};
#'         1 = errors fully spread over support, 0 = point-mass at zero.}
#'   \item{\bold{Misses}}{Number of individuals whose predicted most-probable
#'         alternative differs from their observed choice.}
#'   \item{\bold{Support v}}{Range \eqn{[v_1, v_M]} of the error support
#'         vector.}
#' }
#'
#' @references
#' Golan, A. (2008). Information and entropy econometrics — A review and
#' synthesis. \emph{Foundations and Trends in Econometrics}, 2(1--2), 1--145.
#'
#' @examples
#' set.seed(42)
#' N <- 60; J <- 3
#' df <- data.frame(
#'   id     = rep(seq_len(N), each = J),
#'   alt    = rep(letters[seq_len(J)], N),
#'   price  = runif(N * J),
#'   catch  = runif(N * J),
#'   income = rep(runif(N), each = J),
#'   chosen = FALSE
#' )
#' for (i in seq_len(N)) {
#'   rows <- which(df$id == i)
#'   df$chosen[rows[sample(J, 1L)]] <- TRUE
#' }
#' fit1 <- gce_clogit(chosen ~ price,         data = df, alt = "alt", id = "id")
#' fit2 <- gce_clogit(chosen ~ price + catch,  data = df, alt = "alt", id = "id")
#' gce_table(fit1, fit2,
#'           column.labels  = c("Price only", "Price + Catch"),
#'           dep.var.label  = "Fish choice",
#'           title          = "GCE Conditional Logit — Fishing Data")
#'
#' @export
gce_table <- function(...,
                      column.labels = NULL,
                      dep.var.label = "Dependent variable",
                      title         = NULL,
                      digits        = 4L,
                      stars         = TRUE,
                      out           = NULL) {

  fits <- list(...)

  # ---- validate inputs ------------------------------------------------------
  if (length(fits) == 0L)
    stop("Supply at least one infometrics_clogit object.")
  ok <- vapply(fits, inherits, logical(1L), "infometrics_clogit")
  if (!all(ok))
    stop("All positional arguments must be infometrics_clogit objects ",
         "returned by gce_clogit().")

  n_m    <- length(fits)
  digits <- as.integer(digits)

  if (is.null(column.labels))
    column.labels <- paste0("(", seq_len(n_m), ")")
  if (length(column.labels) != n_m)
    stop("column.labels must have length ", n_m, " (one per model).")

  # ---- union of coefficient names, preserving first-appearance order --------
  all_coefs <- unique(unlist(lapply(fits, function(f) names(f$beta))))

  # ---- helper: significance stars ------------------------------------------
  .st <- function(p) {
    if (!stars || is.na(p) || is.nan(p)) return("")
    if (p < 0.01) "***" else if (p < 0.05) "**" else if (p < 0.10) "*" else ""
  }

  # ---- helper: McFadden pseudo-R² ------------------------------------------
  .mcf <- function(f) {
    lp <- log(pmax(f$p_hat[f$y_mat == 1L], .Machine$double.eps))
    1.0 - sum(lp) / (f$N * log(1.0 / f$J))
  }

  # ---- pre-format support v strings ----------------------------------------
  v_strs <- vapply(fits, function(f)
    paste0("[", round(min(f$v), 3L), ", ", round(max(f$v), 3L), "]"),
    character(1L))

  # ---- column widths --------------------------------------------------------
  #
  # Data column must accommodate:
  #   estimate + stars : e.g. "-0.1234***"  -> digits + 6 chars
  #   SE in parens     : e.g. "(0.1234)"    -> digits + 3 chars
  #   column labels and "(n)" index headers
  #   diagnostic values (support v can be ~16 chars)
  col_nums <- paste0("(", seq_len(n_m), ")")

  val_w <- max(
    digits + 6L,                   # "-X.XXXXn***" where n = optional digit
    max(nchar(column.labels)),
    max(nchar(col_nums)),
    max(nchar(v_strs))
  )

  # Left (label) column: accommodate coef names and all diagnostic labels
  stat_labels <- c("Log-likelihood", "McFadden R²", "S(p_hat)",
                   "S(w_hat)", "Misses", "Support v",
                   "Observations", "Parameters K")
  lbl_w <- max(
    if (length(all_coefs) > 0L) max(nchar(all_coefs)) else 0L,
    max(nchar(stat_labels))
  )

  total_w <- lbl_w + 2L + n_m * val_w + (n_m - 1L) * 2L
  sep_h   <- strrep("=", total_w)
  sep_l   <- strrep("-", total_w)

  # ---- string-alignment helpers --------------------------------------------
  # Left-justify s in a field of width w
  .lj <- function(s, w) formatC(as.character(s), width = -w, flag = "-")
  # Right-justify s in a field of width w
  .rj <- function(s, w) formatC(as.character(s), width =  w)
  # Centre s in a field of width w
  .cj <- function(s, w) {
    s <- as.character(s)
    n <- nchar(s)
    p <- max(0L, w - n)
    paste0(strrep(" ", p %/% 2L), s, strrep(" ", (p + 1L) %/% 2L))
  }

  # Build one full-width table row:
  #   label (left-justified, lbl_w) + 2 spaces + cells (val_w each, 2 spaces apart)
  .row <- function(label, cells, cell_fn = .rj) {
    cell_strs <- vapply(cells, cell_fn, character(1L), w = val_w)
    paste0(.lj(label, lbl_w), "  ", paste(cell_strs, collapse = "  "))
  }

  # ---- number-formatting helpers -------------------------------------------
  # Estimate with significance stars
  .fmte <- function(b, se) {
    pv <- 2.0 * pnorm(-abs(b / se))
    paste0(formatC(round(b, digits), format = "f", digits = digits), .st(pv))
  }
  # Standard error in parentheses
  .fmtse <- function(se)
    paste0("(", formatC(round(se, digits), format = "f", digits = digits), ")")
  # Fixed-decimal scalar
  .ff <- function(x, d)
    formatC(round(as.numeric(x), d), format = "f", digits = d)

  # ---- assemble table lines ------------------------------------------------
  lines <- character(0L)

  if (!is.null(title))
    lines <- c(lines, .cj(title, total_w))

  lines <- c(lines, sep_h)
  lines <- c(lines, .row("", col_nums,      cell_fn = .cj))
  lines <- c(lines, .row("", column.labels, cell_fn = .cj))
  lines <- c(lines, sep_l)
  lines <- c(lines, .lj(paste0("Dep. var.: ", dep.var.label), total_w))
  lines <- c(lines, sep_l)

  # ---- coefficient block ---------------------------------------------------
  for (nm in all_coefs) {
    est_cells <- character(n_m)
    se_cells  <- character(n_m)
    for (i in seq_len(n_m)) {
      b  <- fits[[i]]$beta[nm]
      se <- fits[[i]]$se_beta[nm]
      if (is.na(b)) {
        est_cells[i] <- ""
        se_cells[i]  <- ""
      } else {
        est_cells[i] <- .fmte(b, se)
        se_cells[i]  <- .fmtse(se)
      }
    }
    lines <- c(lines, .row(nm, est_cells))          # estimate row
    lines <- c(lines, .row("",  se_cells))          # SE row (blank label)
  }

  lines <- c(lines, sep_l)

  # ---- diagnostics block ---------------------------------------------------
  lines <- c(lines,
    .row("Log-likelihood",
         vapply(fits, function(f) .ff(f$loglik, 2L), character(1L)),
         cell_fn = .cj))

  lines <- c(lines,
    .row("McFadden R²",
         vapply(fits, function(f) .ff(.mcf(f), 4L), character(1L)),
         cell_fn = .cj))

  lines <- c(lines,
    .row("S(p_hat)",
         vapply(fits, function(f) .ff(f$S_p, 4L), character(1L)),
         cell_fn = .cj))

  lines <- c(lines,
    .row("S(w_hat)",
         vapply(fits, function(f) .ff(f$S_w, 4L), character(1L)),
         cell_fn = .cj))

  lines <- c(lines,
    .row("Misses",
         vapply(fits, function(f) as.character(f$misses), character(1L)),
         cell_fn = .cj))

  lines <- c(lines, .row("Support v", v_strs, cell_fn = .cj))

  lines <- c(lines,
    .row("Observations",
         vapply(fits, function(f) as.character(f$N), character(1L)),
         cell_fn = .cj))

  lines <- c(lines,
    .row("Parameters K",
         vapply(fits, function(f) as.character(f$K), character(1L)),
         cell_fn = .cj))

  lines <- c(lines, sep_h)

  if (stars)
    lines <- c(lines, "Note: *p<0.10; **p<0.05; ***p<0.01")

  # ---- print and optionally write to file ----------------------------------
  cat(paste(lines, collapse = "\n"), "\n")
  if (!is.null(out)) writeLines(lines, con = out)

  invisible(lines)
}

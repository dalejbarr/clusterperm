#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' Get the sign of a difference
#'
#' @param x A vector of length 2, a 2x2 matrix, or a 2x2x2 matrix.
#' @return An integer (1 or -1), indicating effect direction.
#' @details A convenience function that takes the output of \code{\link{model.tables}} applied to an \code{\link{aov}} object. Not exported.
.eff_sign <- function(x) {
  ## get the sign (direction) of the effect
  ## x: a vector of length 2, a 2x2 matrix, or a 2x2x2 matrix
  ##    if matrix provided, calculate interaction
  if (length(dim(x)) == 1L) { # it's a vector
    diff <- x[[1]] - x[[2]]
  } else if (identical(dim(x), c(2L, 2L))) {
      ## calculate two way interaction (a - b) - (c - d)
    diff <- (x[1, 1] - x[1, 2]) - (x[2, 1] - x[2, 2])
  } else {
    if (!identical(dim(x), c(2L, 2L, 2L)))
      stop("'x' must be a vector, 2x2 matrix, or 2x2x2 matrix")
    ## calculate three way interaction
    ## ((a - b) - (c - d)) - ((e - f) - (g - h))
    diff1 <- (x[1, 1, 1] - x[1, 1, 2]) - (x[1, 2, 1] - x[1, 2, 2])
    diff2 <- (x[2, 1, 1] - x[2, 1, 2]) - (x[2, 2, 1] - x[2, 2, 2])
    diff <- diff1 - diff2
  }
  return(sign(diff))
}

#' Run an analysis of variance
#'
#' @param x Data frame (a single bin).
#' @param formula Model formula (passed to \code{\link{aov}}). For within-subjects data, and \code{Error} term should be simple (e.g., \code{Error(id)}.
#' @return A data frame with effect, signed F statistic (indicating direction of the effect), and p value.
.tidy_anova <- function(x, formula) {
  ## run an RM-ANOVA on the data
  mm <- stats::aov(formula, x)
  ## store results in table
  result_tbl <- broom::tidy(mm) %>%
    dplyr::filter(term != "Residuals")

  ## remove Error() from formula and re-fit, otherwise
  ## model.tables will crash
  formstr <- deparse(formula)
  form2 <- stats::as.formula(sub("\\+?\\s*Error\\(.*\\)", "", formstr))
  
  ## get the marginal/cell means
  effs <- stats::model.tables(stats::aov(form2, x), type = "means")[["tables"]]

  ## figure out the signs of each effect
  signs <- sapply(effs[-1], .eff_sign)

  ## combine F statistic with effect sign
  result <- dplyr::inner_join(result_tbl,
             tibble::tibble(term = names(signs), sign = signs), "term")

  tibble::tibble(effect = dplyr::pull(result, term),
                 stat = dplyr::pull(result, statistic) *
                   dplyr::pull(result, sign),
                 p = dplyr::pull(result, p.value))
}

#' Bin-by-bin analysis of variance
#'
#' Run \code{\link{aov}} on each of a series of time bins.
#'
#' @param .data Data frame.
#' @param bin Unquoted name of variable identifying bins.
#' @param formula Model formula, passed to \code{\link{aov}}.
#' @return A data frame containing signed F statistics and p values for each time bin.
#' TODO examples
#' @export
aov_by_bin <- function(.data, bin, formula) {
  tidyr::nest(.data, d = c(-{{bin}})) %>%
    dplyr::mutate(.result = purrr::map(d, .tidy_anova, formula)) %>%
      dplyr::select(-d) %>%
      tidyr::unnest(c(.result))
}

#' Detect clusters and calculate mass statistics
#'
#' @param .data Data frame.
#' @param bin Unquoted name of variable identifying bins.
#' @param stat Unquoted name of variable with signed statistics.
#' @param p Unquoted name of variable with p values.
#' @param alpha Alpha level for each test.
#' @return A data frame listing onset (\code{b0}) and offset (\code{b0}) of each cluster, along with the sign of the effect and cluster mass statistic (\code{cms}).
#' @seealso \code{\link{detect_clusters_by_effect}} for more than one effect.
#' @export
detect_clusters <- function(.data, bin, stat, p, alpha = .05) {
  x2 <- .data %>%
    dplyr::arrange({{bin}})
  sig <- as.integer( (dplyr::pull(x2, {{p}}) < alpha) *
                     sign(dplyr::pull(x2, {{stat}})) )
  runs <- rle(sig)
  run_ix <- which(runs$values != 0L)
  mystat <- dplyr::pull(x2, {{stat}})
  if (length(run_ix) == 0L) {
    tibble::tibble(b0 = NA, b1 = NA, sign = NA_integer_, cms = 0)
  } else {
    bins <- dplyr::pull(x2, {{bin}})
    purrr::map_df(run_ix, function(.x) {
      if (.x == 1L) {
        t0 <- 1L
      } else {
        t0 <- sum(runs$lengths[seq_len(.x - 1L)]) + 1L
      }
      t1 <- t0 + runs$lengths[.x] - 1L
      tibble::tibble(b0 = ifelse(is.na(t0), NA, bins[t0]),
                     b1 = ifelse(is.na(t1), NA, bins[t1]),
                     sign = sign(mystat[t0]),
                     cms = abs(sum(mystat[t0:t1])))
    })
  }
}

#' Detect clusters and calculate mass statistics for multiple effects
#' 
#' @inheritParams detect_clusters
#' @param effect Unquoted name of variable identifying each effect.
#' @details Performs \code{\link{detect_clusters}} for each individual effect in the data frame.
#' @return Data frame with effects, cluster onsets (\code{b0}) and offsets (\code{b1}) and cluster mass statistics (\code{cms}).
#' @export 
detect_clusters_by_effect <- function(.data, effect, bin, stat, p, alpha = .05) {
  .data %>%
    tidyr::nest(data = c(-{{effect}})) %>%
    dplyr::mutate(clust = purrr::map(data,
                                     detect_clusters,
                                     {{bin}}, {{stat}}, {{p}})) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(c(clust))
}

#' Generate null-hypothesis distributions for cluster mass statistics
#'
#' @param n Number of Monte Carlo runs
#' @param .data Data frame.
#' @param bin Unquoted name of variable identifying bins.
#' @param formula Model formula, passed to \code{\link{aov_by_bin}}.
#' @param fn Unquoted name of function to use for exchanging labels.
#' @param ... Arguments to be passed to \code{fn}.
#' @return A data frame with null-hypothesis distributions for each effect in the column \code{nhd}.
#' @details Generates null-hypothesis distributions by the following algorithm: For each of the \code{n} Monte Carlo runs, (1) randomly permute labels in the dataset according to the relabeling function \code{fn}; (2) run \code{\link{aov_by_bin}} on the resulting data, the result of which is passed to (3) \code{\link{detect_clusters_by_effect}}, and (4) store the maximum cluster mass statistic for each effect on each run.
#' @export 
cluster_nhds <- function(n, .data, bin, formula, fn, ...) {
  res <- purrr::rerun(n, {
    fn(.data = .data, ...) %>%
      aov_by_bin({{bin}}, {{formula}}) %>%
      detect_clusters_by_effect(effect, {{bin}}, stat, p) %>%
      dplyr::group_by(effect) %>%
      dplyr::summarize(maxcms = max(cms)) %>%
      dplyr::ungroup()
  })

  res %>%
    dplyr::bind_rows() %>%
    tidyr::nest(data = c(-effect)) %>%
    dplyr::mutate(nhd = purrr::map(data, dplyr::pull, maxcms)) %>%
    dplyr::select(-data)
}

#' Calculate permutation p-value from null-hypothesis distribution
#'
#' @param cms Cluster mass statistic for the original data.
#' @param nhd A vector comprising the null-hypothesis distribution, excluding the original statistic.
#' @return A scalar representing the p-value.
#' @details Calculates the p-value as the number of values in the vector \code{c(cms, nhd)} that are greater than or equal to the original value, divided by the length of \code{nhd} plus one (because the original statistic is included in the NHD).
#' @export
pvalue <- function(cms, nhd) {
  sum(c(cms, nhd) >= cms) / (length(nhd) + 1L)
}

#' Calculate permutation p-values from null-hypothesis distributions
#' @param orig Data frame containing clusters for multiple effects detected on original data.
#' @param nhds Data frame with nhds for each effect.
#' @param key Key on which to join the data frames (i.e., quoted name of the variable identifying effects).
#' @return Original data frame with extra column containing p-values for each effect and cluster.
#' @export
pvalues <- function(orig, nhds, key = "effect") {
  ## TODO: make sure there are nhds for everything
  dplyr::left_join(orig, nhds, key) %>%
    dplyr::mutate(p = purrr::map2_dbl(cms, nhd, pvalue)) %>%
      dplyr::select(-nhd)
}

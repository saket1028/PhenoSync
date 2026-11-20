#' Poole and Rathcke (1979) Peak-Spacing Index
#'
#' Extracts species peak flowering dates from a wide matrix and evaluates
#' whether peak dates are regularly spaced, aggregated, or approximately random
#' across the flowering season.
#'
#' Poole and Rathcke (1979) treat time as continuous: each peak date is placed
#' on an axis running from the start to the end of the growing season, and the
#' gaps between neighbouring peaks (plus the gaps from the season start to the
#' first peak and from the last peak to the season end) are measured as a
#' fraction of the season. The index therefore needs real census dates. If
#' `census_days` is not supplied, censuses are assumed to be one day apart with
#' no gaps, which will give misleading results for weekly or uneven censuses.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are
#'   censuses, columns are species, and values are counts of open flowers.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param census_days Numeric or `Date` vector giving the date of each census
#'   (one value per row of `data`, in increasing order), for example day of
#'   year. If `NULL`, censuses are assumed to be on consecutive days
#'   (1, 2, 3, ...).
#' @param season_start,season_end Start and end of the growing season, on the
#'   same scale as `census_days`. Default to the first and last census.
#' @param season_length Kept for compatibility with earlier versions. If
#'   supplied, `season_end` is set to `season_start + season_length`. Cannot be
#'   used together with `season_end`.
#' @param perform_test Logical. If `TRUE`, calculates the approximate
#'   chi-square test given by Poole and Rathcke (1979).
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with peak dates (on the `census_days` scale), ordered peaks,
#'   gaps between peaks, observed and expected spacing variance, ratio,
#'   optional test statistics, the season used, and an interpretation.
#'
#' @references Poole, R. W., & Rathcke, B. J. (1979). Regularity, randomness,
#'   and aggregation in flowering phenologies. Science, 203, 470-471.
#'
#' @examples
#' data(example_wide)
#' poole_rathcke(example_wide, verbose = FALSE)$ratio
#'
#' # Weekly censuses recorded as day of year
#' m <- matrix(0, nrow = 6, ncol = 3)
#' m[cbind(c(2, 3, 5), 1:3)] <- 5
#' poole_rathcke(m, census_days = c(120, 127, 134, 141, 148, 155),
#'               verbose = FALSE)$ratio
#'
#' @export
poole_rathcke <- function(data,
                          census_col = NULL,
                          census_days = NULL,
                          season_start = NULL,
                          season_end = NULL,
                          season_length = NULL,
                          perform_test = TRUE,
                          verbose = TRUE) {

  # STEP 0: Input checks

  flower_mat <- prepare_census_matrix(data, census_col = census_col, data_arg = "data")

  if (!is.numeric(flower_mat)) {
    stop("Input data must be a numeric matrix or data frame in wide format.")
  }

  if (any(is.na(flower_mat))) {
    stop("Input contains NA values. Please remove or replace missing values before running the function.")
  }

  if (any(flower_mat < 0)) {
    stop("Negative values are not allowed.")
  }

  k <- ncol(flower_mat)
  n_census <- nrow(flower_mat)

  if (k < 2) {
    stop("At least two species are required.")
  }

  if (n_census < 1) {
    stop("The dataset must contain at least one census.")
  }

  # STEP 0b: Census dates and season window

  if (is.null(census_days)) {
    census_days <- seq_len(n_census)
    if (verbose) {
      message("census_days not supplied: assuming one census per day with no gaps.")
    }
  }

  if (inherits(census_days, "Date")) {
    census_days <- as.numeric(census_days)
    if (inherits(season_start, "Date")) season_start <- as.numeric(season_start)
    if (inherits(season_end, "Date")) season_end <- as.numeric(season_end)
  }

  if (!is.numeric(census_days) || length(census_days) != n_census || any(is.na(census_days))) {
    stop("census_days must be a numeric or Date vector with one value per row of data and no missing values.")
  }

  if (n_census > 1 && any(diff(census_days) <= 0)) {
    stop("census_days must be in increasing order with no repeated dates.")
  }

  if (is.null(season_start)) {
    season_start <- census_days[1]
  }

  if (!is.null(season_length)) {
    if (!is.null(season_end)) {
      stop("Supply season_end or season_length, not both.")
    }
    if (!is.numeric(season_length) || length(season_length) != 1 || is.na(season_length) || season_length <= 0) {
      stop("season_length must be a single positive numeric value.")
    }
    season_end <- season_start + season_length
  }

  if (is.null(season_end)) {
    season_end <- census_days[n_census]
  }

  for (nm in c("season_start", "season_end")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1 || is.na(v)) {
      stop(nm, " must be a single numeric value (or Date).")
    }
  }

  if (season_end <= season_start) {
    stop("season_end must be later than season_start.")
  }

  ids <- colnames(flower_mat)
  if (is.null(ids)) {
    ids <- paste0("sp", seq_len(k))
    colnames(flower_mat) <- ids
  }

  # STEP 1: Extract peak dates from the wide flowering matrix
  # Peak date = date of the first census with the maximum flowering count

  peak_census <- rep(NA_integer_, k)
  names(peak_census) <- ids

  for (i in seq_len(k)) {
    x <- flower_mat[, i]
    if (!all(x == 0)) {
      peak_census[i] <- which.max(x)
    }
  }

  # Species with no flowering cannot have a peak date
  valid <- !is.na(peak_census)

  if (sum(valid) < 2) {
    stop("Fewer than two species have valid peak dates.")
  }

  if (any(!valid)) {
    warning(
      "Some species had no flowering and were excluded: ",
      paste(ids[!valid], collapse = ", ")
    )
  }

  peak_census <- peak_census[valid]
  peak_dates <- census_days[peak_census]
  names(peak_dates) <- ids[valid]
  ids_valid <- ids[valid]
  k_valid <- length(peak_dates)

  # Peaks must fall inside the season, otherwise gaps become negative
  outside <- peak_dates < season_start | peak_dates > season_end
  if (any(outside)) {
    stop(
      "Peak dates fall outside the season (", season_start, " to ", season_end,
      ") for: ", paste(ids_valid[outside], collapse = ", "),
      ". Widen season_start/season_end or check census_days."
    )
  }

  # STEP 2: Express peak dates as a fraction of the season and order them

  peak_df <- data.frame(
    species = ids_valid,
    peak_census = peak_census,
    peak_date = peak_dates,
    normalized_peak = (peak_dates - season_start) / (season_end - season_start),
    stringsAsFactors = FALSE
  )

  peak_df <- peak_df[order(peak_df$normalized_peak), , drop = FALSE]
  rownames(peak_df) <- NULL

  y <- peak_df$normalized_peak

  # STEP 3: Gaps between neighbouring peaks, including the season edges

  D <- diff(c(0, y, 1))

  interval_labels <- c(
    paste0("start_to_", peak_df$species[1]),
    if (k_valid > 1) paste0(peak_df$species[1:(k_valid - 1)], "_to_", peak_df$species[2:k_valid]),
    paste0(peak_df$species[k_valid], "_to_end")
  )

  distance_table <- data.frame(
    interval = interval_labels,
    distance = D,
    distance_days = D * (season_end - season_start),
    stringsAsFactors = FALSE
  )

  # STEP 4: Observed variance of the gaps (Poole & Rathcke's P)
  # The gaps always add up to 1, so their mean is 1 / (k + 1)

  D_bar <- mean(D)
  P <- sum((D - D_bar)^2) / (k_valid + 1)

  # STEP 5: Expected variance if peaks were placed at random (eqn 1)

  expected_P <- k_valid / ((k_valid + 1)^2 * (k_valid + 2))
  ratio <- P / expected_P

  # STEP 6: Optional chi-square test, chi2 = k P / E(P) with k degrees of freedom

  chi_square <- NA_real_
  df <- NA_integer_
  p_value <- NA_real_

  if (perform_test) {
    chi_square <- k_valid * ratio
    df <- k_valid
    p_value <- stats::pchisq(chi_square, df = df, lower.tail = FALSE)
  }

  # STEP 7: Interpretation

  tol <- 0.05

  interpretation <- if (ratio < (1 - tol)) {
    "Regular spacing"
  } else if (ratio > (1 + tol)) {
    "Aggregated peaks"
  } else {
    "Approximately random spacing"
  }

  if (verbose) {
    message("Poole & Rathcke index calculated for ", k_valid, " species with valid peak dates.")
    message("Season: ", season_start, " to ", season_end)
    message("Observed P = ", round(P, 5))
    message("Expected E(P) = ", round(expected_P, 5))
    message("P / E(P) ratio = ", round(ratio, 3))
    message("Interpretation: ", interpretation)
  }

  return(list(
    peak_dates = peak_dates,
    ordered_peaks = peak_df,
    distances = D,
    distance_table = distance_table,
    mean_distance = D_bar,
    P = P,
    expected_P = expected_P,
    ratio = ratio,
    chi_square = chi_square,
    df = df,
    p_value = p_value,
    season = c(start = season_start, end = season_end),
    interpretation = interpretation
  ))
}

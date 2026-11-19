#' Mahoro (2002) Flowering Rank Order
#'
#' Calculates Mahoro's flowering rank order (Ri; eqn 1 of Mahoro 2002) using
#' cumulative flowering proportions and ranks across census days up to the
#' population peak. Ri orders plants by earliness: a smaller Ri means earlier
#' flowering. It is a ranking of plants, not a measure of synchrony, so it has
#' no single population value.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are individuals, and values are counts of open flowers.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param ties.method Character. Method passed to [base::rank()] for tied values.
#' @param verbose Logical. If `TRUE`, prints a short progress message.
#'
#' @return A list with `Ri`, `peak_day`, `daily_ranks`, and
#'   `cumulative_proportion`.
#'
#' @examples
#' data(example_wide)
#' mahoro_Ri(example_wide, verbose = FALSE)$Ri
#'
#' @export
mahoro_Ri <- function(data, census_col = NULL, ties.method = "average", verbose = TRUE) {
  
  flower_mat <- prepare_census_matrix(data, census_col = census_col, data_arg = "data")
  
# Input check  
  if (!is.numeric(flower_mat)) {
    stop("Input data must be numeric.")
  }
  
  n <- ncol(flower_mat)
  ids <- colnames(flower_mat)
  
  if (is.null(ids)) {
    ids <- paste0("ind", seq_len(n))
    colnames(flower_mat) <- ids
  }
  
  if (n < 2) {
    stop("At least two individuals are required.")
  }
  
  if (any(flower_mat < 0, na.rm = TRUE)) {
    stop("Negative values are not allowed.")
  }
  
  total_flowers <- colSums(flower_mat, na.rm = TRUE)
  
  if (all(total_flowers == 0)) {
    stop("All individuals have total seasonal flowers = 0.")
  }
  
  if (any(total_flowers == 0)) {
    warning(
      "Some individuals have total seasonal flowers = 0. Their Ri may be difficult to interpret: ",
      paste(ids[total_flowers == 0], collapse = ", ")
    )
  }
  
  # cumulative flower production
  cum_mat <- apply(flower_mat, 2, cumsum)
  
  # if only one column, apply() can simplify awkwardly
  if (is.null(dim(cum_mat))) {
    cum_mat <- matrix(cum_mat, ncol = 1)
    colnames(cum_mat) <- ids
  }
  
  # cumulative proportion for each individual
  x_mat <- matrix(NA_real_, nrow = nrow(flower_mat), ncol = ncol(flower_mat))
  rownames(x_mat) <- rownames(flower_mat)
  colnames(x_mat) <- ids
  
  valid <- total_flowers > 0
  x_mat[, valid] <- sweep(cum_mat[, valid, drop = FALSE], 2, total_flowers[valid], FUN = "/")
  
  # population peak day = census day with maximum total newly opened flowers
  pop_daily_total <- rowSums(flower_mat, na.rm = TRUE)
  peak_day <- which.max(pop_daily_total)
  
  if (length(peak_day) == 0 || is.na(peak_day) || peak_day < 1) {
    stop("Could not determine population peak day.")
  }
  
  Ri_sum <- numeric(n)
  names(Ri_sum) <- ids
  Ri_sum[] <- 0
  
  daily_ranks <- matrix(NA_real_, nrow = peak_day, ncol = n)
  colnames(daily_ranks) <- ids
  rownames(daily_ranks) <- rownames(flower_mat)[1:peak_day]
  
  for (j in 1:peak_day) {
    # Rank 1 = largest share of flowers already opened, i.e. earliest plant
    # (Mahoro 2002: "A smaller rank number indicates earlier blooming")
    ranks <- rank(-x_mat[j, ], ties.method = ties.method, na.last = "keep")
    Ri_sum <- Ri_sum + ranks
    daily_ranks[j, ] <- ranks
  }
  
  # Mahoro (2002) eqn 1: Ri is the sum of daily ranks up to the peak
  Ri <- Ri_sum
  
  if (verbose) {
    message("Mahoro Ri calculated for ", n, " individuals.")
    message("Peak day identified at census day ", peak_day, ".")
  }
  
  return(list(
    Ri = Ri,
    peak_day = peak_day,
    daily_ranks = daily_ranks,
    cumulative_proportion = x_mat
  ))
}


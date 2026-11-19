#' Freitas and Bolmgren (2008) Population Synchrony Index
#'
#' Calculates population-level flowering synchrony from flowering-intensity
#' curves after optional standardization of each individual by its maximum value.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are individuals, and values are flowering intensities.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param standardize Character. Use `"max"` to divide each individual schedule
#'   by its own maximum, or `"none"` if values are already relative intensities.
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with `relative_intensity`, `Ti`, `di`, and `population_S`.
#'
#' @details
#' This function follows the equations for \eqn{d_i} and \eqn{S} in Freitas &
#' Bolmgren (2008, p. 722) exactly. The paper states that \eqn{d_i} and
#' \eqn{S} range from 0 (no synchrony) to 1 (perfect synchrony). With these
#' equations, however, plants with identical flowering schedules only score 1
#' if every census with flowers is at full intensity (a flat schedule) or if
#' flowering lasts a single census. When intensity rises and falls, identical
#' schedules score less than 1: for example, about 0.60 for counts of 1, 2, 3,
#' 2, 1 and about 0.28 for 1, 1, 10, 1, 1. Values should therefore be compared
#' only between populations or species with similarly shaped schedules, and
#' should not be read as "1 = perfect synchrony" unless schedules are flat.
#'
#' @examples
#' data(example_wide)
#' freitas_bolmgren_S(example_wide, verbose = FALSE)$population_S
#'
#' @export
freitas_bolmgren_S <- function(data,
                               census_col = NULL,
                               standardize = c("max", "none"),
                               verbose = TRUE) {
  
  # STEP 0: Input checks
  
  standardize <- match.arg(standardize)
  
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
  
  N <- ncol(flower_mat)
  T_total <- nrow(flower_mat)
  ids <- colnames(flower_mat)
  
  if (is.null(ids)) {
    ids <- paste0("ind", seq_len(N))
    colnames(flower_mat) <- ids
  }
  
  if (N < 2) {
    stop("At least two individuals are required.")
  }
  
  # STEP 1: Convert to relative intensity f_i,t
  
  f_mat <- matrix(NA_real_, nrow = nrow(flower_mat), ncol = ncol(flower_mat))
  rownames(f_mat) <- rownames(flower_mat)
  colnames(f_mat) <- ids
  
  if (standardize == "max") {
    max_vals <- apply(flower_mat, 2, max, na.rm = TRUE)
    
    if (any(max_vals == 0)) {
      warning(
        "Some individuals have all-zero flowering schedules. Their d_i will be NA: ",
        paste(ids[max_vals == 0], collapse = ", ")
      )
    }
    
    valid <- max_vals > 0
    f_mat[, valid] <- sweep(flower_mat[, valid, drop = FALSE], 2, max_vals[valid], FUN = "/")
    
  } else if (standardize == "none") {
    if (any(flower_mat > 1)) {
      warning("standardize = 'none' was selected, but some values are > 1. The index assumes relative intensities from 0 to 1.")
    }
    f_mat <- flower_mat
  }
 
  # STEP 2: Calculate T_i
  # T_i = number of census days where f_i,t > 0

  Ti <- colSums(f_mat > 0, na.rm = TRUE)
  

  # STEP 3: Calculate d_i for each individual
 
  di <- numeric(N)
  names(di) <- ids
  di[] <- NA_real_
  
  for (i in seq_len(N)) {
    
    if (Ti[i] == 0) {
      next
    }
    
    sum_ijt <- 0
    
    for (j in seq_len(N)) {
      if (j != i) {
        sum_ijt <- sum_ijt + sum(sqrt(f_mat[, i] * f_mat[, j]), na.rm = TRUE)
      }
    }
    
    di[i] <- (1 / Ti[i]) * (1 / (N - 1)) * sum_ijt
  }
  
  
  # STEP 4: Population synchrony
 
  S <- mean(di, na.rm = TRUE)
  if (is.nan(S)) S <- NA_real_
  
 

  if (verbose) {
    message("Freitas & Bolmgren synchrony calculated for ", N, " individuals.")
    message("Population synchrony S = ", round(S, 4))
  }
  
  return(list(
    relative_intensity = f_mat,
    Ti = Ti,
    di = di,
    population_S = S
  ))
}


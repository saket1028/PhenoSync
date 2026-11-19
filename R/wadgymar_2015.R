#' Wadgymar et al. (2015) Individual Synchrony Index
#'
#' Calculates an individual-level index describing concentration or spread of
#' flowering intensity across the individual's flowering days.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are individuals, and values are counts of open flowers.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with `individual_results`, `Si_individual`, and `mean_Si`.
#'
#' @examples
#' data(example_wide)
#' wadgymar_Si(example_wide, verbose = FALSE)$mean_Si
#'
#' @export
wadgymar_Si <- function(data, census_col = NULL, verbose = TRUE) {
  
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
  
  n <- ncol(flower_mat)
  ids <- colnames(flower_mat)
  
  if (is.null(ids)) {
    ids <- paste0("ind", seq_len(n))
    colnames(flower_mat) <- ids
  }
  
  if (n < 1) {
    stop("At least one individual is required.")
  }
  
  
  # STEP 1: Calculate Si for each individual
  Si <- numeric(n)
  names(Si) <- ids
  
  Ti <- numeric(n)   # total flowers counted
  Di <- numeric(n)   # schedule duration
  CVi <- numeric(n)  # coefficient of variation
  SSi <- numeric(n)  # sum of squared deviations
  
  for (i in seq_len(n)) {
    
    x <- flower_mat[, i]
    
    # Keep only flowering days for this individual's display schedule
    x_pos <- x[x > 0]
    
    Ti[i] <- sum(x_pos)
    Di[i] <- length(x_pos)
    
    # No flowers at all -> undefined
    if (Di[i] == 0 || Ti[i] == 0) {
      Si[i] <- NA_real_
      CVi[i] <- NA_real_
      SSi[i] <- NA_real_
      next
    }
    
    # If flowering occurred on only one day, assign maximum synchrony
    if (Di[i] == 1) {
      Si[i] <- 1
      CVi[i] <- NA_real_
      SSi[i] <- 0
      next
    }
    
    # Mean daily flower production across flowering days only
    mu <- Ti[i] / Di[i]
    
    # Sum of squared deviations from even daily distribution
    SSi[i] <- sum((x_pos - mu)^2)
    
    # Sample coefficient of variation over flowering days only
    x_sd <- stats::sd(x_pos)
    x_mean <- mean(x_pos)
    CVi[i] <- ifelse(x_mean == 0, NA_real_, x_sd / x_mean)
    
    # Paper's rearranged equation (Wadgymar et al. 2015, eqn 2):
    # Si = sqrt(SSi) / Ti * sqrt(Di / (Di - 1)), equivalent to CVi / sqrt(Di)
    Si[i] <- sqrt(SSi[i]) / Ti[i] * sqrt(Di[i] / (Di[i] - 1))
  }
  
  mean_Si <- mean(Si, na.rm = TRUE)
  if (is.nan(mean_Si)) mean_Si <- NA_real_
  
  results_df <- data.frame(
    individual_id = ids,
    Ti = Ti,
    Di = Di,
    CVi = CVi,
    SSi = SSi,
    Si = Si,
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    message("Wadgymar individual synchrony calculated for ", n, " individuals.")
    message("Mean Si = ", round(mean_Si, 4))
  }
  
  #Step 3: Compile results
  return(list(
    individual_results = results_df,
    Si_individual = Si,
    mean_Si = mean_Si
  ))
}


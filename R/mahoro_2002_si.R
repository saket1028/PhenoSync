#' Mahoro (2002) Flowering Synchrony Index
#'
#' Calculates Mahoro's individual synchrony index from a wide flowering matrix,
#' following equation 2 of Mahoro (2002):
#' \deqn{S_i = \frac{1}{2}\left(2 - \sum_{j=1}^{n} |y_{i,j} - \bar{y}_{i,j}|\right)}
#' where \eqn{y_{i,j}} is the share of individual \eqn{i}'s season total that
#' was counted on census \eqn{j}, and \eqn{\bar{y}_{i,j}} is the mean of the
#' same share across all other individuals. Si is 1 when an individual's
#' schedule matches the average of the others exactly, and 0 when the two
#' never overlap.
#'
#' @references Mahoro, S. (2002) Individual flowering schedule, fruit set, and
#'   flower and seed predation in \emph{Vaccinium hirtum} Thunb. (Ericaceae).
#'   \emph{Canadian Journal of Botany} 80: 82--92.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are individuals, and values are counts of open flowers.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param verbose Logical. If `TRUE`, prints a short progress message.
#'
#' @return A named numeric vector of individual synchrony values.
#'
#' @examples
#' data(example_wide)
#' mahoro_Si(example_wide, verbose = FALSE)
#'
#' @export
mahoro_Si <- function(data, census_col = NULL, verbose = TRUE) {
  
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
  
  if (any(total_flowers == 0)) {
    warning(
      "Some individuals have total seasonal flowers = 0. Their Si will be NA: ",
      paste(ids[total_flowers == 0], collapse = ", ")
    )
  }
  
  # y_i,j = proportion of total seasonal flowers of individual i
  y_mat <- matrix(NA_real_, nrow = nrow(flower_mat), ncol = ncol(flower_mat))
  rownames(y_mat) <- rownames(flower_mat)
  colnames(y_mat) <- ids
  
  valid <- total_flowers > 0
  y_mat[, valid] <- sweep(flower_mat[, valid, drop = FALSE], 2, total_flowers[valid], FUN = "/")
  
  Si <- numeric(n)
  names(Si) <- ids
  Si[] <- NA_real_
  
  for (i in 1:n) {
    
    if (total_flowers[i] == 0) {
      next
    }
    
    others <- setdiff(1:n, i)
    
    # if all "other" individuals have zero flowers, synchrony cannot be computed
    valid_others <- others[total_flowers[others] > 0]
    
    if (length(valid_others) == 0) {
      Si[i] <- NA_real_
      next
    }
    
    # mean proportional schedule of all other flowering individuals
    y_bar_minus_i <- rowMeans(y_mat[, valid_others, drop = FALSE], na.rm = TRUE)
    
    # Mahoro (2002) eqn 2: Si = (1/2) * (2 - sum_j |y_ij - y_bar_ij|)
    Si[i] <- 0.5 * (2 - sum(abs(y_mat[, i] - y_bar_minus_i), na.rm = TRUE))
  }
  
  if (verbose) {
    message("Mahoro Si calculated for ", n, " individuals.")
  }
  
  return(Si)
}


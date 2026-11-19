#' Calculate multiple phenological synchrony indices from one wide dataset
#'
#' Calculates several compatible phenological synchrony indices from the same
#' wide phenology matrix or data frame. This function is intentionally limited
#' to the main wide-data workflow used by most PhenoSync indices.
#'
#' @param data Phenology dataset. The function assumes that one row is one
#'   census day, one column is one individual or species as specified by
#'   `columns_are`, and each cell value is flowering state, flowering count,
#'   abundance, or intensity. If census names or dates are stored in a data
#'   column, specify that column using `census_col`.
#' @param indices Character vector of index function names to calculate, or
#'   `"auto"` to calculate all indices compatible with `columns_are` and
#'   `binary`.
#' @param columns_are Character. Specifies whether columns represent
#'   `"individuals"` or `"species"`. This controls which indices are selected
#'   when `indices = "auto"`.
#' @param binary Logical or `NULL`. Use `TRUE` when the data are
#'   presence/absence, `FALSE` when they contain abundance/intensity values, or
#'   `NULL` to detect this automatically from the matrix values.
#' @param census_col Optional column name or column index containing census
#'   names/dates. When supplied, this column is removed from the phenology
#'   matrix and used as row names before any index is calculated. For ordinary
#'   CSV/Excel imports, use `census_col = 1` or `census_col = "Date"`.
#' @param first_row_colnames Logical. If `TRUE`, raw wide tables where the first
#'   row contains column names are cleaned automatically before analysis. This is
#'   useful when data were imported with `header = FALSE`.
#' @param first_col_rownames Logical. Deprecated/backwards-compatible option.
#'   Prefer `census_col`; set this to `TRUE` only if the first column should be
#'   treated as census labels without naming it explicitly.
#' @param index_args Optional named list of index-specific argument lists. For
#'   example, `list(primack_Z = list(strong_threshold = 0.75))`.
#' @param make_plot Logical. If `TRUE`, draws a base R two-panel summary plot
#'   and stores a recorded base plot in the returned `plot` component when possible.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param ... Additional arguments passed to individual index functions when the
#'   argument names match their formal arguments.
#'
#' @return A list with data structure, column unit, binary status, possible
#'   indices, calculated indices, skipped/failed indices, full results, a
#'   one-row-per-index summary, plot data, an optional recorded base R plot,
#'   the prepared wide data matrix, census labels, and notes.
#'
#' @details
#' `calculate_indices()` does not handle generic long/tidy data and does not
#' calculate specialised index-specific structures such as Bolmgren flower-level
#' data, Dainou inter-sex data, Weis trait-schedule data, or Schoener
#' neighbourhood data. Use the corresponding single-index functions for those
#' specialised inputs.
#'
#' @examples
#' x <- make_pheno_dummy(format = "wide", n_species = 5, n_days = 30, seed = 1)
#' res <- calculate_indices(
#'   x,
#'   columns_are = "species",
#'   binary = FALSE,
#'   indices = c("albert_Si", "pianka_overlap"),
#'   verbose = FALSE,
#'   make_plot = FALSE
#' )
#' res$summary
#'
#' @export
calculate_indices <- function(data,
                              indices = "auto",
                              columns_are = c("individuals", "species"),
                              binary = NULL,
                              census_col = NULL,
                              first_row_colnames = TRUE,
                              first_col_rownames = FALSE,
                              index_args = list(),
                              make_plot = TRUE,
                              verbose = TRUE,
                              ...) {
  columns_are <- match.arg(columns_are)
  validate_binary_argument(binary)
  extra_args <- list(...)

  data_for_analysis <- .prepare_pheno_matrix(
    data = data,
    census_col = census_col,
    first_row_colnames = first_row_colnames,
    first_col_rownames = first_col_rownames
  )

  if (!is.matrix(data_for_analysis) && !is.data.frame(data_for_analysis)) {
    stop("data must be a wide matrix or data frame after preparation.", call. = FALSE)
  }

  binary_state <- detect_binary_values(data_for_analysis, binary)

  if (isTRUE(binary_state)) {
    data_for_analysis <- binary_matrix(data_for_analysis)
  }

  possible <- filter_indices_by_context(
    candidate_indices = wide_candidate_indices(),
    columns_are = columns_are,
    binary = binary_state
  )

  if (identical(indices, "auto")) {
    selected <- possible
  } else {
    selected <- indices
    not_possible <- setdiff(selected, possible)
    if (length(not_possible) > 0 && verbose) {
      warning(
        "The following requested indices were not automatically identified as suitable for this wide-data context: ",
        paste(not_possible, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (length(selected) == 0) {
    empty_summary <- index_summary_from_results(list())
    return(list(
      data_structure = "wide",
      detected_type = "wide",
      columns_are = columns_are,
      binary = binary_state,
      possible_indices = possible,
      calculated_indices = character(0),
      skipped = data.frame(index = character(0), reason = character(0)),
      results = list(),
      summary = empty_summary,
      plot_data = empty_summary,
      plot = NULL,
      prepared_data = data_for_analysis,
      census = rownames(data_for_analysis),
      notes = "Prepared wide phenology data. Rows are census days and columns are individuals or species."
    ))
  }

  run_one <- function(index_name) {
    if (!exists(index_name, mode = "function", inherits = TRUE)) {
      return(list(ok = FALSE, value = NULL, error = paste0("Function not found: ", index_name, ".")))
    }

    fn <- get(index_name, mode = "function", inherits = TRUE)
    fmls <- names(formals(fn))

    base_args <- switch(
      index_name,
      albert_Si = list(data = data_for_analysis),
      augspurger_Xi = list(flower_mat = data_for_analysis),
      freitas_bolmgren_S = list(data = data_for_analysis),
      koenig_rp = list(data = data_for_analysis),
      mahoro_Si = list(data = data_for_analysis),
      mahoro_Ri = list(data = data_for_analysis),
      marquis_SM = list(data = data_for_analysis),
      michalski_durka_CV = list(data = data_for_analysis),
      pianka_overlap = list(data = data_for_analysis),
      poole_rathcke = list(data = data_for_analysis),
      primack_Z = list(data = data_for_analysis),
      shrotri_SI_temp = list(data = data_for_analysis),
      shrotri_SI_abd = list(data = data_for_analysis),
      wadgymar_Si = list(data = data_for_analysis),
      list(data = data_for_analysis)
    )

    if ("verbose" %in% fmls && !("verbose" %in% names(base_args))) {
      base_args$verbose <- FALSE
    }

    idx_args <- index_args[[index_name]]
    if (is.null(idx_args)) idx_args <- list()
    all_args <- c(base_args, extra_args, idx_args)
    all_args <- all_args[names(all_args) %in% fmls]

    res <- tryCatch(
      do.call(fn, all_args),
      error = function(e) e
    )

    if (inherits(res, "error")) {
      return(list(ok = FALSE, value = NULL, error = conditionMessage(res)))
    }

    list(ok = TRUE, value = res, error = NA_character_)
  }

  results <- list()
  skipped <- data.frame(index = character(0), reason = character(0), stringsAsFactors = FALSE)

  for (idx in selected) {
    tmp <- run_one(idx)
    if (isTRUE(tmp$ok)) {
      results[[idx]] <- tmp$value
    } else {
      skipped <- rbind(skipped, data.frame(index = idx, reason = tmp$error, stringsAsFactors = FALSE))
    }
  }

  summary <- index_summary_from_results(
    results,
    columns_are = columns_are,
    binary = binary_state
  )

  index_plot <- NULL
  if (isTRUE(make_plot) && nrow(summary) > 0) {
    plot_index_summary(summary)
    index_plot <- tryCatch(
      grDevices::recordPlot(),
      error = function(e) NULL
    )
  }

  census_names <- rownames(data_for_analysis)

  if (verbose) {
    message("Data structure: wide")
    message("Rows represent census days; columns represent: ", columns_are)
    if (!is.null(binary_state) && !is.na(binary_state)) message("Binary data: ", binary_state)
    message("Possible indices: ", paste(possible, collapse = ", "))
    message("Successfully calculated: ", paste(names(results), collapse = ", "))
    if (nrow(skipped) > 0) {
      message("Skipped/failed: ", paste(skipped$index, collapse = ", "))
    }
  }

  list(
    data_structure = "wide",
    detected_type = "wide",
    columns_are = columns_are,
    binary = binary_state,
    possible_indices = possible,
    calculated_indices = names(results),
    skipped = skipped,
    results = results,
    summary = summary,
    plot_data = summary,
    plot = index_plot,
    prepared_data = data_for_analysis,
    census = census_names,
    notes = "Prepared wide phenology data. Rows are census days and columns are individuals or species."
  )
}

#' Summarise results from multiple indices
#'
#' @param results Named list of PhenoSync index outputs.
#' @param columns_are Character describing whether wide-data columns are
#'   individuals or species.
#' @param binary Logical describing whether data were treated as binary.
#'
#' @return Data frame with one summary estimate per index.
#'
#' @keywords internal
#' @noRd
index_summary_from_results <- function(results, columns_are = NA_character_, binary = NA) {
  if (length(results) == 0) {
    return(data.frame(
      index = character(0),
      estimate = numeric(0),
      range_class = character(0),
      range_panel = character(0),
      columns_are = character(0),
      binary = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  extract <- function(index_name, res) {
    val <- switch(
      index_name,
      albert_Si = res$population_mean,
      augspurger_Xi = res$Z,
      bolmgren_S = res$mean_S,
      dainou_SI = res$population_S,
      freitas_bolmgren_S = res$population_S,
      koenig_rp = res$population_rp,
      mahoro_Si = mean(res, na.rm = TRUE),
      marquis_SM = res$SM_population_mean,
      michalski_durka_CV = res$population_CV,
      pianka_overlap = res$population_mean,
      poole_rathcke = res$ratio,
      primack_Z = res$Z,
      schoener_SI = res$mean_SI,
      shrotri_SI_temp = res$community_mean,
      shrotri_SI_abd = res$community_mean,
      wadgymar_Si = res$mean_Si,
      weis_sync = res$S,
      NA_real_
    )
    as.numeric(val)[1]
  }

  # Indices that rank plants rather than measure synchrony have no single
  # population value: mahoro_Ri is a rank, so its average is always (n + 1) / 2.
  # Their full results are still returned, but they are left out of the summary.
  no_summary <- c("mahoro_Ri")
  summary_names <- setdiff(names(results), no_summary)
  if (length(summary_names) == 0) {
    return(index_summary_from_results(list()))
  }

  out <- data.frame(
    index = summary_names,
    estimate = vapply(summary_names, function(nm) extract(nm, results[[nm]]), numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  meta <- pheno_index_rules()
  out$range_class <- meta$range_class[match(out$index, meta$index)]
  out$range_class[is.na(out$range_class)] <- "other"
  out$range_panel <- ifelse(
    out$range_class == "bounded_0_1",
    "Indices bounded 0 to 1",
    "Indices on other scales"
  )
  out$columns_are <- columns_are
  out$binary <- binary
  out
}

# Internal dictionary of synchrony indices.
# This function is intentionally not exported; use define_index() instead.
pheno_index_database <- function() {
  list(
    albert_Si = list(
      name = "Albert (2001) pairwise flowering overlap index",
      function_name = "albert_Si",
      citation = "Albert, M. J., Escudero, A., & Iriondo, J. M. (2001). Female reproductive success of narrow endemic Erodium paularense in contrasting microhabitats. Ecology, 82(6), 1734-1747. https://doi.org/10.1890/0012-9658(2001)082[1734:FRSONE]2.0.CO;2",
      aspect = "Pairwise temporal overlap / timing similarity based on presence-absence.",
      formula = "O_ij = a_ij / b_ij; S_i = mean_j(O_ij), where a_ij is shared flowering days and b_ij is union flowering days.",
      input = "Wide matrix/data frame: rows = census days; columns = individuals/species; values = counts or intensities. Values > 0 are treated as presence.",
      output = "List: individual_Si and population_mean.",
      range = "Usually 0 to 1. Higher values indicate greater overlap."
    ),
    augspurger_Xi = list(
      name = "Augspurger (1983) individual and population synchrony",
      function_name = "augspurger_Xi",
      citation = "Augspurger, C. K. (1983). Phenology, flowering synchrony, and fruit set of six neotropical shrubs. Biotropica, 257-267. https://doi.org/10.2307/2387650",
      aspect = "Individual synchrony with the rest of the population based on flowering duration and overlap.",
      formula = "X_i = sum_j(e_ij) / ((n - 1) f_i); Z = mean_i(X_i).",
      input = "Wide matrix/data frame: rows = census days; columns = individuals; values = flowering counts. Values > 0 are treated as presence.",
      output = "List: Xi_table, Z, presence_matrix, overlap_matrix.",
      range = "Usually 0 to 1. Higher Z indicates stronger population synchrony."
    ),
    bolmgren_S = list(
      name = "Bolmgren (1998) within-individual synchrony",
      function_name = "bolmgren_S",
      citation = "Bolmgren, K. (1998). The use of synchronization measures in studies of plant reproductive phenology. Oikos, 82(2), 411-415. https://doi.org/10.2307/3546986",
      aspect = "Within-individual spread of flower opening relative to floral persistence.",
      formula = "S_i = mean floral persistence_i / SD(onset dates_i).",
      input = "Long flower-level data frame: one row per flower with plant_id, onset_day, and end_day.",
      output = "List: individual_results and mean_S.",
      range = "Unbounded positive scale; Inf possible when all onset dates are identical."
    ),
    dainou_SI = list(
      name = "Dainou et al. (2012) specifically designed to measure synchrony between males and females in sexually dimorphic species",
      function_name = "dainou_SI",
      citation = "Daïnou, K., Laurenty, E., Mahy, G., Hardy, O. J., Brostaux, Y., Tagg, N., & Doucet, J. L. (2012). Phenological patterns in a natural population of a tropical timber tree species, Milicia excelsa (Moraceae): Evidence of isolation by time and its interaction with feeding strategies of dispersers. American Journal of Botany, 99(9), 1453-1463. https://doi.org/10.3732/ajb.1200147",
      aspect = "Synchrony of male flowering with availability of female flowering.",
      formula = "d_j = mean_t(p_F,t) over days when male j flowers; S = mean_j(d_j).",
      input = "Data frame with id column, sex column, and numeric census-day columns. Values > 0 are treated as flowering presence.",
      output = "List: male_dj, female_prop_by_day, population_S.",
      range = "0 to 1. Higher values indicate stronger male-female synchrony."
    ),
    freitas_bolmgren_S = list(
      name = "Freitas & Bolmgren (2008) population synchrony",
      function_name = "freitas_bolmgren_S",
      citation = "Freitas, L., & Bolmgren, K. (2008). Synchrony is more than overlap: measuring phenological synchronization considering time length and intensity. Brazilian Journal of Botany, 31(4), 721-724. https://doi.org/10.1590/S0100-84042008000400017",
      aspect = "Distributional overlap among flowering-intensity curves.",
      formula = "d_i = (1/T_i)(1/(N-1)) sum_{j != i} sum_t sqrt(f_i,t f_j,t); S = mean_i(d_i).",
      input = "Wide matrix/data frame with flowering intensities. Can standardize each column by its maximum.",
      output = "List: relative_intensity, Ti, di, population_S.",
      range = "Usually 0 to 1 after max standardization."
    ),
    koenig_rp = list(
      name = "Koenig et al. (2003) mean pairwise correlation synchrony",
      function_name = "koenig_rp",
      citation = "Koenig, W. D., Kelly, D., Sork, V. L., Duncan, R. P., Elkinton, J. S., Peltonen, M. S., & Westfall, R. D. (2003). Dissecting components of population‐level variation in seed production and the evolution of masting behavior. Oikos, 102(3), 581–591. https://doi.org/10.1034/j.1600-0706.2003.12272.x",
      aspect = "Correlation in temporal abundance/intensity fluctuations.",
      formula = "r_p = mean of all pairwise Pearson correlations among flowering schedules.",
      input = "Wide numeric matrix/data frame; rows = census days; columns = individuals/species.",
      output = "List: pairwise_correlation, pairwise_values, individual_mean_r, population_rp.",
      range = "-1 to 1. Higher values indicate more correlated schedules."
    ),
    mahoro_Si = list(
      name = "Mahoro (2002) proportional-schedule synchrony",
      function_name = "mahoro_Si",
      citation = "Mahoro, S. (2002). Individual flowering schedule, fruit set, and flower and seed predation in Vaccinium hirtum Thunb. (Ericaceae). Canadian Journal of Botany, 80(1), 82–92. https://doi.org/10.1139/b01-136",
      aspect = "Similarity between each individual's proportional flowering schedule and the rest of the population.",
      formula = "S_i = (1/2) * (2 - sum_t |y_i,t - y_bar_-i,t|), where y_i,t is proportional seasonal flowering (Mahoro 2002, eqn 2).",
      input = "Wide numeric matrix/data frame; rows = census days; columns = individuals.",
      output = "Named numeric vector of individual Si values.",
      range = "0 to 1. Higher values indicate greater synchrony."
    ),
    mahoro_Ri = list(
      name = "Mahoro (2002) rank-order index",
      function_name = "mahoro_Ri",
      citation = "Mahoro, S. (2002). Individual flowering schedule, fruit set, and flower and seed predation in Vaccinium hirtum Thunb. (Ericaceae). Canadian Journal of Botany, 80(1), 82–92. https://doi.org/10.1139/b01-136",
      aspect = "Rank order of cumulative flowering up to the population peak day.",
      formula = "R_i is derived from daily ranks of cumulative proportional flowering up to the population peak day.",
      input = "Wide numeric matrix/data frame; rows = census days; columns = individuals.",
      output = "List: Ri, peak_day, daily_ranks, cumulative_proportion.",
      range = "Relative/rank-based scale depending on implementation."
    ),
    marquis_SM = list(
      name = "Marquis (1988) synchrony measure",
      function_name = "marquis_SM",
      citation = "Marquis, R. J. (1988). Phenological Variation in the Neotropical Understory Shrub Piper Arielanum: Causes and Consequences. Ecology, 69(5), 1552–1565. https://doi.org/10.2307/1941653",
      aspect = "Individual synchrony weighted by the number/proportion of other individuals flowering each day.",
      formula = "SM_i is calculated over days when individual i flowers using the population flowering proportion p_t.",
      input = "Wide numeric matrix/data frame; rows = census days; columns = individuals. Values > 0 are treated as presence.",
      output = "List: SM_individual, SM_population_mean, p_t, flowering_presence, total_output.",
      range = "Usually 0 to 1."
    ),
    michalski_durka_CV = list(
      name = "Michalski & Durka (2007) coefficient-of-variation index",
      function_name = "michalski_durka_CV",
      citation = "Michalski, S. G., & Durka, W. (2007). Synchronous Pulsed Flowering: Analysis of the Flowering Phenology in Juncus (Juncaceae). Annals of Botany, 100(6), 1271–1285. https://doi.org/10.1093/aob/mcm206",
      aspect = "Temporal spread/dispersion of flowering intensity.",
      formula = "CV = SD(x_t) / mean(x_t), calculated for population and individual flowering profiles.",
      input = "Wide numeric matrix/data frame; rows = census days; columns = individuals/species.",
      output = "List: population_profile, population_CV, individual_CV, mean_individual_CV.",
      range = "Non-negative; larger values indicate more uneven flowering over time."
    ),
    pianka_overlap = list(
      name = "Pianka (1974) niche-overlap index applied to phenology by Fleming & Partridge, 1984 (https://doi.org/10.1007/BF00384266)",
      function_name = "pianka_overlap",
      citation = "Pianka, E. R. (1974). Niche Overlap and Diffuse Competition. Proceedings of the National Academy of Sciences, 71(5), 2141–2145. https://doi.org/10.1073/pnas.71.5.2141",
      aspect = "Distributional overlap of proportional flowering schedules.",
      formula = "O_jk = sum_t p_j,t p_k,t / sqrt(sum_t p_j,t^2 * sum_t p_k,t^2).",
      input = "Wide numeric matrix/data frame; rows = census days; columns = individuals/species.",
      output = "List: proportional_matrix, pairwise_overlap, individual_mean_overlap, population_mean.",
      range = "0 to 1. Higher values indicate greater overlap."
    ),
    poole_rathcke = list(
      name = "Poole & Rathcke (1979) peak-spacing index",
      function_name = "poole_rathcke",
      citation = "Poole, R. W., & Rathcke, B. J. (1979). Regularity, Randomness, and Aggregation in Flowering Phenologies. Science, 203(4379), 470–471. https://doi.org/10.1126/science.203.4379.470",
      aspect = "Spacing, aggregation, or regularity of peak flowering dates.",
      formula = "P = variance of normalized distances among ordered peak dates and season boundaries; ratio = P / E(P).",
      input = "Wide numeric matrix/data frame; rows = census days; columns = species. Peak date is extracted as the first maximum day.",
      output = "List: peak_dates, ordered_peaks, distances, P, expected_P, ratio, p_value, interpretation.",
      range = "Ratio < 1 suggests regular spacing; ratio > 1 suggests aggregation."
    ),
    primack_Z = list(
      name = "Primack (1980) strong-flowering overlap index",
      function_name = "primack_Z",
      citation = "Primack, R. B. (1980). Variation in the Phenology of Natural Populations of Montane Shrubs in New Zealand. The Journal of Ecology, 68(3), 849–862. https://doi.org/10.2307/2259460",
      aspect = "Overlap among strong flowering periods.",
      formula = "c_ij = shared strong-flowering days / shorter strong-flowering duration; Z = mean pairwise c_ij.",
      input = "Wide numeric matrix/data frame. Strong flowering is defined relative to each individual's maximum flowering intensity.",
      output = "List: pairwise_table, Z, strong_matrix.",
      range = "0 to 1. Higher values indicate greater strong-period overlap."
    ),
    schoener_SI = list(
      name = "Schoener-style neighbourhood index applied to plant phenology by Hall et al., 2018 (https://doi.org/10.1080/15230430.2018.1548866)",
      function_name = "schoener_SI",
      citation = "Schoener, T. W. (1970). Nonsynchronous Spatial Overlap of Lizards in Patchy Habitats. Ecology, 51(3), 408–418. https://doi.org/10.2307/1935376",
      aspect = "Spatial neighbourhood synchrony; compares focal individuals with neighbours within a distance threshold.",
      formula = "Neighbourhood SI is calculated from temporal overlap/similarity among individuals within the specified spatial neighbourhood.",
      input = "Metadata with IDs and coordinates plus matching wide flowering matrix. Can use x/y or lon/lat.",
      output = "List: summary, SI_individual, mean_SI, distance_matrix, coordinates_xy, details.",
      range = "Usually 0 to 1 depending on implementation."
    ),
    shrotri_SI_temp = list(
      name = "Shrotri et al. (2025) temporal synchrony index",
      function_name = "shrotri_SI_temp",
      citation = "Shrotri, S., Kaur, S., Dhargalkar, R., Pathura Valappil, N., & Gowda, V. (2025). Mass Flowering and Flowering Asynchrony Characterize a Seasonal Herbaceous Community in the Western Ghats. Biotropica, 57(5), e70080. https://doi.org/10.1111/btp.70080",
      aspect = "Pairwise temporal overlap based on presence/absence.",
      formula = "SI_temp = F0 / ((F1 + F2) - F0), where F0 is shared flowering duration and F1,F2 are individual/species flowering durations.",
      input = "Wide numeric matrix/data frame; values > 0 are treated as presence.",
      output = "List: pairwise_SItemp, individual_SItemp, community_mean.",
      range = "0 to 1. Higher values indicate greater temporal overlap."
    ),
    shrotri_SI_abd = list(
      name = "Shrotri et al. (2025) abundance-weighted synchrony index",
      function_name = "shrotri_SI_abd",
      citation = "Shrotri, S., Kaur, S., Dhargalkar, R., Pathura Valappil, N., & Gowda, V. (2025). Mass Flowering and Flowering Asynchrony Characterize a Seasonal Herbaceous Community in the Western Ghats. Biotropica, 57(5), e70080. https://doi.org/10.1111/btp.70080",
      aspect = "Abundance-weighted co-flowering similarity.",
      formula = "SI_abd = mean_t(min(x_t, y_t) / max(x_t, y_t)) across shared flowering days.",
      input = "Wide numeric matrix/data frame with abundance/intensity values.",
      output = "List: pairwise_SIabd, individual_SIabd, community_mean.",
      range = "0 to 1. Higher values indicate similar abundance during co-flowering."
    ),
    wadgymar_Si = list(
      name = "Wadgymar et al. (2015) individual flowering-schedule synchrony",
      function_name = "wadgymar_Si",
      citation = "Wadgymar, S. M., Austen, E. J., Cumming, M. N., & Weis, A. E. (2015). Simultaneous pulsed flowering in a temperate legume: Causes and consequences of multimodality in the shape of floral display schedules. Journal of Ecology, 103(2), 316–327. https://doi.org/10.1111/1365-2745.12362",
      aspect = "Within-individual daily display unevenness/spread across flowering days.",
      formula = "S_i = sqrt(SS_i / T_i) * sqrt(D_i / (D_i - 1)), where SS_i is squared deviation from even daily flowering.",
      input = "Wide numeric matrix/data frame; rows = census days; columns = individuals.",
      output = "List: individual_results, Si_individual, mean_Si.",
      range = "Implementation bounded to 0 to 1."
    ),
    weis_sync = list(
      name = "Weis et al. (2014) phenotypic assortment and synchrony",
      function_name = "weis_sync",
      citation = "Weis, A. E., Nardone, E., & Fox, G. A. (2014). The strength of assortative mating for flowering date and its basis in individual variation in flowering schedule. Journal of Evolutionary Biology, 27(10), 2138–2151. https://doi.org/10.1111/jeb.12465",
      aspect = "Flowering synchrony and trait-based assortative mating from mating matrix structure.",
      formula = "PHI = M %*% P; S = lambda_1 / sum(lambda); rho = z_trait' PHI z_trait.",
      input = "Data frame/matrix with individuals in rows, one numeric trait column, and flowering-schedule columns.",
      output = "List: M_matrix, P_matrix, PHI_matrix, eigenvalues, S, rho, optional bootstrap intervals.",
      range = "S and rho depend on mating matrix and standardized trait structure."
    )
  )
}

#' Define and interpret a phenological synchrony index
#'
#' Returns a short definition, formula, citation, expected input structure, and
#' output structure for one or more synchrony indices implemented in PhenoSync.
#'
#' @param index Character index/function name. If `NULL`, all available indices
#'   are listed.
#' @param return Character. One of `"text"`, `"data.frame"`, or `"list"`.
#' @param fuzzy Logical. If `TRUE`, partial/fuzzy matching is allowed when it
#'   gives exactly one match.
#'
#' @return For `return = "text"`, prints an explanation and invisibly returns the
#'   matching list entry. For `return = "data.frame"`, returns a data frame. For
#'   `return = "list"`, returns a list.
#'
#' @examples
#' define_index()
#' define_index("pianka")
#' define_index("weis_sync", return = "data.frame")
#'
#' @export
define_index <- function(index = NULL,
                         return = c("text", "data.frame", "list"),
                         fuzzy = TRUE) {
  return <- match.arg(return)
  db <- pheno_index_database()

  if (is.null(index)) {
    if (return == "list") return(db)
    out <- data.frame(
      function_name = vapply(db, `[[`, character(1), "function_name"),
      name = vapply(db, `[[`, character(1), "name"),
      aspect = vapply(db, `[[`, character(1), "aspect"),
      input = vapply(db, `[[`, character(1), "input"),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    if (return == "data.frame") return(out)
    print(out, row.names = FALSE)
    return(invisible(out))
  }

  key <- index
  if (!key %in% names(db)) {
    if (fuzzy) {
      pattern <- tolower(index)
      hits <- names(db)[grepl(pattern, tolower(names(db))) |
                          grepl(pattern, tolower(vapply(db, `[[`, character(1), "name"))) |
                          grepl(pattern, tolower(vapply(db, `[[`, character(1), "function_name")))]
      if (length(hits) == 1) {
        key <- hits
      } else if (length(hits) > 1) {
        stop("Multiple matching indices found: ", paste(hits, collapse = ", "),
             ". Use the exact function name.", call. = FALSE)
      } else {
        stop("Unknown index: ", index, ". Run define_index() to list available indices.", call. = FALSE)
      }
    } else {
      stop("Unknown index: ", index, ". Run define_index() to list available indices.", call. = FALSE)
    }
  }

  item <- db[[key]]

  if (return == "list") {
    return(item)
  }

  if (return == "data.frame") {
    return(data.frame(
      field = names(item),
      value = unlist(item, use.names = FALSE),
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }

  cat("\n", item$name, "\n", sep = "")
  cat(strrep("=", nchar(item$name)), "\n\n", sep = "")
  cat("Function: ", item$function_name, "\n", sep = "")
  cat("Citation: ", item$citation, "\n\n", sep = "")
  cat("Definition / aspect captured:\n", item$aspect, "\n\n", sep = "")
  cat("Formula:\n", item$formula, "\n\n", sep = "")
  cat("Input data structure:\n", item$input, "\n\n", sep = "")
  cat("Output data structure:\n", item$output, "\n\n", sep = "")
  cat("Range / interpretation:\n", item$range, "\n\n", sep = "")
  invisible(item)
}

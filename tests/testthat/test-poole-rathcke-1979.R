test_that("Poole and Rathcke 1979 returns peak spacing output", {
  data(example_wide)
  res <- poole_rathcke(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c(
    "peak_dates", "ordered_peaks", "distances", "distance_table", "mean_distance",
    "P", "expected_P", "ratio", "chi_square", "df", "p_value", "season",
    "interpretation"
  ))
  expect_true(is.numeric(res$ratio))
  expect_true(res$ratio >= 0)
})

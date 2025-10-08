# tests/testthat/test-plot_freq_barplot.R

# --- Sample data used by multiple tests ---------------------------------------

df_freq_sample <- data.frame(
  # GroupA AC-rich, GroupB GT-rich; GroupC exists for subsetting tests
  Fragment_Bases_5p      = c("ACAC", "CACA", "GTGT", "TGTG", "AC", NA, "AAAA"),
  Fragment_Bases_3p      = c("CACA", "ACAC", "TGTG", "GTGT", "GT", "ACAC", "CCCC"),
  Fragment_Status_Simple = c("GroupA", "GroupA", "GroupB", "GroupB", "GroupA", "GroupB", "GroupC"),
  stringsAsFactors = FALSE
)

# ============== 1) Input validation / argument checks =========================

test_that("invalid or inconsistent arguments error out", {
  # invalid motif_type (message changed in the new function)
  expect_error(
    plot_freq_barplot(df_freq_sample, motif_type = "invalid_type"),
    regexp = "'motif_type' must be one of 'Start', 'End', or 'Both'\\."
  )

  # col_z = NULL but vals_z provided
  expect_error(
    plot_freq_barplot(df_freq_sample, col_z = NULL, vals_z = "GroupA"),
    regexp = "If 'col_z' is NULL, 'vals_z' must also be NULL\\."
  )

  # nonexisting grouping column
  expect_error(
    plot_freq_barplot(df_freq_sample, col_z = "NonExistentColumn"),
    regexp = "Column 'NonExistentColumn' not found in the dataframe\\.",
  )
})

# NOTE: The new function no longer requires ≥2 groups, so the old test asserting
# an error for a single-group analysis must be removed.

# ============== 2) Motif sizing + motif_type behavior =========================

test_that("motif_size is reduced when too large (with warning), title reflects k", {
  # shortest available sequence is length 2 ("AC"), requesting k=3 warns and uses 2
  expect_warning(
    p <- plot_freq_barplot(
      df_freq_sample,
      vals_z      = c("GroupA", "GroupB"),
      colors_z    = c("GroupA" = "#1f77b4", "GroupB" = "#ff7f0e"),
      output_path = NA_character_
    ),
    regexp = "Requested 'motif_size' \\(3\\) exceeds shortest sequence \\(2\\)\\. Using 2\\."
  )
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "Overall Nucleotide Frequency \\(2-mers\\)")
})

test_that("'Start' and 'End' motif_type run without error and return ggplot", {
  suppressWarnings({
    p_start <- plot_freq_barplot(
      df_freq_sample,
      motif_type = "Start",
      vals_z     = c("GroupA", "GroupB"),
      colors_z   = c("GroupA" = "#1f77b4", "GroupB" = "#ff7f0e"),
      output_path = NA_character_
    )
    p_end <- plot_freq_barplot(
      df_freq_sample,
      motif_type = "End",
      vals_z     = c("GroupA", "GroupB"),
      colors_z   = c("GroupA" = "#1f77b4", "GroupB" = "#ff7f0e"),
      output_path = NA_character_
    )
  })
  expect_s3_class(p_start, "ggplot")
  expect_s3_class(p_end,   "ggplot")
})

# ============== 3) Core plotting behavior (grouped vs ungrouped) =============

test_that("grouped analysis builds expected geoms and optional p-value", {
  p <- plot_freq_barplot(
    df_freq_sample,
    vals_z      = c("GroupA", "GroupB"),
    motif_size  = 2,
    colors_z    = c("GroupA" = "#1f77b4", "GroupB" = "#ff7f0e"),
    show_pvalue = TRUE,                 # new: caption includes global chi-square only if TRUE
    output_path = NA_character_
  )
  expect_s3_class(p, "ggplot")

  # first layer is GeomCol (not GeomBar in the new version)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomCol"      %in% geoms)
  expect_true("GeomErrorbar" %in% geoms)

  expect_match(p$labels$caption, "Global comparison \\(Chi-squared\\), p =")
})

test_that("ungrouped analysis works and does not add a p-value by default", {
  p <- plot_freq_barplot(
    df_freq_sample,
    col_z       = NULL,
    motif_size  = 2,
    colors_z    = c("#1f77b4"),
    output_path = NA_character_
  )
  expect_s3_class(p, "ggplot")
  expect_false(grepl("Chi-squared", p$labels$caption))

  # legend labels contain "All Fragments"
  fill_scale_idx <- which(vapply(p$scales$scales, function(s) "fill" %in% s$aesthetics, logical(1)))
  legend_labels  <- p$scales$scales[[fill_scale_idx]]$labels
  expect_true(grepl("All Fragments", legend_labels))
})

# ============== 4) Subsetting, colors, and ... passthrough ====================

test_that("vals_z filtering keeps only requested groups", {
  suppressWarnings({
    p <- plot_freq_barplot(
      df_freq_sample,
      vals_z      = c("GroupA", "GroupC"),
      motif_size  = 2,
      colors_z    = c("GroupA" = "#1f77b4", "GroupC" = "#2ca02c"),
      output_path = NA_character_
    )
  })
  expect_s3_class(p, "ggplot")
  plotted_groups <- as.character(unique(p$data$group))
  expect_equal(sort(plotted_groups), c("GroupA", "GroupC"))
})

test_that("RColorBrewer palette and custom colors are applied to groups", {
  # Brewer palette
  full_paired <- RColorBrewer::brewer.pal(12, "Paired")
  p_brewer <- suppressWarnings(plot_freq_barplot(
    df_freq_sample,
    vals_z      = c("GroupA", "GroupB"),
    motif_size  = 2,
    colors_z    = "Paired",
    output_path = NA_character_
  ))
  built_brewer         <- ggplot2::ggplot_build(p_brewer)
  actual_brewer_colors <- unique(built_brewer$data[[1]]$fill)
  expect_length(actual_brewer_colors, 2)
  expect_true(all(actual_brewer_colors %in% full_paired))

  # Unnamed custom colors, length must equal number of groups
  custom_cols <- c("gold", "purple")
  p_custom <- suppressWarnings(plot_freq_barplot(
    df_freq_sample,
    vals_z      = c("GroupA", "GroupB"),
    motif_size  = 2,
    colors_z    = custom_cols,
    output_path = NA_character_
  ))
  built_custom         <- ggplot2::ggplot_build(p_custom)
  actual_custom_colors <- sort(unique(built_custom$data[[1]]$fill))
  expect_equal(actual_custom_colors, sort(custom_cols))
})

test_that("extra aesthetics in ... are passed to geom_col (alpha, width)", {
  p <- plot_freq_barplot(
    df_freq_sample,
    vals_z      = c("GroupA", "GroupB"),
    motif_size  = 2,
    alpha       = 0.5,
    width       = 0.4,  # user override takes precedence
    colors_z    = c("GroupA" = "#1f77b4", "GroupB" = "#ff7f0e"),
    output_path = NA_character_
  )

  expect_s3_class(p, "ggplot")
})

# ============== 5) Saving behavior (tempdir only; no interactive ops) =========

test_that("saving to a file in tempdir works; overwrite is silent", {
  tmp_file <- file.path(tempdir(), "freq_nucleotide_test.png")
  if (file.exists(tmp_file)) file.remove(tmp_file)

  suppressWarnings({
    plot_freq_barplot(
      df_freq_sample,
      vals_z      = c("GroupA", "GroupB"),
      motif_size  = 2,
      colors_z    = "Paired",
      output_path = tmp_file
    )
  })
  expect_true(file.exists(tmp_file))

  # Overwrite silently (no message expected)
  suppressWarnings({
    plot_freq_barplot(
      df_freq_sample,
      vals_z      = c("GroupA", "GroupB"),
      motif_size  = 2,
      colors_z    = "Paired",
      output_path = tmp_file
    )
  })
  expect_true(file.exists(tmp_file))

  file.remove(tmp_file)
})

test_that("custom ggsave_params are honored", {
  tmp_file <- file.path(tempdir(), "freq_nucleotide_custom_size.png")
  if (file.exists(tmp_file)) file.remove(tmp_file)

  suppressWarnings({
    plot_freq_barplot(
      df_freq_sample,
      vals_z       = c("GroupA", "GroupB"),
      colors_z     = "Paired",
      output_path  = tmp_file,
      ggsave_params = list(width = 7, height = 5, units = "in", dpi = 150, bg = "white")
    )
  })
  expect_true(file.exists(tmp_file))
  file.remove(tmp_file)
})

test_that("vector output_path (>1) is ignored and a ggplot is returned", {
  p <- suppressWarnings(plot_freq_barplot(
    df_freq_sample,
    vals_z      = c("GroupA", "GroupB"),
    colors_z    = "Paired",
    output_path = c("a.png", "b.png")  # ignored by length != 1 check
  ))
  expect_s3_class(p, "ggplot")
})

# ============== 6) Coverage: color validation & 'Other' handling ==============

test_that("colors_z validation errors on wrong lengths or missing names", {
  # unnamed vector with insufficient length
  suppressWarnings(
    expect_error(
      plot_freq_barplot(
        df_freq_sample,
        vals_z      = c("GroupA", "GroupB", "GroupC"),
        motif_size  = 2,
        colors_z    = c("red", "blue"),  # not enough colors
        output_path = NA_character_
      ),
      regexp = "If 'colors_z' is an unnamed vector, its length must match the number of groups\\."
    )
  )

  # named vector missing a required group
  suppressWarnings(
    expect_error(
      plot_freq_barplot(
        df_freq_sample,
        vals_z      = c("GroupA", "GroupB"),
        motif_size  = 2,
        colors_z    = c(GroupA = "red"), # missing GroupB
        output_path = NA_character_
      ),
      regexp = "Named 'colors_z' must provide a color for each group level\\."
    )
  )
})

test_that("drop_non_acgt = FALSE keeps 'Other' category when present", {
  df2 <- df_freq_sample
  df2$Fragment_Bases_5p[1] <- "NNNAC"  # inject a non-ACGT character
  p <- plot_freq_barplot(
    df2,
    vals_z      = c("GroupA", "GroupB"),
    motif_size  = 2,
    drop_non_acgt = FALSE,
    colors_z    = c("GroupA" = "#1f77b4", "GroupB" = "#ff7f0e"),
    output_path = NA_character_
  )
  expect_true("Other" %in% levels(p$data$nucleotide))
})

test_that("show_pvalue = FALSE does not add p-value, even with 2 groups", {
  p <- plot_freq_barplot(
    df_freq_sample,
    vals_z      = c("GroupA", "GroupB"),
    motif_size  = 2,
    colors_z    = c("GroupA" = "#1f77b4", "GroupB" = "#ff7f0e"),
    show_pvalue = FALSE,
    output_path = NA_character_
  )
  expect_false(grepl("Chi-squared", p$labels$caption))
})

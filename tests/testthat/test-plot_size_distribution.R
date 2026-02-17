# tests/testthat/test-plot_size_distribution.R

# --- Sample data used by multiple tests ---------------------------------------

set.seed(42)
df_sample <- data.frame(
  Fragment_Size = c(rnorm(100, 150, 20), rnorm(100, 320, 30), rnorm(50, 480, 40), NA),
  Fragment_Status_Simple = c(rep("Mono", 100), rep("Di", 100), rep("Tri", 50), "Mono"),
  Other_Info = sample(letters, 251, replace = TRUE),
  stringsAsFactors = FALSE
)
# Add an NA to the grouping column as well.
df_sample$Fragment_Status_Simple[1] <- NA


# ============== 1) Input validation / argument checks =========================

test_that("errors for invalid or missing columns", {
  # Missing size column
  expect_error(
    plot_size_distribution(df_sample, size_col = "NonExistentSizeCol"),
    regexp = "Size column 'NonExistentSizeCol' not found in the dataframe\\."
  )

  # Missing grouping column (use fixed string match)
  expect_error(
    plot_size_distribution(df_sample, col_z = "NonExistentGroupCol"),
    regexp = "Column 'NonExistentGroupCol' not found in the dataframe.",
    fixed = TRUE
  )

  # Non-numeric size column
  df_bad_type <- df_sample
  df_bad_type$Fragment_Size <- as.character(df_bad_type$Fragment_Size)
  expect_error(
    plot_size_distribution(df_bad_type, size_col = "Fragment_Size"),
    regexp = "Size column 'Fragment_Size' must be numeric\\."
  )
})

test_that("errors for logical inconsistencies in arguments", {
  # col_z = NULL but vals_z is provided
  expect_error(
    plot_size_distribution(df_sample, col_z = NULL, vals_z = c("Mono", "Di")),
    regexp = "If 'col_z' is NULL, 'vals_z' must also be NULL\\."
  )

  # Both histogram and density disabled
  expect_error(
    plot_size_distribution(df_sample, show_histogram = FALSE, show_density = FALSE),
    regexp = "At least one of 'show_histogram' or 'show_density' must be TRUE\\."
  )
})

test_that("errors when filtering leaves no data", {
  expect_error(
    plot_size_distribution(df_sample, vals_z = "NonExistentValue"),
    regexp = "No data remains after filtering. Check 'col_z' and 'vals_z'\\."
  )
})


# ============== 2) Core functionality and outputs ============================

test_that("returns a ggplot object with default settings", {
  p <- plot_size_distribution(df_sample)
  expect_s3_class(p, "ggplot")

  # Default includes density layer
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomDensity" %in% geoms)
})

test_that("ungrouped analysis (col_z = NULL) works", {
  p <- plot_size_distribution(df_sample, col_z = NULL)
  expect_s3_class(p, "ggplot")

  # One group drawn
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$group)), 1)

  # Legend labels: find factor column with "(N=...)" labels and assert pooled label
  fac_cols <- names(p$data)[vapply(p$data, is.factor, logical(1))]
  all_levels <- unlist(lapply(fac_cols, function(nm) levels(p$data[[nm]])))
  expect_true(any(grepl("^All Fragments \\(N=\\d+\\)$", all_levels)))
})

test_that("filtering by 'vals_z' keeps only requested groups", {
  groups_to_show <- c("Mono", "Di")
  p <- plot_size_distribution(df_sample, vals_z = groups_to_show)
  expect_s3_class(p, "ggplot")

  # Legend labels come from factor levels in the data ("Group (N=...)")
  fac_cols <- names(p$data)[vapply(p$data, is.factor, logical(1))]
  all_levels <- unlist(lapply(fac_cols, function(nm) levels(p$data[[nm]])))
  expect_true(all(grepl("^(Mono|Di) \\(N=\\d+\\)$", all_levels)))
  expect_false(any(grepl("^Tri \\(N=\\d+\\)$", all_levels)))
})


# ============== 3) Specific parameterization behaviors =======================

test_that("plot layers are correctly added or removed", {
  # Histogram only
  p_hist <- plot_size_distribution(df_sample, show_histogram = TRUE, show_density = FALSE)
  geoms_hist <- vapply(p_hist$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomBar" %in% geoms_hist) # geom_histogram
  expect_false("GeomDensity" %in% geoms_hist)
  expect_equal(p_hist$labels$y, "Density") # y label for histogram is density (via after_stat)

  # Both histogram and density
  p_both <- plot_size_distribution(df_sample, show_histogram = TRUE, show_density = TRUE)
  geoms_both <- vapply(p_both$layers, function(l) class(l$geom)[1], character(1))
  expect_true(all(c("GeomBar", "GeomDensity") %in% geoms_both))

  # Turn off nucleosome peaks
  p_no_peaks <- plot_size_distribution(df_sample, show_nuc_peaks = FALSE)
  geoms_no_peaks <- vapply(p_no_peaks$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomVline" %in% geoms_no_peaks)
})

test_that("custom arguments are forwarded to geom layers", {
  # Density linewidth override (robust across ggplot2 versions using size/linewidth)
  p_density <- plot_size_distribution(df_sample, density_args = list(linewidth = 2))
  line_param <- p_density$layers[[1]]$aes_params$linewidth
  if (is.null(line_param)) line_param <- p_density$layers[[1]]$aes_params$size
  expect_equal(line_param, 2)

  # Histogram-only alpha
  p_hist <- plot_size_distribution(
    df_sample,
    show_histogram = TRUE,
    show_density   = FALSE,
    histo_args     = list(alpha = 0.5)
  )
  expect_equal(p_hist$layers[[1]]$aes_params$alpha, 0.5)
})

test_that("coloring schemes are applied correctly (manual and Brewer)", {
  # Helper to convert any color spec to hex
  to_hex <- function(cols) {
    m <- grDevices::col2rgb(cols)
    grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
  }

  # Manual unnamed colors (order matches groups in vals_z)
  custom_colors <- c("red", "green", "blue")
  p_custom <- plot_size_distribution(
    df_sample,
    vals_z   = c("Mono", "Di", "Tri"),
    colors_z = custom_colors
  )
  built_custom <- ggplot2::ggplot_build(p_custom)
  used_custom <- unique(built_custom$data[[1]]$colour)
  # Compare hex-to-hex to be robust across devices/versions
  expect_setequal(tolower(to_hex(used_custom)), tolower(to_hex(custom_colors)))

  # RColorBrewer palette
  p_brewer <- plot_size_distribution(
    df_sample,
    vals_z   = c("Mono", "Di", "Tri"),
    colors_z = "Set1"
  )
  built_brewer <- ggplot2::ggplot_build(p_brewer)
  used_brewer <- unique(built_brewer$data[[1]]$colour)
  expected_brewer <- RColorBrewer::brewer.pal(3, "Set1")
  expect_setequal(tolower(to_hex(used_brewer)), tolower(to_hex(expected_brewer)))
})

test_that("axis limits are applied via coord_cartesian", {
  custom_limits <- c(100, 500)
  p <- plot_size_distribution(df_sample, x_limits = custom_limits)
  coord_limits <- ggplot2::ggplot_build(p)$layout$coord$limits$x
  expect_equal(coord_limits, custom_limits)
})

test_that("groups with <2 points are dropped only in density mode", {
  # Ensure at least one group has n >= 2 to avoid full drop in density mode
  tiny <- data.frame(
    Fragment_Size = c(170, 330, 500, 180, 175),
    Fragment_Status_Simple = c("A", "B", "C", "Tiny", "A")
  )

  # Density mode: 'Tiny' (n=1) should be removed, 'A' (n=2) retained
  p_dens <- plot_size_distribution(tiny, show_density = TRUE, show_histogram = FALSE)
  fac_cols_d <- names(p_dens$data)[vapply(p_dens$data, is.factor, logical(1))]
  levels_d <- unique(unlist(lapply(fac_cols_d, function(nm) levels(p_dens$data[[nm]]))))
  expect_false(any(grepl("^Tiny \\(N=", levels_d)))
  expect_true(any(grepl("^A \\(N=", levels_d)))

  # Histogram-only: 'Tiny' should remain
  p_hist <- plot_size_distribution(tiny, show_density = FALSE, show_histogram = TRUE)
  fac_cols_h <- names(p_hist$data)[vapply(p_hist$data, is.factor, logical(1))]
  levels_h <- unique(unlist(lapply(fac_cols_h, function(nm) levels(p_hist$data[[nm]]))))
  expect_true(any(grepl("^Tiny \\(N=", levels_h)))
})

test_that("named color vector using base group names maps to '(N=...)' labels", {
  to_hex <- function(cols) {
    m <- grDevices::col2rgb(cols)
    grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
  }

  cols_named <- c(Mono = "black", Di = "grey50", Tri = "grey80")
  p <- plot_size_distribution(
    df_sample,
    vals_z   = c("Mono", "Di", "Tri"),
    colors_z = cols_named
  )
  built <- ggplot2::ggplot_build(p)
  used <- unique(built$data[[1]]$colour)
  expected <- c("black", "grey50", "grey80")
  expect_true(all(tolower(to_hex(used)) %in% tolower(to_hex(expected))))
})

test_that("custom title is honored", {
  p <- plot_size_distribution(df_sample, title = "My Custom Title")
  expect_identical(p$labels$title, "My Custom Title")
})


# ============== 4) Saving behavior (tempdir only; no interactive ops) =========

test_that("no output_path returns a ggplot object (no saving branch)", {
  p <- suppressWarnings(plot_size_distribution(df_fragments = df_sample))
  expect_s3_class(p, "ggplot")
})

test_that("valid output_path saves the plot and honors ggsave_params", {
  temp_dir <- tempdir()
  out_file <- file.path(temp_dir, "size_dist_save.png")
  if (file.exists(out_file)) file.remove(out_file)

  suppressWarnings({
    plot_size_distribution(
      df_fragments       = df_sample,
      show_histogram     = TRUE, # exercise histogram path too
      histogram_binwidth = 10,
      output_path        = out_file,
      ggsave_params      = list(width = 7, height = 5, units = "in", dpi = 150, bg = "white")
    )
  })
  expect_true(file.exists(out_file))
  file.remove(out_file)
})

test_that("overwriting an existing file is silent", {
  temp_dir <- tempdir()
  out_file <- file.path(temp_dir, "size_dist_overwrite.png")

  # First save
  suppressWarnings({
    plot_size_distribution(df_fragments = df_sample, output_path = out_file)
  })
  expect_true(file.exists(out_file))

  # Second save (no message expected)
  suppressWarnings({
    plot_size_distribution(df_fragments = df_sample, output_path = out_file)
  })
  expect_true(file.exists(out_file))

  file.remove(out_file)
})

test_that("invalid output_path inputs are ignored and a ggplot is returned", {
  # length > 1
  p_vec <- suppressWarnings(plot_size_distribution(
    df_fragments = df_sample,
    output_path  = c("path1", "path2")
  ))
  expect_s3_class(p_vec, "ggplot")

  # non-character scalar
  p_num <- suppressWarnings(plot_size_distribution(
    df_fragments = df_sample,
    output_path  = 123
  ))
  expect_s3_class(p_num, "ggplot")
})

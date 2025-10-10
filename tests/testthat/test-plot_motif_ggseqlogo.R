.capture_warnings <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(
    expr,
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  msgs
}

# --- Sample data used by multiple tests ---------------------------------------

df_motif_sample <- data.frame(
  Fragment_Bases_5p = c(
    "ACGT", "ACGG", "AGTA", "CCTT",        # Short
    "TGCAT", "TGAAC", "TGCAC", NA,         # Long
    "GATTACA", "GATTANA",                  # Mixed (with 'N' later in string)
    "GC"                                   # Short entry (too short for k=3)
  ),
  Fragment_Bases_3p = c(
    "TGCA", "CCGG", "TACT", "AAGG",        # Short
    "ATGCA", "GTTCA", "GTGCA", "A",        # Long (one too short)
    "TACAATG", "TANAATG",                  # Mixed (one has N)
    "CG"                                   # Short
  ),
  Fragment_Status_Simple = c(
    rep("Short", 4),
    rep("Long", 4),
    rep("Mixed", 2),
    "Short"
  ),
  stringsAsFactors = FALSE
)

# Helpers ----------------------------------------------------------------------

.find_geom_idx <- function(p, geom_class) {
  which(vapply(p$layers, function(l) inherits(l$geom, geom_class), logical(1)))
}

.to_hex <- function(cols) {
  if (length(cols) == 0) return(character(0))
  m <- grDevices::col2rgb(cols)
  apply(m, 2L, function(v) grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255))
}

# ============== 1) Input validation / argument checks =========================

test_that("errors for invalid arguments", {
  # col_z = NULL but vals_z provided
  expect_error(
    plot_qqseqlogo_meme(df_motif_sample, col_z = NULL, vals_z = "Short"),
    regexp = "If 'col_z' is NULL, 'vals_z' must also be NULL\\."
  )

  # Non-existent grouping column (use fixed match to avoid regex quirks)
  expect_error(
    plot_qqseqlogo_meme(df_motif_sample, col_z = "NonExistentCol"),
    regexp = "Column 'NonExistentCol' not found in the dataframe.",
    fixed  = TRUE
  )

  # Invalid motif_type
  expect_error(
    plot_qqseqlogo_meme(df_motif_sample, motif_type = "InvalidType"),
    regexp = "'motif_type' must be 'Start', 'End', or 'Both'\\."
  )
})

# ============== 2) Data filtering and grouping ===============================

test_that("ungrouped and grouped analyses return ggplot and facet labels", {
  # Ungrouped (default 'Both'); use k=2 for deterministic counts
  p_ungrouped <- suppressWarnings(
    plot_qqseqlogo_meme(df_motif_sample, col_z = NULL, motif_size = 2)
  )
  expect_s3_class(p_ungrouped, "ggplot")

  built <- ggplot2::ggplot_build(p_ungrouped)
  layout_df <- built$layout$layout
  group_col <- grep("seq_group", names(layout_df), value = TRUE)
  expect_true(length(group_col) == 1)
  # Expect label form, not a hard-coded count
  expect_true(grepl("^All Fragments \\(N=\\d+\\)$", as.character(layout_df[[group_col]])))

  # Grouped with two groups; accept either exact N or implementation-dependent N
  p_grouped <- suppressWarnings(
    plot_qqseqlogo_meme(df_motif_sample, vals_z = c("Short", "Long"), motif_size = 3)
  )
  expect_s3_class(p_grouped, "ggplot")

  built_g <- ggplot2::ggplot_build(p_grouped)
  layout_g <- built_g$layout$layout
  group_col_g <- grep("seq_group", names(layout_g), value = TRUE)
  expect_true(length(group_col_g) == 1)

  labs <- as.character(unique(layout_g[[group_col_g]]))
  # Ensure both groups appear with counts, but don't pin the exact N
  expect_true(any(grepl("^Short \\(N=\\d+\\)$", labs)))
  expect_true(any(grepl("^Long \\(N=\\d+\\)$",  labs)))
})

# ============== 3) Motif extraction and handling =============================

test_that("motif_type 'Start', 'End', and 'Both' produce expected structure", {
  # Start (k=2)
  p_start <- plot_qqseqlogo_meme(df_motif_sample, motif_type = "Start", motif_size = 2, vals_z = "Short")
  expect_s3_class(p_start, "ggplot")
  b_start <- ggplot2::ggplot_build(p_start)
  k_start <- round(max(b_start$data[[1]]$x))
  expect_equal(k_start, 2)
  expect_false("GeomVline" %in% sapply(p_start$layers, function(l) class(l$geom)[1]))

  # End (k=2)
  p_end <- plot_qqseqlogo_meme(df_motif_sample, motif_type = "End", motif_size = 2, vals_z = "Short")
  expect_s3_class(p_end, "ggplot")
  b_end <- ggplot2::ggplot_build(p_end)
  k_end <- round(max(b_end$data[[1]]$x))
  expect_equal(k_end, 2)

  # Both (k=2 per side ⇒ 2 + 1 + 2)
  p_both <- suppressWarnings(plot_qqseqlogo_meme(df_motif_sample, motif_size = 2, vals_z = "Short"))
  expect_s3_class(p_both, "ggplot")
  b_both <- ggplot2::ggplot_build(p_both)
  k_both <- round(max(b_both$data[[1]]$x))
  expect_equal(k_both, 5)
  expect_true("GeomVline" %in% sapply(p_both$layers, function(l) class(l$geom)[1]))
})

test_that("warnings: motif_size capped when too large; short sequences per-group noted", {
  # Start: we should see the global capping message
  wrn_start <- .capture_warnings(
    plot_qqseqlogo_meme(df_motif_sample, motif_type = "Start", motif_size = 3)
  )
  expect_true(
    any(grepl("Requested 'motif_size'", wrn_start)),
    info = paste("Warnings seen (Start):", paste(wrn_start, collapse = " | "))
  )

  # End: depending on data, either capping happens OR (in other datasets) short-sequence removal
  wrn_end <- .capture_warnings(
    plot_qqseqlogo_meme(df_motif_sample, motif_type = "End", motif_size = 3, vals_z = "Long")
  )
  expect_true(
    any(grepl("Requested 'motif_size'", wrn_end)) ||
      any(grepl("removed \\d+ 3' sequences shorter than motif_size", wrn_end)),
    info = paste("Warnings seen (End):", paste(wrn_end, collapse = " | "))
  )
})

test_that("behaves gracefully when nothing remains to plot (returns empty ggplot)", {
  df_only_n <- data.frame(
    Fragment_Bases_5p = "NNN",
    Fragment_Bases_3p = "NNN",
    Fragment_Status_Simple = "A"
  )
  p_empty <- suppressWarnings(plot_qqseqlogo_meme(df_only_n))
  expect_s3_class(p_empty, "ggplot")
  expect_identical(p_empty$labels$title, "No data to display")
})

# ============== 4) Axis labels and annotation =================================

test_that("custom axis labels reflect motif_type", {
  # Start, k=3 (will cap if needed)
  p_start <- suppressWarnings(
    plot_qqseqlogo_meme(df_motif_sample, motif_type = "Start", motif_size = 3, vals_z = "Short")
  )
  b_start <- ggplot2::ggplot_build(p_start)
  idx_text <- .find_geom_idx(p_start, "GeomText")
  lab_start <- b_start$data[[idx_text]]$label
  expect_equal(lab_start, as.character(1:2))

  # End, k=3 (will cap if needed)
  p_end <- suppressWarnings(
    plot_qqseqlogo_meme(df_motif_sample, motif_type = "End", motif_size = 3, vals_z = "Short")
  )
  b_end <- ggplot2::ggplot_build(p_end)
  idx_text_e <- .find_geom_idx(p_end, "GeomText")
  lab_end <- b_end$data[[idx_text_e]]$label
  expect_equal(lab_end, as.character((-2):-1))

  # Both, k=2 per side
  p_both <- plot_qqseqlogo_meme(df_motif_sample, motif_type = "Both", motif_size = 2, vals_z = "Short")
  b_both <- ggplot2::ggplot_build(p_both)
  idx_text_b <- .find_geom_idx(p_both, "GeomText")
  lab_both <- b_both$data[[idx_text_b]]$label
  expect_equal(lab_both, as.character(c(1:2, "", (-2):-1)))
})

# ============== 5) Color schemes ==============================================

test_that("named nucleotide colors are applied; '-' is ignored (logo separator)", {
  custom_cols <- c(A = "blue", C = "red", G = "green", T = "purple", "-" = "grey")
  p <- suppressWarnings(plot_qqseqlogo_meme(df_motif_sample, motif_size = 1, colors_z = custom_cols))
  b <- ggplot2::ggplot_build(p)
  used <- unique(b$data[[1]]$fill)
  exp <- c("blue", "red", "green", "purple")
  expect_setequal(.to_hex(used), .to_hex(exp))
})

test_that("RColorBrewer palette name uses four A/C/G/T colors", {
  p <- suppressWarnings(plot_qqseqlogo_meme(df_motif_sample, motif_size = 1, colors_z = "Dark2"))
  b <- ggplot2::ggplot_build(p)
  used <- unique(b$data[[1]]$fill)
  exp <- RColorBrewer::brewer.pal(4, "Dark2")
  expect_setequal(.to_hex(used), .to_hex(exp))
})

test_that("named vector must include A/C/G/T; unnamed must have 4 colors; invalid palette name errors meaningfully", {
  expect_error(
    suppressWarnings(plot_qqseqlogo_meme(df_motif_sample, motif_size = 1, colors_z = c(A="blue", C="red"))),
    regexp = "Custom colors are incomplete: missing G, T\\."
  )

  expect_error(
    suppressWarnings(plot_qqseqlogo_meme(df_motif_sample, motif_size = 1, colors_z = c("red","green","blue"))),
    regexp = "An unnamed vector must contain exactly 4 colors \\(received: 3\\)\\."
  )

  # Non-Brewer single string falls through to "unnamed vector of length 1" case
  expect_error(
    suppressWarnings(plot_qqseqlogo_meme(df_motif_sample, motif_size = 1, colors_z = "NotAPalette")),
    regexp = "An unnamed vector must contain exactly 4 colors \\(received: 1\\)\\."
  )
})

test_that("user-supplied col_scheme in ... takes precedence over colors_z", {
  cs <- ggseqlogo::make_col_scheme(
    chars = c("A","C","G","T"),
    cols  = c("#111111", "#222222", "#333333", "#444444")
  )
  p <- plot_qqseqlogo_meme(
    df_motif_sample,
    motif_size = 1,
    colors_z   = "Set1",
    col_scheme = cs
  )
  b <- ggplot2::ggplot_build(p)
  used <- sort(unique(.to_hex(b$data[[1]]$fill)))
  expect_setequal(used, sort(.to_hex(c("#111111", "#222222", "#333333", "#444444"))))
})

# ============== 6) Saving behavior (tempdir only; no interactive ops) =========

test_that("no output_path returns a ggplot object (no saving branch)", {
  p <- suppressWarnings(plot_qqseqlogo_meme(df_fragments = df_motif_sample))
  expect_s3_class(p, "ggplot")
})

test_that("valid output_path saves the plot and honors ggsave_params", {
  temp_dir <- tempdir()
  out_file <- file.path(temp_dir, "qqseqlogo_save.png")
  if (file.exists(out_file)) file.remove(out_file)

  suppressWarnings({
    plot_qqseqlogo_meme(
      df_fragments  = df_motif_sample,
      motif_type    = "Both",
      output_path   = out_file,
      ggsave_params = list(width = 7, height = 5, units = "in", dpi = 150, bg = "white")
    )
  })
  expect_true(file.exists(out_file))
  file.remove(out_file)
})

test_that("overwriting an existing file is silent (no message required)", {
  temp_dir <- tempdir()
  out_file <- file.path(temp_dir, "qqseqlogo_overwrite.png")

  suppressWarnings({
    plot_qqseqlogo_meme(df_fragments = df_motif_sample, output_path = out_file)
  })
  expect_true(file.exists(out_file))

  suppressWarnings({
    plot_qqseqlogo_meme(df_fragments = df_motif_sample, output_path = out_file)
  })
  expect_true(file.exists(out_file))

  file.remove(out_file)
})

test_that("invalid output_path inputs are ignored and a ggplot is returned", {
  p_vec <- suppressWarnings(plot_qqseqlogo_meme(
    df_fragments = df_motif_sample,
    output_path  = c("path1", "path2")
  ))
  expect_s3_class(p_vec, "ggplot")

  p_num <- suppressWarnings(plot_qqseqlogo_meme(
    df_fragments = df_motif_sample,
    output_path  = 123
  ))
  expect_s3_class(p_num, "ggplot")
})

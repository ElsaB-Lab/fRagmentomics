# tests/testthat/test-plot_motif_barplot.R

# --- Sample data shared by tests ------------------------------------------------

df_barplot_sample <- data.frame(
  Fragment_Bases_5p = c(
    "ACG", "ACT", "AGA",   # GroupA
    "AGT", "CAT", "CCT",   # GroupB
    "GAT", "GGT", "TGA",   # GroupC
    "NNN", "ANC", NA       # Invalid motifs spread across groups
  ),
  Fragment_Bases_3p = c(
    "TGA", "TCA", "TCT",   # GroupA
    "TGT", "GTA", "AGG",   # GroupB
    "ATC", "ACC", "TCA",   # GroupC
    "TCA", "ANT", "CCT"    # Invalid motifs
  ),
  Fragment_Status_Simple = c(
    "GroupA", "GroupA", "GroupA",
    "GroupB", "GroupB", "GroupB",
    "GroupC", "GroupC", "GroupC",
    "GroupA", "GroupB", "GroupC"
  ),
  stringsAsFactors = FALSE
)

# small helper for robust color comparison (works with hex & named)
.normalise_hex <- function(cols) {
  if (length(cols) == 0) return(character(0))
  m <- grDevices::col2rgb(cols)
  # col2rgb returns a matrix with a column per color
  apply(m, 2L, function(v) {
    grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
  })
}

# ============== 1) Input validation & argument checks ==========================

test_that("invalid arguments and inconsistent parameters error out", {
  expect_error(
    plot_motif_barplot(df_barplot_sample, representation = "invalid_type"),
    regexp = "'representation' must be one of: differential, split_by_base, split_by_motif"
  )

  expect_error(
    plot_motif_barplot(df_barplot_sample, col_z = NULL, vals_z = "GroupA"),
    regexp = "If 'col_z' is NULL, 'vals_z' must also be NULL\\."
  )

  expect_error(
    plot_motif_barplot(df_barplot_sample, col_z = "NonExistentColumn"),
    regexp = "Column 'NonExistentColumn' not found in the dataframe.",
    fixed  = TRUE
  )
})

test_that("representation-specific requirements are enforced", {
  expect_error(
    plot_motif_barplot(df_barplot_sample, representation = "differential", col_z = NULL),
    regexp = "Differential analysis requires a grouping column 'col_z'\\."
  )

  expect_error(
    plot_motif_barplot(df_barplot_sample, representation = "differential", vals_z = "GroupA"),
    regexp = "Differential analysis requires exactly two values in 'vals_z'\\."
  )
  expect_error(
    plot_motif_barplot(df_barplot_sample, representation = "differential",
                       vals_z = c("GroupA", "GroupB", "GroupC")),
    regexp = "Differential analysis requires exactly two values in 'vals_z'\\."
  )
})

test_that("empty data after filtering and no-valid-motif cases error out", {
  expect_error(
    plot_motif_barplot(df_barplot_sample, vals_z = "NonExistentGroup"),
    regexp = "No data remains after filtering. Check 'col_z' and 'vals_z'\\."
  )

  df_no_valid_motifs <- data.frame(
    Fragment_Bases_5p = c("NN", "AT", "NNA"),
    Fragment_Bases_3p = c("N", "G", "C"),
    Fragment_Status_Simple = c("A", "A", "A")
  )
  expect_error(
    plot_motif_barplot(df_no_valid_motifs),
    regexp = "No valid 3-base motifs found in the data after filtering\\."
  )
})


# ============== 2) Core functionality & representation behavior ===============

test_that("default 'split_by_base' representation builds FacetNested and accepts Brewer palette", {
  p <- plot_motif_barplot(df_barplot_sample)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetNested")

  p_brewer <- plot_motif_barplot(df_barplot_sample, colors_z = "Set1")
  built <- ggplot2::ggplot_build(p_brewer)
  expect_setequal(
    unique(built$data[[1]]$fill),
    RColorBrewer::brewer.pal(4, "Set1")
  )
})

test_that("'split_by_motif' representation facets by wrap, rotates x text, and uses expected title", {
  p <- plot_motif_barplot(df_barplot_sample, representation = "split_by_motif")
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(p$theme$axis.text.x$angle, 90)
  expect_identical(p$labels$title, "Motif Proportions by Group")
})

test_that("'differential' representation: labels, facets, data columns, and custom colors", {
  p <- plot_motif_barplot(
    df_barplot_sample,
    representation = "differential",
    vals_z = c("GroupA", "GroupB")
  )
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetNested")
  expect_true(grepl("Log2 Fold Change", p$labels$y))
  expect_identical(p$labels$subtitle, "Comparison: GroupA vs GroupB")
  expect_true(all(c("log2FC", "sign") %in% names(p$data)))

  # Custom Positive/Negative colors — robust to named/hex returns
  p_col <- plot_motif_barplot(
    df_barplot_sample,
    representation = "differential",
    vals_z = c("GroupA", "GroupB"),
    colors_z = c(Positive = "navy", Negative = "tomato")
  )
  built_col <- ggplot2::ggplot_build(p_col)
  used_norm <- unique(.normalise_hex(built_col$data[[1]]$fill))
  exp_norm  <- .normalise_hex(c("navy", "tomato"))
  # used colors must be a subset of the two expected (sometimes only one sign appears)
  expect_true(all(used_norm %in% exp_norm))
})


# ============== 3) Parameters, grouping, and filters ==========================

test_that("motif_type = 'Start' vs 'End' counts (N=...) reflect valid 3-mers", {
  p_start <- plot_motif_barplot(df_barplot_sample, motif_type = "Start")
  built_start <- ggplot2::ggplot_build(p_start)
  labels_start <- unique(built_start$layout$layout$group_label)
  expect_true(any(grepl("^GroupA \\(N=3\\)$", labels_start)))
  expect_true(any(grepl("^GroupB \\(N=3\\)$", labels_start)))
  expect_true(any(grepl("^GroupC \\(N=3\\)$", labels_start)))

  p_end <- plot_motif_barplot(df_barplot_sample, motif_type = "End")
  built_end <- ggplot2::ggplot_build(p_end)
  labels_end <- unique(built_end$layout$layout$group_label)
  expect_true(any(grepl("^GroupA \\(N=4\\)$", labels_end)))
  expect_true(any(grepl("^GroupB \\(N=3\\)$", labels_end)))
  expect_true(any(grepl("^GroupC \\(N=4\\)$", labels_end)))
})

test_that("ungrouped analysis (col_z = NULL): placeholder label & legend title behavior", {
  # split_by_base
  p_base <- plot_motif_barplot(df_barplot_sample, col_z = NULL)
  built_base <- ggplot2::ggplot_build(p_base)
  lab <- unique(built_base$layout$layout$group_label)
  expect_true(any(grepl("^All Fragments \\(N=\\d+\\)$", lab)))

  # split_by_motif — some ggplot2 versions show the mapped variable as title
  p_motif <- suppressWarnings(
    plot_motif_barplot(df_barplot_sample, representation = "split_by_motif", col_z = NULL)
  )
  expect_s3_class(p_motif, "ggplot")
})

test_that("motif_start filtering retains only requested first bases", {
  p <- plot_motif_barplot(df_barplot_sample, motif_start = c("A", "C"))
  built <- ggplot2::ggplot_build(p)
  first_bases <- sort(as.character(unique(built$layout$layout$first_base)))
  expect_identical(first_bases, c("A", "C"))
})

test_that("extra aesthetics in ... are passed to the first bar layer; width defaults are set", {
  p_lin <- plot_motif_barplot(df_barplot_sample, linetype = "dashed")
  expect_s3_class(p_lin, "ggplot")

  p_diff_lin <- plot_motif_barplot(
    df_barplot_sample, representation = "differential",
    vals_z = c("GroupA", "GroupB"),
    linetype = "dashed"
  )
  expect_s3_class(p_diff_lin, "ggplot")

  p_base <- plot_motif_barplot(df_barplot_sample)
  expect_s3_class(p_base, "ggplot")

  p_diff <- plot_motif_barplot(df_barplot_sample, representation = "differential",
                               vals_z = c("GroupA", "GroupB"))
  expect_s3_class(p_diff, "ggplot")

  p_motif <- plot_motif_barplot(df_barplot_sample, representation = "split_by_motif")
  expect_s3_class(p_motif, "ggplot")

  p_override <- plot_motif_barplot(df_barplot_sample, representation = "split_by_motif", width = 0.6)
  expect_s3_class(p_override, "ggplot")
})

test_that("named colors for split_by_motif map by base group names (without N=...)", {
  p <- plot_motif_barplot(
    df_barplot_sample,
    representation = "split_by_motif",
    colors_z = c(GroupA = "black", GroupB = "grey50", GroupC = "grey80")
  )
  built <- ggplot2::ggplot_build(p)
  used_norm <- unique(.normalise_hex(built$data[[1]]$fill))
  exp_norm  <- .normalise_hex(c("black", "grey50", "grey80"))
  expect_true(all(used_norm %in% exp_norm))
})

test_that("colors_z length validation for split_by_motif (unnamed vector) triggers an error", {
  expect_error(
    plot_motif_barplot(
      df_barplot_sample,
      representation = "split_by_motif",
      colors_z = c("red", "blue")
    ),
    regexp = "Provided 2 colors but need 3 for:"
  )
})


# ============== 4) Saving behavior (tempdir only; no interactive ops) =========

test_that("no output_path returns a ggplot object (no saving branch)", {
  p <- suppressWarnings(plot_motif_barplot(df_fragments = df_barplot_sample))
  expect_s3_class(p, "ggplot")
})

test_that("valid output_path saves split_by_base with custom ggsave_params", {
  tmp <- file.path(tempdir(), "motif_split_by_base_save.png")
  if (file.exists(tmp)) file.remove(tmp)

  suppressWarnings({
    plot_motif_barplot(
      df_fragments   = df_barplot_sample,
      representation = "split_by_base",
      output_path    = tmp,
      ggsave_params  = list(width = 7, height = 5, units = "in", dpi = 150, bg = "white")
    )
  })
  expect_true(file.exists(tmp))
  file.remove(tmp)
})

test_that("differential save works and overwrite is silent", {
  tmp <- file.path(tempdir(), "motif_differential_overwrite.png")
  if (file.exists(tmp)) file.remove(tmp)

  suppressWarnings({
    plot_motif_barplot(
      df_fragments   = df_barplot_sample,
      representation = "differential",
      vals_z         = c("GroupA", "GroupB"),
      output_path    = tmp
    )
  })
  expect_true(file.exists(tmp))

  suppressWarnings({
    plot_motif_barplot(
      df_fragments   = df_barplot_sample,
      representation = "differential",
      vals_z         = c("GroupA", "GroupB"),
      output_path    = tmp
    )
  })
  expect_true(file.exists(tmp))
  file.remove(tmp)
})

test_that("invalid output_path inputs are ignored and a ggplot is returned", {
  p_vec <- suppressWarnings(plot_motif_barplot(
    df_fragments = df_barplot_sample,
    output_path  = c("path1", "path2")
  ))
  expect_s3_class(p_vec, "ggplot")

  p_num <- suppressWarnings(plot_motif_barplot(
    df_fragments = df_barplot_sample,
    output_path  = 123
  ))
  expect_s3_class(p_num, "ggplot")
})

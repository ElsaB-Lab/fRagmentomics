# --- Sample Data for Testing ---

df_freq_sample <- data.frame(
    # GroupA is designed to be AC-rich, GroupB is GT-rich.
    Fragment_Bases_5p = c("ACAC", "CACA", "GTGT", "TGTG", "AC", NA, "AAAA"),
    Fragment_Bases_3p = c("CACA", "ACAC", "TGTG", "GTGT", "GT", "ACAC", "CCCC"),
    Fragment_Status_Simple = c("GroupA", "GroupA", "GroupB", "GroupB", "GroupA", "GroupB", "GroupC")
)

# --- Test File for plot_freq_barplot ---

# ============== 1. Input Validation and Error Handling Tests ==============

test_that("Function stops for invalid or inconsistent arguments", {
    # Error for invalid 'motif_type'.
    expect_error(
        plot_freq_barplot(df_freq_sample, motif_type = "invalid_type"),
        regexp = "motif_type must be one of 'Start', 'End', or 'Both'."
    )

    # Error when col_z is NULL but vals_z is not.
    expect_error(
        plot_freq_barplot(df_freq_sample, col_z = NULL, vals_z = "GroupA"),
        regexp = "If 'col_z' is NULL, 'vals_z' must also be NULL."
    )

    # Error when col_z does not exist.
    expect_error(
        plot_freq_barplot(df_freq_sample, col_z = "NonExistentColumn"),
        regexp = "Column 'NonExistentColumn' not found in the dataframe.",
        fixed = TRUE
    )
})

test_that("Function stops for issues related to data filtering and grouping", {
    # Error when filtering with vals_z results in an empty dataframe.
    expect_error(
        plot_freq_barplot(df_freq_sample, vals_z = "NonExistentGroup"),
        regexp = "No data remains after filtering. Check 'col_z' and 'vals_z'."
    )

    # Error when a grouped analysis is requested with fewer than two groups.
    expect_error(
        plot_freq_barplot(df_freq_sample, vals_z = "GroupA"),
        regexp = "Grouped analysis requires at least two groups for comparison."
    )
})


# ============== 2. Data Handling and Warning Tests ==============

test_that("Motif size is correctly adjusted with a warning", {
    # The shortest sequence in df_freq_sample is "AC" (length 2).
    # Requesting motif_size = 3 (the default) should trigger a warning.
    expect_warning(
        p <- plot_freq_barplot(df_freq_sample, vals_z = c("GroupA", "GroupB"), colors_z = "paired"),
        regexp = "Requested 'motif_size' \\(3\\) is larger than the shortest available sequence \\(2\\). Using maximum possible size: 2."
    )

    # The plot title should reflect the adjusted motif size.
    expect_match(p$labels$title, "Overall Nucleotide Frequency \\(2-mers\\)")
})

test_that("Different 'motif_type' options work correctly", {
    suppressWarnings({
        p_start <- plot_freq_barplot(df_freq_sample, motif_type = "Start", vals_z = c("GroupA", "GroupB"))
        p_end <- plot_freq_barplot(df_freq_sample, motif_type = "End", vals_z = c("GroupA", "GroupB"))
    })

    expect_s3_class(p_start, "ggplot")
    expect_s3_class(p_end, "ggplot")
})


# ============== 3. Core Functionality Tests (Grouped vs. Ungrouped) ==============

test_that("Grouped analysis (default) works and includes statistical test", {
    p <- plot_freq_barplot(df_freq_sample, vals_z = c("GroupA", "GroupB"), motif_size = 2, colors_z = "paired")
    expect_s3_class(p, "ggplot")

    # Check for essential plot layers.
    geoms <- sapply(p$layers, function(l) class(l$geom)[1])
    expect_true("GeomBar" %in% geoms)
    expect_true("GeomErrorbar" %in% geoms)

    # The plot caption should contain the result of the Chi-squared test.
    expect_match(p$labels$caption, "Global comparison \\(Chi-squared test\\), p-value")

    # Find the fill scale dynamically instead of assuming it's the first one
    fill_scale_idx <- which(sapply(p$scales$scales, function(s) "fill" %in% s$aesthetics))
    legend_labels <- p$scales$scales[[fill_scale_idx]]$labels

    expect_true(all(grepl("GroupA \\(N=\\d+\\)|GroupB \\(N=\\d+\\)", legend_labels)))
})

test_that("Ungrouped analysis (col_z = NULL) works correctly", {
    # Run the function with col_z = NULL to trigger ungrouped analysis.
    p <- plot_freq_barplot(df_freq_sample, col_z = NULL, motif_size = 2, colors_z = "paired")
    expect_s3_class(p, "ggplot")

    # The caption should NOT contain the Chi-squared test result.
    expect_false(grepl("Chi-squared test", p$labels$caption))

    # Find the fill scale dynamically
    fill_scale_idx <- which(sapply(p$scales$scales, function(s) "fill" %in% s$aesthetics))
    legend_labels <- p$scales$scales[[fill_scale_idx]]$labels

    expect_true(grepl("All Fragments", legend_labels))
})


# ============== 4. Customization and Parameter Tests ==============

test_that("Filtering with 'vals_z' works as expected", {
    # Use vals_z to select two out of the three available groups.
    suppressWarnings({
        p <- plot_freq_barplot(df_freq_sample, vals_z = c("GroupA", "GroupC"), motif_size = 2, colors_z = "paired")
    })
    expect_s3_class(p, "ggplot")

    # Check the data used for plotting to ensure only the selected groups are present.
    plotted_groups <- as.character(unique(p$data$group))
    expect_equal(sort(plotted_groups), c("GroupA", "GroupC"))
})

test_that("Coloring schemes ('colors_z') are applied correctly", {
    # Generate the full "Paired" palette as a reliable reference
    full_paired_palette <- RColorBrewer::brewer.pal(12, "Paired")

    # Run the plot function (suppress warnings from prop.test/brewer.pal)
    p_brewer <- suppressWarnings(plot_freq_barplot(df_freq_sample,
        vals_z = c("GroupA", "GroupB"),
        motif_size = 2, colors_z = "Paired"
    ))

    # Build the plot and get the actual colors used
    built_brewer <- ggplot_build(p_brewer)
    actual_brewer_colors <- unique(built_brewer$data[[1]]$fill)

    # Check that exactly 2 colors were used for the 2 groups.
    expect_length(actual_brewer_colors, 2)
    # Check that the colors used are from the correct "Paired" palette.
    expect_true(all(actual_brewer_colors %in% full_paired_palette))

    # Test using a custom vector of colors
    custom_colors <- c("gold", "purple")
    p_custom <- suppressWarnings(plot_freq_barplot(df_freq_sample,
        vals_z = c("GroupA", "GroupB"),
        motif_size = 2, colors_z = custom_colors
    ))
    built_custom <- ggplot_build(p_custom)
    actual_custom_colors <- unique(built_custom$data[[1]]$fill)
    expect_equal(sort(actual_custom_colors), sort(custom_colors))
})

test_that("Additional arguments (...) are passed to geom_bar", {
    # Pass 'alpha' as an additional argument.
    p <- plot_freq_barplot(df_freq_sample,
        vals_z = c("GroupA", "GroupB"),
        motif_size = 2, alpha = 0.5, colors_z = "paired"
    )
    # Check that the alpha aesthetic is set in the GeomBar layer.
    expect_equal(p$layers[[1]]$aes_params$alpha, 0.5)
})


# ============== 5. File Saving Tests ==============

test_that("File saving functionality works correctly under various conditions", {
  # Use a temporary directory to avoid creating files in the project
  temp_dir <- tempdir()

  # --- Test 1: Directory creation ---
  output_subdir <- file.path(temp_dir, "new_test_dir")
  # Ensure the directory does not exist initially
  if (dir.exists(output_subdir)) unlink(output_subdir, recursive = TRUE)
  expect_false(dir.exists(output_subdir))

  suppressWarnings({
    plot_freq_barplot(
      df_freq_sample,
      vals_z = c("GroupA", "GroupB"),
      output_folder = output_subdir,
      sample_id = "dir_creation_test"
    )
  })

  # Check that the directory was created...
  expect_true(dir.exists(output_subdir))
  # ...and that the file was saved inside it.
  expected_file_1 <- file.path(output_subdir, "dir_creation_test_nucleotide_frequency.png")
  expect_true(file.exists(expected_file_1))
  # Clean up the created directory
  unlink(output_subdir, recursive = TRUE)


  # --- Test 2: Correct filename generation (with and without sample_id) ---
  # With sample_id
  suppressWarnings({
    plot_freq_barplot(
      df_freq_sample,
      vals_z = c("GroupA", "GroupB"),
      output_folder = temp_dir,
      sample_id = "sample123"
    )
  })
  expected_file_2 <- file.path(temp_dir, "sample123_nucleotide_frequency.png")
  expect_true(file.exists(expected_file_2))
  file.remove(expected_file_2)

  # Without sample_id (NA)
  suppressWarnings({
    plot_freq_barplot(
      df_freq_sample,
      vals_z = c("GroupA", "GroupB"),
      output_folder = temp_dir,
      sample_id = NA
    )
  })
  expected_file_3 <- file.path(temp_dir, "nucleotide_frequency.png")
  expect_true(file.exists(expected_file_3))


  # --- Test 3: Overwrite message ---
  # Create a dummy file first
  file.create(expected_file_3)
  expect_true(file.exists(expected_file_3))

  # Expect a message warning about the overwrite
  suppressWarnings({
    expect_message(
      plot_freq_barplot(
        df_freq_sample,
        vals_z = c("GroupA", "GroupB"),
        output_folder = temp_dir,
        sample_id = NA
      ),
      regexp = "already exists and will be overwritten"
    )
  })
  file.remove(expected_file_3)
})

  # --- Test 4: Custom ggsave_params are used ---
test_that("File saving input validation works", {
  # Error if output_folder is not a single string
  expect_error(
    suppressWarnings(plot_freq_barplot(df_freq_sample, output_folder = c("path1", "path2"))),
    regexp = "'output_folder' must be a single character string."
  )

  # Error if sample_id is not a single string (and not NA)
  expect_error(
    suppressWarnings(plot_freq_barplot(df_freq_sample, output_folder = tempdir(), sample_id = 123)),
    regexp = "'sample_id' must be a single character string."
  )
})

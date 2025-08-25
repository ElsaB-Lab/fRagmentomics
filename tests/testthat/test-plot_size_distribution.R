# --- Test File for plot_size_distribution ---

# Sample data
df_sample <- data.frame(
    Fragment_Size = c(rnorm(100, 150, 20), rnorm(100, 320, 30), rnorm(50, 480, 40), NA),
    Fragment_Status_Simple = c(rep("Mono", 100), rep("Di", 100), rep("Tri", 50), "Mono"),
    Other_Info = sample(letters, 251, replace = TRUE)
)

# Add an NA to the grouping column as well.
df_sample$Fragment_Status_Simple[1] <- NA

# ============== 1. Input Validation Tests ==============
# These tests ensure that the function correctly stops and provides
# informative error messages for invalid inputs.

test_that("Function throws errors for invalid or missing columns", {
    # Test for error when the size column does not exist.
    expect_error(
        plot_size_distribution(df_sample, size_col = "NonExistentSizeCol"),
        regexp = "Size column NonExistentSizeCol not found in the dataframe."
    )

    # Test for error when the grouping column does not exist.
    expect_error(
        plot_size_distribution(df_sample, col_z = "NonExistentGroupCol"),
        regexp = "Column NonExistentGroupCol not found in the dataframe."
    )

    # Test for error when the size column is not numeric.
    df_bad_type <- df_sample
    df_bad_type$Fragment_Size <- as.character(df_bad_type$Fragment_Size)
    expect_error(
        plot_size_distribution(df_bad_type, size_col = "Fragment_Size"),
        regexp = "Size column Fragment_Size must be numeric."
    )
})

test_that("Function throws errors for logical inconsistencies in arguments", {
    # Test for error when col_z is NULL but vals_z is provided.
    expect_error(
        plot_size_distribution(df_sample, col_z = NULL, vals_z = c("Mono", "Di")),
        regexp = "If 'col_z' is NULL, 'vals_z' must also be NULL."
    )

    # Test for error when both histogram and density plots are disabled.
    expect_error(
        plot_size_distribution(df_sample, show_histogram = FALSE, show_density = FALSE),
        regexp = "At least one of 'show_histogram' or 'show_density' must be TRUE."
    )
})

test_that("Function throws errors for empty data after filtering", {
    # Test for error when filtering with vals_z results in an empty dataframe.
    expect_error(
        plot_size_distribution(df_sample, vals_z = "NonExistentValue"),
        regexp = "No data remains after filtering. Check 'col_z' and 'vals_z'."
    )
})


# ============== 2. Core Functionality and Output Tests ==============
# These tests check that the function produces the correct output (a ggplot object)
# under various standard configurations.

test_that("Function returns a ggplot object with default settings", {
    # Generate a plot with default arguments.
    p <- plot_size_distribution(df_sample)
    # The output should be a ggplot object.
    expect_s3_class(p, "ggplot")
    # The plot should contain a density layer ('GeomDensity') by default.
    expect_true("GeomDensity" %in% sapply(p$layers, function(l) class(l$geom)[1]))
})

test_that("Ungrouped analysis (col_z = NULL) works correctly", {
    # Generate a plot without grouping.
    p <- plot_size_distribution(df_sample, col_z = NULL)
    expect_s3_class(p, "ggplot")

    # Check if the legend label is correct for the placeholder group.
    # We build the plot to inspect its components.
    plot_data <- ggplot_build(p)$data[[1]]
    # There should only be one group.
    expect_equal(length(unique(plot_data$group)), 1)
    # The label should contain the placeholder text "All Fragments" and the count.
    expect_true(grepl("All Fragments \\(N=\\d+\\)", p$data$placeholder_group[1]))
})

test_that("Filtering by `vals_z` works as expected", {
    # Filter to show only "Mono" and "Di" groups.
    groups_to_show <- c("Mono", "Di")
    p <- plot_size_distribution(df_sample, vals_z = groups_to_show)
    expect_s3_class(p, "ggplot")

    # Build the plot to check the labels in the legend.
    built_p <- ggplot_build(p)
    # Extract the labels from the color scale.
    scale_idx <- which(sapply(built_p$plot$scales$scales, function(s) "colour" %in% s$aesthetics))
    legend_labels <- built_p$plot$scales$scales[[scale_idx]]$get_labels()

    # Check that only the specified groups are present in the legend.
    expect_equal(length(legend_labels), 2)
    expect_true(all(grepl("Mono|Di", legend_labels)))
    # Ensure the "Tri" group is not present.
    expect_false(any(grepl("Tri", legend_labels)))
})

# ============== 3. Specific Parameterization Tests ==============
# These tests verify that individual parameters correctly modify the plot output.

test_that("Plot layers are correctly added or removed", {
    # Test with histogram only.
    p_hist <- plot_size_distribution(df_sample, show_histogram = TRUE, show_density = FALSE)
    geoms_hist <- sapply(p_hist$layers, function(l) class(l$geom)[1])
    expect_true("GeomBar" %in% geoms_hist) # geom_histogram creates GeomBar
    expect_false("GeomDensity" %in% geoms_hist)
    # The y-axis label should be "Density".
    expect_equal(p_hist$labels$y, "Density")

    # Test with both histogram and density.
    p_both <- plot_size_distribution(df_sample, show_histogram = TRUE, show_density = TRUE)
    geoms_both <- sapply(p_both$layers, function(l) class(l$geom)[1])
    expect_true("GeomBar" %in% geoms_both)
    expect_true("GeomDensity" %in% geoms_both)

    # Test turning off nucleosome peaks.
    p_no_peaks <- plot_size_distribution(df_sample, show_nuc_peaks = FALSE)
    geoms_no_peaks <- sapply(p_no_peaks$layers, function(l) class(l$geom)[1])
    # There should be no vertical line layer ('GeomVline').
    expect_false("GeomVline" %in% geoms_no_peaks)
})

test_that("Custom arguments for geoms are passed correctly", {
    # Test custom density arguments (e.g., changing line width).
    p_density <- plot_size_distribution(df_sample, density_args = list(linewidth = 2))
    # The linewidth parameter should be set to 2 in the GeomDensity layer.
    expect_equal(p_density$layers[[1]]$aes_params$linewidth, 2)

    # Test custom histogram arguments (e.g., changing alpha).
    p_hist <- plot_size_distribution(df_sample, show_histogram = TRUE, show_density = FALSE, histo_args = list(alpha = 0.5))
    # The alpha parameter should be set to 0.5 in the GeomBar layer.
    expect_equal(p_hist$layers[[1]]$aes_params$alpha, 0.5)
})

test_that("Coloring schemes are applied correctly", {
    # Test using a vector of custom colors.
    custom_colors <- c("red", "green", "blue")
    # We must specify vals_z to ensure a consistent order for colors.
    p_custom_cols <- plot_size_distribution(df_sample, vals_z = c("Mono", "Di", "Tri"), colors_z = custom_colors)
    # Check if the manual color scale was applied with the correct colors.
    scale_values <- ggplot_build(p_custom_cols)$plot$scales$scales[[1]]$palette(3)
    expect_equal(tolower(scale_values), custom_colors)

    # Test using an RColorBrewer palette name.
    p_brewer <- plot_size_distribution(df_sample, vals_z = c("Mono", "Di", "Tri"), colors_z = "Set1")
    # Get the expected colors from RColorBrewer.
    expected_brewer_colors <- RColorBrewer::brewer.pal(3, "Set1")
    scale_values_brewer <- ggplot_build(p_brewer)$plot$scales$scales[[1]]$palette(3)
    expect_equal(scale_values_brewer, expected_brewer_colors)
})

test_that("Axis limits are set correctly", {
    # Define custom x-axis limits.
    custom_limits <- c(100, 500)
    p <- plot_size_distribution(df_sample, x_limits = custom_limits)
    # Check if the coordinate system limits match the provided limits.
    coord_limits <- ggplot_build(p)$layout$coord$limits$x
    expect_equal(coord_limits, custom_limits)
})

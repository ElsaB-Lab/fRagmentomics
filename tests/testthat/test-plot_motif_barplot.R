# --- Sample Data for Testing ---

df_barplot_sample <- data.frame(
    Fragment_Bases_5p = c(
        "ACG", "ACT", "AGA", # GroupA
        "AGT", "CAT", "CCT", # GroupB
        "GAT", "GGT", "TGA", # GroupC
        "NNN", "ANC", NA # Invalid motifs spread across groups
    ),
    Fragment_Bases_3p = c(
        "TGA", "TCA", "TCT", # GroupA
        "TGT", "GTA", "AGG", # GroupB
        "ATC", "ACC", "TCA", # GroupC
        "TCA", "ANT", "CCT" # Invalid motifs
    ),
    Fragment_Status_Simple = c(
        "GroupA", "GroupA", "GroupA",
        "GroupB", "GroupB", "GroupB",
        "GroupC", "GroupC", "GroupC",
        "GroupA", "GroupB", "GroupC"
    )
)

# --- Test File for plot_motif_barplot ---

# ============== 1. Input Validation and Error Handling Tests ==============

test_that("Function stops for invalid arguments and inconsistent parameters", {
    # Error for invalid 'representation' string.
    expect_error(
        plot_motif_barplot(df_barplot_sample, representation = "invalid_type"),
        regexp = "'representation' must be one of: differential, split_by_base, split_by_motif"
    )

    # Error when col_z is NULL but vals_z is not.
    expect_error(
        plot_motif_barplot(df_barplot_sample, col_z = NULL, vals_z = "GroupA"),
        regexp = "If 'col_z' is NULL, 'vals_z' must also be NULL."
    )

    # Error when col_z does not exist in the dataframe.
    expect_error(
        plot_motif_barplot(df_barplot_sample, col_z = "NonExistentColumn"),
        regexp = "Column NonExistentColumn not found in the dataframe."
    )
})

test_that("Function stops for representation-specific requirement violations", {
    # Differential analysis requires a grouping column.
    expect_error(
        plot_motif_barplot(df_barplot_sample, representation = "differential", col_z = NULL),
        regexp = "Differential analysis requires a grouping column 'col_z'."
    )

    # Differential analysis requires exactly two groups.
    expect_error(
        plot_motif_barplot(df_barplot_sample, representation = "differential", vals_z = "GroupA"),
        regexp = "Differential analysis requires exactly two values in 'vals_z'."
    )
    expect_error(
        plot_motif_barplot(df_barplot_sample, representation = "differential", vals_z = c("GroupA", "GroupB", "GroupC")),
        regexp = "Differential analysis requires exactly two values in 'vals_z'."
    )
})

test_that("Function stops when data becomes empty after filtering", {
    # Error when vals_z filters out all data.
    expect_error(
        plot_motif_barplot(df_barplot_sample, vals_z = "NonExistentGroup"),
        regexp = "No data remains after filtering. Check 'col_z' and 'vals_z'."
    )

    # Error when no valid 3-base motifs are found in the data.
    df_no_valid_motifs <- data.frame(
        Fragment_Bases_5p = c("NN", "AT", "NNA"),
        Fragment_Bases_3p = c("N", "G", "C"),
        Fragment_Status_Simple = c("A", "A", "A")
    )
    expect_error(
        plot_motif_barplot(df_no_valid_motifs),
        regexp = "No valid 3-base motifs found in the data."
    )
})


# ============== 2. Core Functionality & Representation Tests ==============

test_that("Default representation 'split_by_base' works correctly", {
    # Test with default settings.
    p <- plot_motif_barplot(df_barplot_sample)
    expect_s3_class(p, "ggplot")
    # Check for the correct faceting function from ggh4x.
    expect_s3_class(p$facet, "FacetNested")

    # Test with an RColorBrewer palette name.
    p_brewer <- plot_motif_barplot(df_barplot_sample, colors_z = "Set1")
    # Build plot to inspect scales.
    built_p <- ggplot_build(p_brewer)
    expected_colors <- RColorBrewer::brewer.pal(4, "Set1")
    # The fill scale should use the specified Brewer palette.
    expect_setequal(unique(built_p$data[[1]]$fill), expected_colors)
})

test_that("Representation 'split_by_motif' works correctly", {
    p <- plot_motif_barplot(df_barplot_sample, representation = "split_by_motif")
    expect_s3_class(p, "ggplot")
    # Check that faceting is by FacetWrap.
    expect_s3_class(p$facet, "FacetWrap")
    # Check that x-axis text is rotated.
    expect_equal(p$theme$axis.text.x$angle, 90)
    # Check that plot title is correct for this representation.
    expect_equal(p$labels$title, "Motif Proportions by Group")
})

test_that("Representation 'differential' works correctly", {
    p <- plot_motif_barplot(df_barplot_sample, representation = "differential", vals_z = c("GroupA", "GroupB"))
    expect_s3_class(p, "ggplot")
    # Check y-axis label for Log2 Fold Change.
    expect_match(p$labels$y, "Log2 Fold Change")
    # Check subtitle for correct comparison groups.
    expect_equal(p$labels$subtitle, "Comparison: GroupA vs GroupB")
    # Check for correct faceting.
    expect_s3_class(p$facet, "FacetNested")

    # Ensure the data contains log2FC and sign columns.
    expect_true(all(c("log2FC", "sign") %in% names(p$data)))
})


# ============== 3. Parameter and Data Handling Tests ==============

test_that("Analysis of different motif types ('Start', 'End', 'Both') works", {
    # Test 'Start' only.
    p_start <- plot_motif_barplot(df_barplot_sample, motif_type = "Start")
    expect_s3_class(p_start, "ggplot")
    # Check that the number of motifs is as expected from the 5p column.
    # 3 in A, 3 in B, 3 in C are valid. Total N=9.
    group_labels <- ggplot_build(p_start)$layout$layout$group_label
    expect_true(all(grepl("N=3", group_labels)))

    # Test 'End' only.
    p_end <- plot_motif_barplot(df_barplot_sample, motif_type = "End")
    expect_s3_class(p_end, "ggplot")
    # 4 in A, 4 in B, 3 in C are valid. Total N=10.
    group_labels_end <- ggplot_build(p_end)$layout$layout$group_label
    expect_true(grepl("N=4", group_labels_end[1])) # Group A
    expect_true(grepl("N=4", group_labels_end[2])) # Group B
    expect_true(grepl("N=3", group_labels_end[3])) # Group C
})

test_that("Ungrouped analysis (col_z = NULL) works correctly", {
    # Test ungrouped 'split_by_base'.
    p_base <- plot_motif_barplot(df_barplot_sample, col_z = NULL)
    expect_s3_class(p_base, "ggplot")
    # Check that the facet label shows the placeholder text.
    facet_label <- ggplot_build(p_base)$layout$layout$group_label
    expect_true(grepl("All Fragments", unique(facet_label)))

    # Test ungrouped 'split_by_motif'.
    p_motif <- plot_motif_barplot(df_barplot_sample, representation = "split_by_motif", col_z = NULL)
    expect_s3_class(p_motif, "ggplot")
    # Check that the legend title is the placeholder column name.
    expect_equal(p_motif$labels$fill, "placeholder_group")
})

test_that("Filtering by 'motif_start' works", {
    # Filter motifs to only those starting with 'A' or 'C'.
    p <- plot_motif_barplot(df_barplot_sample, motif_start = c("A", "C"))
    expect_s3_class(p, "ggplot")
    # Build the plot to check the facets.
    built_p <- ggplot_build(p)
    first_base_facets <- unique(built_p$layout$layout$first_base)
    # Only 'A' and 'C' should be present as first_base facets.
    expect_equal(sort(as.character(first_base_facets)), c("A", "C"))
})

test_that("Additional arguments (...) are passed to geom_bar", {
    # Pass 'linetype' as an additional argument.
    p <- plot_motif_barplot(df_barplot_sample, linetype = "dashed")
    # Check that the linetype aesthetic is set in the geom layer.
    expect_equal(p$layers[[1]]$aes_params$linetype, "dashed")

    # Test with the 'differential' plot as well.
    p_diff <- plot_motif_barplot(df_barplot_sample,
        representation = "differential",
        vals_z = c("GroupA", "GroupB"), linetype = "dashed"
    )
    expect_equal(p_diff$layers[[1]]$aes_params$linetype, "dashed")
})

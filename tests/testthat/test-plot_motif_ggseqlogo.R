# --- Sample Data for Testing ---
df_motif_sample <- data.frame(
    # 5' end motifs of varying lengths, including NA and 'N'.
    Fragment_Bases_5p = c(
        "ACGT", "ACGG", "AGTA", "CCTT", # Group 'Short'
        "TGCAT", "TGAAC", "TGCAC", NA, # Group 'Long'
        "GATTACA", "GATTANA", # Group 'Mixed' (with 'N')
        "GC" # Too short for motif_size=3
    ),
    # 3' end motifs.
    Fragment_Bases_3p = c(
        "TGCA", "CCGG", "TACT", "AAGG", # Group 'Short'
        "ATGCA", "GTTCA", "GTGCA", "A", # Group 'Long' (one is too short)
        "TACAATG", "TANAATG", # Group 'Mixed' (with 'N')
        "CG" # Too short
    ),
    # A grouping variable, including an NA to test filtering.
    Fragment_Status_Simple = c(
        rep("Short", 4),
        rep("Long", 4),
        rep("Mixed", 2),
        "Short"
    )
)

# --- Test File for plot_qqseqlogo_meme ---

# ============== 1. Input Validation Tests ==============
# These tests ensure the function stops with informative errors for invalid inputs.

test_that("Function throws errors for invalid arguments", {
    # Error when col_z is NULL but vals_z is not.
    expect_error(
        plot_qqseqlogo_meme(df_motif_sample, col_z = NULL, vals_z = "Short"),
        regexp = "If 'col_z' is NULL, 'vals_z' must also be NULL."
    )

    # Error when the grouping column doesn't exist.
    expect_error(
        plot_qqseqlogo_meme(df_motif_sample, col_z = "NonExistentCol"),
        regexp = "Column NonExistentCol not found in the dataframe."
    )

    # Error for an invalid motif_type.
    expect_error(
        plot_qqseqlogo_meme(df_motif_sample, motif_type = "InvalidType"),
        regexp = "motif_type must be one of 'Start', 'End', or 'Both'."
    )
})

# ============== 2. Data Filtering and Grouping Tests ==============

test_that("Data is correctly filtered and grouped", {
    # --- Test 1: Ungrouped analysis (col_z = NULL) ---
    suppressWarnings({
        p_ungrouped <- plot_qqseqlogo_meme(df_motif_sample, col_z = NULL, motif_size = 2)
    })
    
    # Verify the object is a ggplot
    expect_s3_class(p_ungrouped, "ggplot")
    
    # Get the facet label from the plot's layout definition using the correct column name 'seq_group'.
    plot_layout <- ggplot_build(p_ungrouped)$layout$layout
    facet_label <- as.character(plot_layout$seq_group)
    
    # Check that the label is correct.
    expect_equal(facet_label, "All Fragments (N=11)")


    # --- Test 2: Grouped analysis with two groups ---
    suppressWarnings({
        p_grouped <- plot_qqseqlogo_meme(
            df_motif_sample,
            vals_z = c("Short", "Long"),
            motif_size = 3
        )
    })
    
    # Verify the object is a ggplot
    expect_s3_class(p_grouped, "ggplot")
    
    # Use the same layout method for the multi-facet plot.
    plot_layout_grouped <- ggplot_build(p_grouped)$layout$layout
    facet_labels <- as.character(plot_layout_grouped$seq_group)
    
    # We expect two facet labels.
    expect_equal(length(facet_labels), 2)
    
    # Check if the labels are correct (including counts).
    expected_labels <- c("Short (N=5)", "Long (N=4)")
    expect_true(all(sort(facet_labels) == sort(expected_labels)))
})

# ============== 3. Motif Extraction and Handling Tests ==============

test_that("Motif extraction works for 'Start', 'End', and 'Both'", {
    # Test 'Start' motif type.
    p_start <- plot_qqseqlogo_meme(df_motif_sample, motif_type = "Start", motif_size = 2, vals_z = "Short")
    expect_s3_class(p_start, "ggplot")

    # Build the plot to inspect its data
    p_start_build <- ggplot_build(p_start)
    # The motif length is the ROUNDED max x-position in the plot data
    motif_len <- round(max(p_start_build$data[[1]]$x))
    expect_equal(motif_len, 2)

    # Check that no vertical line is present.
    geoms <- sapply(p_start$layers, function(l) class(l$geom)[1])
    expect_false("GeomVline" %in% geoms)

    # ---

    # Test 'End' motif type.
    p_end <- plot_qqseqlogo_meme(df_motif_sample, motif_type = "End", motif_size = 2, vals_z = "Short")
    expect_s3_class(p_end, "ggplot")
    
    p_end_build <- ggplot_build(p_end)
    motif_len <- round(max(p_end_build$data[[1]]$x))
    expect_equal(motif_len, 2)

    # ---

    # Test 'Both' motif type.
    suppressWarnings({
        p_both <- plot_qqseqlogo_meme(df_motif_sample, motif_size = 3, vals_z = "Short")
    })
    
    expect_s3_class(p_both, "ggplot")

    p_both_build <- ggplot_build(p_both)
    # Motif length for 'Both' is motif_size (5') + 1 (hyphen) + motif_size (3')
    motif_len <- round(max(p_both_build$data[[1]]$x))
    expect_equal(motif_len, 3 + 1 + 3)

    # Check that a vertical line is present for 'Both'.
    geoms <- sapply(p_both$layers, function(l) class(l$geom)[1])
    expect_true("GeomVline" %in% geoms)
})

test_that("Function issues warnings correctly", {
    # Warning for truncating motif_size when sequences are too short.
    # Group 'Short' has one sequence "GC" of length 2. motif_size=3 should trigger a warning.
    expect_warning(
        plot_qqseqlogo_meme(df_motif_sample, motif_type = "Start", vals_z = "Short", motif_size = 3),
        regexp = "5' motif_size \\(3\\) truncated to 2"
    )
    expect_warning(
        plot_qqseqlogo_meme(df_motif_sample, motif_type = "End", vals_z = "Long", motif_size = 4),
        regexp = "3' motif_size \\(4\\) truncated to 1"
    )

    # Warning for removing motifs containing 'N'.
    # The 'Mixed' group contains 'N', so it should trigger this warning.
    expect_warning(
        plot_qqseqlogo_meme(df_motif_sample, vals_z = "Mixed", motif_size = 5),
        regexp = "motifs containing 'N' were found and removed"
    )
})

test_that("Function handles cases with no data to plot", {
    # Create a dataframe where all motifs will be filtered out.
    df_only_n <- data.frame(
        Fragment_Bases_5p = "NNN",
        Fragment_Bases_3p = "NNN",
        Fragment_Status_Simple = "A"
    )

    # Expect a message and a blank ggplot object.
    # Suppress the warning about 'N' removal, as it's expected behavior.
    expect_message(
        suppressWarnings({
            p_empty <- plot_qqseqlogo_meme(df_only_n)
        }),
        regexp = "No data available to plot after filtering."
    )
    
    # The returned object should still be a ggplot.
    expect_s3_class(p_empty, "ggplot")
    # It should have a specific title indicating no data.
    expect_equal(p_empty$labels$title, "No data to display")
})

# ============== 4. Color Scheme Tests ==============

test_that("Color schemes are applied correctly", {
    # Test using a named vector of colors.
    custom_colors <- c("A" = "blue", "C" = "red", "G" = "green", "T" = "purple", "-" = "grey")
    p_custom <- suppressWarnings(plot_qqseqlogo_meme(df_motif_sample, motif_size = 1, colors_z = custom_colors))
    
    p_custom_build <- ggplot_build(p_custom)
    used_colors <- unique(p_custom_build$data[[1]]$fill)
    
    # Expect only 4 colors, as ggseqlogo ignores the "-" and its "grey" color.
    expected_colors <- c("blue", "red", "green", "purple")
    
    expect_equal(sort(used_colors), sort(expected_colors))
    
    # ---

    # Test using an RColorBrewer palette name.
    p_brewer <- suppressWarnings(plot_qqseqlogo_meme(df_motif_sample, motif_size = 1, colors_z = "Dark2"))
    p_brewer_build <- ggplot_build(p_brewer)
    
    used_colors_brewer <- unique(p_brewer_build$data[[1]]$fill)
    
    # FIX: The function gets 5 colors (for A,C,G,T,-) but ggseqlogo only uses the last 4 (for A,C,G,T).
    # We replicate that logic here for a correct comparison.
    full_palette <- RColorBrewer::brewer.pal(5, "Dark2")
    expected_cols_brewer <- full_palette[2:5] # Colors for A, C, G, T
    
    expect_equal(sort(used_colors_brewer), sort(expected_cols_brewer))
    
    # ---

    # Error when Brewer palette doesn't have enough colors.
    ten_chars <- strsplit("ACGT-XYZPQ", "")[[1]]
    df_extended_alphabet <- data.frame(
        Fragment_Bases_5p = ten_chars,
        Fragment_Bases_3p = ten_chars,
        Fragment_Status_Simple = "A"
    )
    expect_error(
        plot_qqseqlogo_meme(df_extended_alphabet, motif_size = 1, colors_z = "Set1"), 
        regexp = "Palette 'Set1' only has 9 colors, but 10 are needed."
    )
    
    # ---

    # Error when named vector is missing a nucleotide definition.
    missing_colors <- c("A" = "blue", "C" = "red") 
    expect_error(
        suppressWarnings(plot_qqseqlogo_meme(df_motif_sample, motif_size = 1, colors_z = missing_colors)),
        regexp = "Custom colors provided, but missing definitions for nucleotide\\(s\\): -, G, T"
    )
    
    # ---

    # Error for an invalid colors_z argument type.
    expect_error(
        suppressWarnings(plot_qqseqlogo_meme(df_motif_sample, colors_z = 123)),
        regexp = "'colors_z' must be NULL, a valid RColorBrewer palette name, or a named character vector."
    )
})

# ============== 5. Plot Theming and Annotation Tests ==============

# These tests verify the final visual components of the plot.
test_that("Custom axis labels are generated correctly for each motif_type", {
    # Test labels for 'Start'.
    # Suppress warnings because motif_size=3 is larger than some sequences.
    p_start <- suppressWarnings({
        plot_qqseqlogo_meme(df_motif_sample, motif_type = "Start", motif_size = 3, vals_z = "Short")
    })
    # The geom_text layer (Layer 2) contains the custom labels.
    label_data_start <- ggplot_build(p_start)$data[[2]]
    expect_equal(label_data_start$label, as.character(1:3))

    # Test labels for 'End'.
    p_end <- suppressWarnings({
        plot_qqseqlogo_meme(df_motif_sample, motif_type = "End", motif_size = 3, vals_z = "Short")
    })
    label_data_end <- ggplot_build(p_end)$data[[2]]
    expect_equal(label_data_end$label, as.character((-3):-1))

    # Test labels for 'Both'.
    p_both <- plot_qqseqlogo_meme(df_motif_sample, motif_type = "Both", motif_size = 2, vals_z = "Short")
    # The labels are in Layer 2 (geom_text), not Layer 3 (geom_vline).
    label_data_both <- ggplot_build(p_both)$data[[2]]
    expected_labels <- as.character(c(1:2, "", (-2):-1))
    expect_equal(label_data_both$label, expected_labels)
})

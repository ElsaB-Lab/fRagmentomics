#' Plot 3-base motif proportions with various representations
#'
#' @description Creates a bar plot to visualize the proportion of 3-base motifs at fragment ends. Supports grouped analysis and
#' three different visual representations: hierarchical faceting by base, log2 fold change, or side-by-side motifs.
#'
#' @param df_fragments The input dataframe containing fragment sequence data.
#' @param end_motif_5p Character string. Column name for 5' end sequences.
#' @param end_motif_3p Character string. Column name for 3' end sequences.
#' @param motif_type Character string. Which ends to analyze: 'Start', 'End', or 'Both'.
#' @param motif_start Optional character vector ('A','C','G','T') to filter motifs by their starting base.
#' @param col_z Character string. Column name for grouping. If NULL, no grouping is applied.
#' @param vals_z A character vector of group names from 'col_z' to include. If NULL, all unique groups in 'col_z' are used.
#' @param representation Character string. The type of plot to generate.
#'   - '"split_by_base"' (default): A proportion plot with hierarchical axes, splitting motifs by each base position into facets.
#'   - '"differential"': A log2 fold change plot comparing two groups.
#'   - '"split_by_motif"': A proportion plot with motifs on the x-axis, with bars for different groups placed side-by-side.
#' @param ... Additional arguments passed on to 'ggplot2::geom_bar()'.
#' @param colors_z For the "split_by_base" plot, a character vector of 4 colors for
#'   A, C, G, T, or a single string naming an RColorBrewer palette. For other plots,
#'   a suitable palette is chosen automatically.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr %>% filter select mutate all_of group_by summarise ungroup n count bind_rows left_join pull arrange
#' @importFrom stringr str_sub str_length str_detect
#' @importFrom tidyr drop_na pivot_wider
#' @importFrom ggh4x facet_nested
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom stats setNames
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' ## --- Create a dataset for demonstration ---
#' # Set a seed for reproducibility
#' set.seed(42)
#'
#' # Helper function to generate random DNA sequences with a bias
#' generate_biased_dna <- function(n_seq, len, prob) {
#'     bases <- c("A", "C", "G", "T")
#'     replicate(n_seq, paste(sample(bases, len, replace = TRUE, prob = prob), collapse = ""))
#' }
#'
#' # Create 50 "MUT" fragments with a high proportion of motifs starting with 'C'
#' df_mut <- data.frame(
#'     Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'     Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'     Fragment_Status_Simple = "MUT"
#' )
#'
#' # Create 50 "WT" fragments with a high proportion of motifs starting with 'G'
#' df_wt <- data.frame(
#'     Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
#'     Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
#'     Fragment_Status_Simple = "WT"
#' )
#'
#' # Combine into a single dataframe
#' example_df <- rbind(df_mut, df_wt)
#'
#' ## --- Function Calls for Each Representation ---
#'
#' # 1. Hierarchical Plot (representation = "split_by_base")
#' # This is the default. It creates nested facets for each base position.
#' p1 <- plot_motif_barplot(
#'     df_fragments = example_df,
#'     representation = "split_by_base"
#' )
#' print(p1)
#'
#' # You can also filter this plot to show only motifs starting with certain bases.
#' p1_filtered <- plot_motif_barplot(
#'     df_fragments = example_df,
#'     representation = "split_by_base",
#'     motif_start = c("C", "G")
#' )
#' print(p1_filtered)
#'
#' # 2. Differential Plot (representation = "differential")
#' # This shows the log2 fold change in motif proportions between two groups.
#' # It requires exactly two groups specified in `vals_z`.
#' p2 <- plot_motif_barplot(
#'     df_fragments = example_df,
#'     representation = "differential",
#'     vals_z = c("MUT", "WT")
#' )
#' print(p2)
#'
#' # 3. Side-by-side Motif Plot (representation = "split_by_motif")
#' # This creates a more traditional bar plot with motifs on the x-axis and
#' # bars for each group shown side-by-side.
#' p3 <- plot_motif_barplot(
#'     df_fragments = example_df,
#'     representation = "split_by_motif"
#' )
#' print(p3)
plot_motif_barplot <- function(df_fragments,
                               end_motif_5p = "Fragment_Bases_5p",
                               end_motif_3p = "Fragment_Bases_3p",
                               motif_type = "Both",
                               motif_start = NULL,
                               col_z = "Fragment_Status_Simple",
                               vals_z = NULL,
                               representation = "split_by_base",
                               ...,
                               colors_z = c("#FD96A9", "#E88B00", "#0D539E", "#6CAE75")) {
    # --- 1. Input Validation and Setup ---
    motif_size <- 3 # This plot is specifically designed for 3-base motifs.

    # Validate representation argument
    valid_representations <- c("differential", "split_by_base", "split_by_motif")
    if (!representation %in% valid_representations) {
        stop(sprintf(
            "'representation' must be one of: %s",
            paste(valid_representations, collapse = ", ")
        ))
    }
    # Store original col_z for later use in labels
    original_col_z <- col_z

    # Enforce consistent logic for grouping parameters
    if (is.null(col_z) && !is.null(vals_z)) stop("If 'col_z' is NULL, 'vals_z' must also be NULL.")
    if (!is.null(col_z) && !col_z %in% names(df_fragments)) {
        stop(sprintf("Column '%s' not found in the dataframe.", col_z))
    }
    if (representation == "differential" && is.null(col_z)) stop("Differential analysis requires a grouping column 'col_z'.")

    # Handle the case for no grouping (col_z is NULL)
    is_grouped_analysis <- !is.null(col_z)
    if (!is_grouped_analysis) {
        df_fragments$placeholder_group <- "All Fragments"
        col_z <- "placeholder_group"
    }

    # Filter the dataframe based on specified groups
    df_filtered <- df_fragments %>% filter(!is.na(.data[[col_z]]))
    if (is_grouped_analysis && !is.null(vals_z)) {
        df_filtered <- df_filtered %>% filter(.data[[col_z]] %in% vals_z)
    }
    if (nrow(df_filtered) == 0) stop("No data remains after filtering. Check 'col_z' and 'vals_z'.")

    # Determine final group levels for plotting
    if (is.null(vals_z)) {
        vals_z <- sort(unique(df_filtered[[col_z]]))
    }

    # Validate group counts for the requested analysis type
    if (representation == "differential" && length(vals_z) != 2) {
        stop("Differential analysis requires exactly two values in 'vals_z'.")
    }

    # --- 2. Data Preparation & Motif Aggregation ---

    # Extract motifs based on the 'motif_type' parameter
    motifs_list <- list()
    if (motif_type %in% c("Start", "Both")) {
        motifs_list$start <- df_filtered %>%
            select(group = all_of(col_z), motif = all_of(end_motif_5p))
    }
    if (motif_type %in% c("End", "Both")) {
        motifs_list$end <- df_filtered %>%
            select(group = all_of(col_z), motif = all_of(end_motif_3p))
    }

    # Combine, validate, and process motifs
    analysis_df <- bind_rows(motifs_list) %>%
        drop_na(motif) %>%
        mutate(motif = str_sub(motif, 1, motif_size)) %>%
        filter(str_detect(motif, paste0("^[ACGT]{", motif_size, "}$")))

    if (nrow(analysis_df) == 0) stop("No valid 3-base motifs found in the data.")

    # Calculate motif proportions per group
    proportions_df <- analysis_df %>%
        group_by(group) %>%
        count(motif, name = "count") %>%
        mutate(proportion = count / sum(count)) %>%
        ungroup()

    # --- 3. Data Transformation for Plotting ---
    base_levels <- c("A", "C", "G", "T")
    plot_data <- proportions_df %>%
        mutate(
            first_base = factor(str_sub(motif, 1, 1), levels = base_levels),
            second_base = factor(str_sub(motif, 2, 2), levels = base_levels),
            third_base = factor(str_sub(motif, 3, 3), levels = base_levels),
            # Ensure motifs are sorted lexicographically for the x-axis
            motif = factor(motif, levels = sort(unique(motif)))
        )

    if (!is.null(motif_start)) {
        plot_data <- plot_data %>% filter(first_base %in% motif_start)
    }

    # --- 4. Plot Generation (based on 'representation') ---
    plot_xlab <- "" # Initialize axis label

    if (representation == "differential") {
        # --- 4a. Differential Plot (Log2 Fold Change) ---
        pseudocount <- 1e-9
        diff_data <- plot_data %>%
            select(motif, group, proportion, first_base, second_base, third_base) %>%
            pivot_wider(names_from = group, values_from = proportion, values_fill = 0) %>%
            mutate(
                log2FC = log2((.data[[vals_z[1]]] + pseudocount) / (.data[[vals_z[2]]] + pseudocount)),
                sign = ifelse(log2FC >= 0, "Positive", "Negative")
            )

        final_plot <- ggplot(diff_data, aes(x = third_base, y = log2FC, fill = sign)) +
            geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2, alpha = 0.8, ...) +
            geom_hline(yintercept = 0, color = "darkgrey", linewidth = 1) +
            ggh4x::facet_nested(~ first_base + second_base, scales = "free_x") +
            scale_fill_manual(name = "Direction", values = c("Positive" = "#66C2A5FF", "Negative" = "#E78AC3FF"))

        plot_xlab <- "Third Base"
    } else {
        # --- 4b. Common setup for proportion plots ---
        total_counts <- analysis_df %>%
            count(group, name = "total_n")

        plot_data <- plot_data %>%
            left_join(total_counts, by = "group") %>%
            mutate(group_label = paste0(group, " (N=", total_n, ")"))

        ordered_levels <- plot_data %>%
            mutate(group = factor(group, levels = vals_z)) %>%
            arrange(group) %>%
            pull(group_label) %>%
            unique()

        plot_data$group_label <- factor(plot_data$group_label, levels = ordered_levels)

        if (representation == "split_by_base") {
            # --- Proportion Plot with hierarchical split ---
            if (length(colors_z) == 1 && colors_z %in% row.names(RColorBrewer::brewer.pal.info)) {
                colors_z <- RColorBrewer::brewer.pal(n = 4, name = colors_z)
            }

            final_plot <- ggplot(plot_data, aes(x = third_base, y = proportion, fill = second_base)) +
                geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2, alpha = 0.8, ...) +
                ggh4x::facet_nested(~ group_label + first_base + second_base, scales = "free_x") +
                scale_fill_manual(values = setNames(colors_z, base_levels), name = "2nd Base")

            plot_xlab <- "Third Base"
        } else if (representation == "split_by_motif") {
            # --- Proportion Plot with motifs on X-axis ---
            final_plot <- ggplot(plot_data, aes(x = motif, y = proportion, fill = group_label)) +
                geom_bar(stat = "identity", position = position_dodge(preserve = "single"), color = "black", linewidth = 0.2, alpha = 0.8, ...) +
                facet_wrap(~first_base, scales = "free_x", nrow = 1) +
                scale_fill_brewer(palette = "Set2") # A robust palette for groups

            plot_xlab <- "Motif"
        }
    }

    # --- 5. Common Plot Styling ---
    plot_title <- switch(representation,
        "differential" = "Differential Motif Usage",
        "split_by_base" = "Motif Proportions at Fragment Ends",
        "split_by_motif" = "Motif Proportions by Group"
    )
    plot_subtitle <- switch(representation,
        "differential" = paste("Comparison:", vals_z[1], "vs", vals_z[2]),
        "split_by_base" = "Hierarchical view by base position",
        "split_by_motif" = "Grouped by motif"
    )
    plot_ylab <- switch(representation,
        "differential" = paste("Log2 Fold Change (", vals_z[1], "/", vals_z[2], ")"),
        "Proportion" # Default for both proportion plots
    )

    final_plot <- final_plot +
        labs(title = plot_title, subtitle = plot_subtitle, y = plot_ylab, x = plot_xlab) +
        theme_bw(base_size = 11) +
        theme(
            axis.text.x = element_text(
                angle = if (representation == "split_by_motif") 90 else 0,
                vjust = 0.5,
                hjust = if (representation == "split_by_motif") 1 else 0.5
            ),
            strip.background = element_rect(fill = "grey90", color = "grey50"),
            strip.text = element_text(face = "bold"),
            panel.spacing = unit(0.2, "lines"),
            legend.position = "right",
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linetype = "dashed", color = "grey85"),
            panel.grid.minor = element_blank()
        )

    # Explicitly set the fill legend title for the split_by_motif case.
    if (representation == "split_by_motif") {
        legend_title <- if (is.null(original_col_z)) "placeholder_group" else original_col_z
        final_plot <- final_plot + labs(fill = legend_title)
    }

    return(final_plot)
}

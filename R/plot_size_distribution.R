#' Plot Fragment Size Distribution
#'
#' @description Generates a plot visualizing the distribution of fragment lengths. It allows for grouping by a
#' categorical variable and can represent the distribution as a histogram, a density plot, or an overlay of both.
#' It also displays the sample size (N) for each group in the legend.
#'
#' @param df_fragments The input dataframe containing the fragment data.
#' @param size_col A character string specifying the name of the numeric column that contains the fragment lengths.
#' @param col_z A character string specifying the name of the column to use for grouping the data. If NULL, no grouping is applied.
#' @param vals_z An optional character vector to filter and display only specific groups from 'col_z'. If NULL, all groups are used.
#' @param histo_args A named list of additional arguments to pass to 'ggplot2::geom_histogram()'.
#' @param density_args A named list of additional arguments to pass to 'ggplot2::geom_density()'.
#' @param colors_z A character vector of colors for the groups, or a single string naming an RColorBrewer palette.
#' @param show_histogram A logical value. If TRUE, a histogram layer is added.
#' @param show_density A logical value. If TRUE, a density plot layer is added.
#' @param x_limits An optional numeric vector of length 2 to set the x-axis limits (e.g., c(0, 700)).
#' @param histogram_binwidth A numeric value specifying the bin width for the histogram.
#' @param show_nuc_peaks A logical value. If TRUE, adds vertical lines for nucleosome peaks.
#'
#' @return A 'ggplot' object representing the size distribution plot.
#'
#' @importFrom dplyr %>% filter count
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density labs theme_bw theme element_text after_stat
#' @importFrom ggplot2 scale_color_manual scale_fill_manual geom_vline coord_cartesian
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
#' @examples
#' ## --- Create a dataset for demonstration ---
#' # Set a seed for reproducibility
#' set.seed(42)
#'
#' # Generate fragment sizes for two groups with different distributions
#' # "MUT" group: N=100, shorter fragments
#' mut_sizes <- rnorm(100, mean = 150, sd = 20)
#'
#' # "WT" group: N=150, centered around the mononucleosome peak
#' wt_sizes <- rnorm(150, mean = 170, sd = 25)
#'
#' # Add some larger, dinucleosomal fragments to both groups
#' di_nuc_sizes <- rnorm(30, mean = 330, sd = 30)
#'
#' # Combine into a single dataframe
#' example_df_size <- data.frame(
#'     Fragment_Size = c(mut_sizes, wt_sizes, di_nuc_sizes),
#'     Fragment_Status_Simple = c(
#'         rep("MUT", 100),
#'         rep("WT", 150),
#'         sample(c("MUT", "WT"), 30, replace = TRUE)
#'     )
#' )
#' # Ensure all fragment sizes are positive
#' example_df_size <- example_df_size[example_df_size$Fragment_Size > 0, ]
#'
#' ## --- Function Calls ---
#'
#' # 1. Default plot: A grouped density plot with nucleosome peaks shown.
#' p1 <- plot_size_distribution(example_df_size)
#' print(p1)
#'
#' # 2. Histogram plot: Show distributions as histograms instead of density curves.
#' #    We add transparency (alpha) so overlapping bars are visible.
#' p2 <- plot_size_distribution(
#'     df_fragments = example_df_size,
#'     show_histogram = TRUE,
#'     show_density = FALSE,
#'     histo_args = list(alpha = 0.6)
#' )
#' print(p2)
#'
#' # 3. Combined plot: Overlay both density curves and histograms.
#' p3 <- plot_size_distribution(
#'     df_fragments = example_df_size,
#'     show_histogram = TRUE,
#'     show_density = TRUE,
#'     histo_args = list(alpha = 0.4)
#' )
#' print(p3)
#'
#' # 4. Ungrouped and customized plot: Analyze all fragments together,
#' #    zoom in on the x-axis, and hide the nucleosome peak lines.
#' p4 <- plot_size_distribution(
#'     df_fragments = example_df_size,
#'     col_z = NULL,
#'     x_limits = c(50, 400),
#'     show_nuc_peaks = FALSE
#' )
#' print(p4)
plot_size_distribution <- function(df_fragments,
                                   size_col = "Fragment_Size",
                                   col_z = "Fragment_Status_Simple",
                                   vals_z = NULL,
                                   histo_args = list(),
                                   density_args = list(),
                                   colors_z = NULL,
                                   show_histogram = FALSE,
                                   show_density = TRUE,
                                   x_limits = c(0, 750),
                                   histogram_binwidth = 5,
                                   show_nuc_peaks = TRUE) {
    # --- 1. Input Validation and Setup ---
    if (is.null(col_z) && !is.null(vals_z)) stop("If 'col_z' is NULL, 'vals_z' must also be NULL.")
    if (!is.null(col_z) && !col_z %in% names(df_fragments)) {
        stop(sprintf("Column '%s' not found in the dataframe.", col_z))
    }
    if (!size_col %in% names(df_fragments)) {
        stop(sprintf("Size column '%s' not found in the dataframe.", size_col))
    }
    if (!is.numeric(df_fragments[[size_col]])) {
        stop(sprintf("Size column '%s' must be numeric.", size_col))
    }
    if (!show_histogram && !show_density) stop("At least one of 'show_histogram' or 'show_density' must be TRUE.")

    # --- 2. Data Preparation and Legend Labels ---
    is_grouped_analysis <- !is.null(col_z)
    if (!is_grouped_analysis) {
        df_fragments$placeholder_group <- "All Fragments"
        col_z <- "placeholder_group"
    }
    df_filtered <- df_fragments %>% filter(!is.na(.data[[size_col]]) & !is.na(.data[[col_z]]))
    if (is_grouped_analysis && !is.null(vals_z)) {
        df_filtered <- df_filtered %>% filter(.data[[col_z]] %in% vals_z)
    }
    if (nrow(df_filtered) == 0) stop("No data remains after filtering. Check 'col_z' and 'vals_z'.")
    if (is.null(vals_z)) {
        vals_z <- sort(unique(df_filtered[[col_z]]))
    }
    group_counts <- df_filtered %>% count(.data[[col_z]], name = "n")
    group_totals <- group_counts %>%
        mutate(group = factor(.data[[col_z]], levels = vals_z)) %>%
        arrange(group)
    new_group_labels <- paste0(group_totals$group, " (N=", group_totals$n, ")")
    df_filtered[[col_z]] <- factor(df_filtered[[col_z]], levels = vals_z, labels = new_group_labels)

    # --- 3. Plot Construction ---
    p <- ggplot(df_filtered, aes(x = .data[[size_col]], color = .data[[col_z]], fill = .data[[col_z]]))

    if (show_histogram) {
        # Combine default arguments with user-provided arguments
        histo_defaults <- list(mapping = aes(y = after_stat(density)), position = "identity", binwidth = histogram_binwidth)
        final_histo_args <- c(histo_defaults, histo_args)
        p <- p + do.call(geom_histogram, final_histo_args)
    }

    if (show_density) {
        # Combine default arguments with user-provided arguments
        density_defaults <- list(fill = NA, na.rm = TRUE)
        final_density_args <- c(density_defaults, density_args)
        p <- p + do.call(geom_density, final_density_args)
    }

    if (show_nuc_peaks) {
        nuc_peaks <- c(mono = 167, di = 334, tri = 501)
        p <- p + geom_vline(xintercept = nuc_peaks, linetype = "dashed", color = "grey30")
    }

    # --- 4. Color Scaling ---
    if (!is.null(colors_z)) {
        if (length(colors_z) == 1 && colors_z %in% row.names(RColorBrewer::brewer.pal.info)) {
            colors_z <- RColorBrewer::brewer.pal(n = length(vals_z), name = colors_z)
        }
        p <- p + scale_color_manual(values = colors_z) + scale_fill_manual(values = colors_z)
    }

    # --- 5. Labels and Theming ---
    y_label <- if (show_histogram || show_density) "Density" else "Count"
    p <- p + labs(title = "Fragment Size Distribution", x = "Fragment Size (bp)", y = y_label, color = "Group", fill = "Group") +
        theme_bw(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_text(face = "bold"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom"
        )
    if (!is.null(x_limits)) {
        p <- p + coord_cartesian(xlim = x_limits)
    }

    return(p)
}

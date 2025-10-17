#' Plot Fragment Size Distribution
#'
#' @description
#' Generates a plot visualizing the distribution of fragment lengths. Supports grouping by a
#' categorical variable and can represent the distribution as a histogram, a density plot,
#' or an overlay of both. The legend displays the sample size (N) per group. Groups with
#' fewer than 2 observations are automatically removed when `show_density = TRUE`.
#'
#' @param df_fragments The input dataframe containing the fragment data.
#' @param size_col A character string specifying the name of the numeric column that contains the fragment lengths.
#' @param col_z A character string specifying the name of the column to use for grouping the data. If NULL, no grouping is applied.
#' @param vals_z An optional character vector to filter and display only specific groups from `col_z`. If NULL, all groups are used.
#' @param histo_args A named list of additional arguments passed to [ggplot2::geom_histogram()].
#' @param density_args A named list of additional arguments passed to [ggplot2::geom_density()].
#' @param colors_z A character vector of colors for the groups (named or unnamed), or a single string naming an RColorBrewer palette.
#' @param show_histogram Logical. If TRUE, a histogram layer is added.
#' @param show_density Logical. If TRUE, a density plot layer is added.
#' @param x_limits Optional numeric vector of length 2 to set the x-axis limits (e.g., `c(0, 700)`).
#' @param histogram_binwidth Numeric value specifying the bin width for the histogram.
#' @param show_nuc_peaks Logical. If TRUE, adds dashed vertical lines for nucleosome peaks (mono/di/tri).
#' @param title Character or NA. Plot title; if NULL/NA/'NA'/empty, a default title is used.
#' @param output_path Character or NA. If provided and non-empty, the plot is saved to this file.
#' @param ggsave_params A named list of arguments passed to [ggplot2::ggsave()].
#'
#' @return A `ggplot` object representing the size distribution plot (invisibly `NULL` if saved).
#'
#' @importFrom dplyr %>% filter count
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density labs theme_bw theme element_text after_stat
#' @importFrom ggplot2 scale_color_manual scale_fill_manual geom_vline coord_cartesian scale_y_continuous
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom utils modifyList
#' @export
#'
#' @examples
#' ## --- Create a dataset for demonstration ---
#' set.seed(42)
#'
#' # Generate fragment sizes for two groups with different distributions
#' # 'MUT' group: N=100, shorter fragments
#' mut_sizes <- rnorm(100, mean = 150, sd = 20)
#'
#' # 'WT' group: N=150, centered around the mononucleosome peak
#' wt_sizes <- rnorm(150, mean = 170, sd = 25)
#'
#' # Add some larger, dinucleosomal fragments to both groups
#' di_nuc_sizes <- rnorm(30, mean = 330, sd = 30)
#'
#' # Combine into a single dataframe
#' example_df_size <- data.frame(
#'     Fragment_Size = c(mut_sizes, wt_sizes, di_nuc_sizes),
#'     Fragment_Status_Simple = c(
#'         rep('MUT', 100),
#'         rep('WT', 150),
#'         sample(c('MUT', 'WT'), 30, replace = TRUE)
#'     )
#' )
#' # Ensure all fragment sizes are positive
#' example_df_size <- example_df_size[example_df_size$Fragment_Size > 0, ]
#'
#' ## --- Plotting Examples ---
#'
#' # 1) Default plot: grouped density with nucleosome peaks.
#' p1 <- plot_size_distribution(example_df_size)
#' print(p1)
#'
#' # 2) Histogram only: add transparency so overlapping bars are visible.
#' p2 <- plot_size_distribution(
#'     df_fragments   = example_df_size,
#'     show_histogram = TRUE,
#'     show_density   = FALSE,
#'     histo_args     = list(alpha = 0.6)
#' )
#' print(p2)
#'
#' # 3) Combined: overlay density curves and histograms.
#' p3 <- plot_size_distribution(
#'     df_fragments   = example_df_size,
#'     show_histogram = TRUE,
#'     show_density   = TRUE,
#'     histo_args     = list(alpha = 0.4)
#' )
#' print(p3)
#'
#' # 4) Ungrouped + zoomed x-axis + no nucleosome peaks.
#' p4 <- plot_size_distribution(
#'     df_fragments   = example_df_size,
#'     col_z          = NULL,
#'     x_limits       = c(50, 400),
#'     show_nuc_peaks = FALSE,
#'     title          = 'All fragments (zoomed)'
#' )
#' print(p4)
#'
#' # 5) Custom colors using an RColorBrewer palette.
#' p5 <- plot_size_distribution(
#'     df_fragments = example_df_size,
#'     colors_z     = 'Set2',
#'     title        = 'Fragment size (Set2 palette)'
#' )
#' print(p5)
#'
#' # 6) Save to file (commented for CRAN):
#' # out_png <- file.path(tempdir(), 'size_distribution.png')
#' # plot_size_distribution(
#' #   df_fragments  = example_df_size,
#' #   output_path   = out_png,
#' #   ggsave_params = list(width = 8, height = 6, units = 'in', dpi = 300, bg = 'white'),
#' #   title         = 'Size distribution (saved)'
#' # )
plot_size_distribution <- function(df_fragments, size_col = "Fragment_Size", col_z = "Fragment_Status_Simple",
    vals_z = NULL, histo_args = list(), density_args = list(), colors_z = NULL, show_histogram = FALSE,
    show_density = TRUE, x_limits = c(0, 750), histogram_binwidth = 5, show_nuc_peaks = TRUE,
    title = NULL, output_path = NA_character_, ggsave_params = list(width = 10, height = 7,
        units = "in", dpi = 300, bg = "white")) {
    # ---- basic checks ----
    if (is.null(col_z) && !is.null(vals_z))
        stop("If 'col_z' is NULL, 'vals_z' must also be NULL.")
    if (!is.null(col_z) && !col_z %in% names(df_fragments)) {
        stop(sprintf("Column '%s' not found in the dataframe.", col_z))
    }
    if (!size_col %in% names(df_fragments)) {
        stop(sprintf("Size column '%s' not found in the dataframe.", size_col))
    }
    if (!is.numeric(df_fragments[[size_col]])) {
        stop(sprintf("Size column '%s' must be numeric.", size_col))
    }
    if (!show_histogram && !show_density) {
        stop("At least one of 'show_histogram' or 'show_density' must be TRUE.")
    }

    # ---- helpers ----
    .clean_args <- function(defaults, user) {
        # only add defaults when user didn't provide the key
        for (nm in names(defaults)) {
            if (!nm %in% names(user))
                user[[nm]] <- defaults[[nm]]
        }
        user
    }
    .resolve_group_palette <- function(keys, colors_z) {
        # keys are legend labels (e.g., 'MUT (N=123)')
        if (is.null(colors_z)) {
            return(NULL)
        }

        # Brewer palette name
        if (is.character(colors_z) && length(colors_z) == 1 && colors_z %in% rownames(RColorBrewer::brewer.pal.info)) {
            maxc <- RColorBrewer::brewer.pal.info[colors_z, "maxcolors"]
            n <- length(keys)
            if (is.na(maxc) || maxc < n)
                stop(sprintf("Palette '%s' must support at least %d colors.", colors_z, n))
            vals <- RColorBrewer::brewer.pal(n, colors_z)
            names(vals) <- keys
            return(vals)
        }

        # unnamed vector → take in order
        if (is.null(names(colors_z))) {
            if (length(colors_z) < length(keys)) {
                stop(sprintf("Provided %d colors but need %d for: %s", length(colors_z), length(keys), paste(keys, collapse = ", ")))
            }
            vals <- colors_z[seq_along(keys)]
            names(vals) <- keys
            return(vals)
        }

        # named vector → try exact match, then fallback to base group names
        # without (N=…)
        if (all(keys %in% names(colors_z))) {
            return(colors_z[keys])
        }
        base_keys <- sub(" \\(N=.*\\)$", "", keys)
        if (all(base_keys %in% names(colors_z))) {
            vals <- colors_z[base_keys]
            names(vals) <- keys
            return(vals)
        }
        missing <- setdiff(keys, names(colors_z))
        stop(sprintf("Missing color(s) for: %s", paste(missing, collapse = ", ")))
    }

    # ---- grouping ---- Save col_z for legend
    original_col_z <- col_z
    legend_name <- if (is.null(original_col_z))
        NULL else original_col_z

    is_grouped <- !is.null(col_z)
    if (!is_grouped) {
        df_fragments$..group.. <- "All Fragments"
        col_z <- "..group.."
    }
    df_filtered <- dplyr::filter(df_fragments, !is.na(.data[[size_col]]) & is.finite(.data[[size_col]]) &
        !is.na(.data[[col_z]]))

    if (is_grouped && !is.null(vals_z)) {
        df_filtered <- dplyr::filter(df_filtered, .data[[col_z]] %in% vals_z)
    }
    if (nrow(df_filtered) == 0)
        stop("No data remains after filtering. Check 'col_z' and 'vals_z'.")

    # remove groups with < 2 points for density (to avoid errors)
    if (is_grouped && show_density) {
        group_counts <- df_filtered %>%
            dplyr::count(.data[[col_z]], name = "n")
        keep <- group_counts %>%
            dplyr::filter(n >= 2) %>%
            dplyr::pull(.data[[col_z]])
        if (length(keep) < nrow(group_counts)) {
            message("Note: Groups with fewer than 2 data points were removed as they cannot be plotted (density).")
        }
        df_filtered <- dplyr::filter(df_filtered, .data[[col_z]] %in% keep)
    }
    if (nrow(df_filtered) == 0)
        stop("No data remains after density QC. Check your filters.")

    # determine final group order
    if (is.null(vals_z)) {
        vals_z <- sort(unique(as.character(df_filtered[[col_z]])))
    } else {
        vals_z <- intersect(vals_z, unique(as.character(df_filtered[[col_z]])))
    }

    # ---- legend labels ----
    if (is.null(original_col_z)) {
        # None grouped ()
        n_total <- nrow(df_filtered)
        single_lab <- paste0("All Fragments (N=", n_total, ")")
        df_filtered[[col_z]] <- factor(df_filtered[[col_z]], levels = "All Fragments",
            labels = single_lab)
    } else {
        # Grouped
        group_counts_final <- df_filtered %>%
            dplyr::count(.data[[col_z]], name = "n")
        group_totals <- group_counts_final %>%
            dplyr::mutate(group = factor(.data[[col_z]], levels = vals_z)) %>%
            dplyr::arrange(group)
        new_labels <- paste0(group_totals$group, " (N=", group_totals$n, ")")
        df_filtered[[col_z]] <- factor(df_filtered[[col_z]], levels = vals_z, labels = new_labels)
    }

    # ---- base plot ----
    final_plot <- ggplot2::ggplot(df_filtered, ggplot2::aes(x = .data[[size_col]],
        color = .data[[col_z]]))

    # histogram layer
    if (show_histogram) {
        histo_defaults <- list(mapping = ggplot2::aes(y = ggplot2::after_stat(density),
            fill = .data[[col_z]]), position = "identity", binwidth = histogram_binwidth)
        histo_args <- .clean_args(histo_defaults, histo_args)
        final_plot <- final_plot + do.call(ggplot2::geom_histogram, histo_args)
    }

    # density layer
    if (show_density) {
        dens_defaults <- list(na.rm = TRUE)
        density_args <- .clean_args(dens_defaults, density_args)
        final_plot <- final_plot + do.call(ggplot2::geom_density, density_args)
    }

    # nucleosome peaks
    if (show_nuc_peaks) {
        nuc_peaks <- c(mono = 167, di = 334, tri = 501)
        final_plot <- final_plot + ggplot2::geom_vline(xintercept = nuc_peaks, linetype = "dashed",
            color = "grey30")
    }

    # ---- colors ----
    group_keys <- levels(df_filtered[[col_z]])
    pal <- .resolve_group_palette(group_keys, colors_z)
    if (!is.null(pal)) {
        final_plot <- final_plot + ggplot2::scale_color_manual(values = pal, name = legend_name)
        if (show_histogram)
            final_plot <- final_plot + ggplot2::scale_fill_manual(values = pal, name = legend_name)
    } else {
        # force consistent legend titles
        final_plot <- final_plot + ggplot2::scale_color_discrete(name = legend_name)
        if (show_histogram)
            final_plot <- final_plot + ggplot2::scale_fill_discrete(name = legend_name)
    }

    # ---- labels, theme, axes ----
    y_lab <- if (show_histogram || show_density)
        "Density" else "Count"
    auto_title <- "Fragment Size Distribution"
    final_plot <- final_plot + ggplot2::labs(title = if (is.null(title) || is.na(title) ||
        !nzchar(title) || identical(title, "NA"))
        auto_title else title, x = "Fragment Size (bp)", y = y_lab) + ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            legend.title = ggplot2::element_text(face = "bold"), legend.position = "bottom",
            panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))  # 0% bottom gap, +10% headroom
)

    if (!is.null(x_limits)) {
        final_plot <- final_plot + ggplot2::coord_cartesian(xlim = x_limits)
    }

    # ---- save if requested ----
    if (!is.null(output_path) && is.character(output_path) && length(output_path) ==
        1L && !is.na(output_path) && nzchar(output_path)) {
        out_dir <- dirname(output_path)
        if (!dir.exists(out_dir))
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        # Merge defaults with user; user wins
        defaults <- list(width = 8, height = 6, units = "in", dpi = 300, bg = "white")
        gp <- utils::modifyList(defaults, ggsave_params)

        ggplot2::ggsave(filename = output_path, plot = final_plot, width = gp$width,
            height = gp$height, units = gp$units, dpi = gp$dpi, bg = gp$bg)
        return(invisible(NULL))
    }

    final_plot
}

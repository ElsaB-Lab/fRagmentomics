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
#'   - 'split_by_base' (default): A proportion plot with hierarchical axes, splitting motifs by each base position into facets.
#'   - 'differential': A log2 fold change plot comparing two groups.
#'   - 'split_by_motif': A proportion plot with motifs on the x-axis, with bars for different groups placed side-by-side.
#' @param ... Additional arguments passed on to 'ggplot2::geom_bar()'.
#' @param colors_z Colors for the representation:
#'   - For 'split_by_base': 4 colors for A/C/G/T, or a single RColorBrewer palette name.
#'   - For 'differential': 2 colors for 'Positive'/'Negative' (named vector or palette).
#'   - For 'split_by_motif': colors per group (palette, unnamed vector, or a vector named by group names).
#' @param title Character or NA. Plot title; if NULL/NA/'NA'/empty, a default title is used.
#' @param output_path Character or NA. If provided and non-empty, the plot is saved to this file.
#' @param ggsave_params A named list of arguments passed to 'ggplot2::ggsave()'.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr %>% filter select mutate all_of group_by summarise ungroup n count bind_rows left_join pull arrange
#' @importFrom stringr str_sub str_length str_detect
#' @importFrom tidyr drop_na pivot_wider
#' @importFrom ggh4x facet_nested
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grDevices hcl.colors
#' @importFrom grid unit
#' @importFrom stats setNames
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' ## --- Create a dataset for demonstration ---
#' set.seed(42)
#'
#' # Helper function to generate random DNA sequences with a bias
#' generate_biased_dna <- function(n_seq, len, prob) {
#'     bases <- c('A', 'C', 'G', 'T')
#'     replicate(n_seq, paste(sample(bases, len, replace = TRUE, prob = prob), collapse = ''))
#' }
#'
#' # Create 50 'MUT' fragments with a high proportion of motifs starting with 'C'
#' df_mut <- data.frame(
#'     Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'     Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'     Fragment_Status_Simple = 'MUT'
#' )
#'
#' # Create 50 'WT' fragments with a high proportion of motifs starting with 'G'
#' df_wt <- data.frame(
#'     Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
#'     Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
#'     Fragment_Status_Simple = 'WT'
#' )
#'
#' # Combine into a single dataframe
#' example_df <- rbind(df_mut, df_wt)
#'
#' ## --- Function Calls for Each Representation ---
#'
#' # 1. Hierarchical Plot (representation = 'split_by_base')
#' # This is the default. It creates nested facets for each base position.
#' p1 <- plot_motif_barplot(
#'     df_fragments   = example_df,
#'     representation = 'split_by_base'
#' )
#' print(p1)
#'
#' # You can also filter this plot to show only motifs starting with certain bases.
#' p1_filtered <- plot_motif_barplot(
#'     df_fragments   = example_df,
#'     representation = 'split_by_base',
#'     motif_start    = c('C', 'G'),
#'     title          = 'Motifs starting with C/G'
#' )
#' print(p1_filtered)
#'
#' # Optional: customize colors for the 2nd base (A/C/G/T) in split_by_base
#' p1_colors <- plot_motif_barplot(
#'     df_fragments   = example_df,
#'     representation = 'split_by_base',
#'     colors_z       = c(A = '#FD96A9', C = '#E88B00', G = '#0D539E', T = '#6CAE75')
#' )
#' print(p1_colors)
#'
#' # 2. Differential Plot (representation = 'differential')
#' # Shows log2 fold-change in motif proportions between two groups (needs exactly two groups).
#' p2 <- plot_motif_barplot(
#'     df_fragments   = example_df,
#'     representation = 'differential',
#'     vals_z         = c('MUT', 'WT'),
#'     colors_z       = c(Positive = '#66C2A5', Negative = '#E78AC3'),
#'     title          = 'MUT vs WT (log2FC)'
#' )
#' print(p2)
#'
#' # 3. Side-by-side Motif Plot (representation = 'split_by_motif')
#' # Motifs on the x-axis; bars for each group shown side-by-side.
#' p3 <- plot_motif_barplot(
#'     df_fragments   = example_df,
#'     representation = 'split_by_motif',
#'     colors_z       = 'Set2' # or a vector named by group names
#' )
#' print(p3)
#'
#' # 4. Save the default hierarchical plot (commented for CRAN)
#' # out_file1 <- file.path(tempdir(), 'motif_split_by_base.png')
#' # plot_motif_barplot(
#' #   df_fragments   = example_df,
#' #   representation = 'split_by_base',
#' #   title          = 'Motif proportions (hierarchical)',
#' #   output_path    = out_file1,
#' #   ggsave_params  = list(width = 8, height = 6, units = 'in', dpi = 300, bg = 'white')
#' # )
#'
#' # 5. Save the differential plot with custom dimensions (commented for CRAN)
#' # out_file2 <- file.path(tempdir(), 'motif_differential.png')
#' # plot_motif_barplot(
#' #   df_fragments   = example_df,
#' #   representation = 'differential',
#' #   vals_z         = c('MUT', 'WT'),
#' #   title          = 'Differential motif usage',
#' #   output_path    = out_file2,
#' #   ggsave_params  = list(width = 12, height = 8, units = 'in', dpi = 300, bg = 'white')
#' # )
plot_motif_barplot <- function(df_fragments, end_motif_5p = "Fragment_Bases_5p",
    end_motif_3p = "Fragment_Bases_3p", motif_type = "Both", motif_start = NULL,
    col_z = "Fragment_Status_Simple", vals_z = NULL, representation = "split_by_base",
    ..., colors_z = NULL, title = NULL, output_path = NA_character_, ggsave_params = list(width = 10,
        height = 6, units = "in", dpi = 300, bg = "white")) {
    motif_size <- 3  # fixed for this function

    # ---- Helpers ----
    .clean_dots <- function(dots) {
        if (!length(dots)) {
            return(dots)
        }
        keep <- (lengths(dots) > 0) & (!vapply(dots, is.null, logical(1)))
        dots <- dots[keep]
        if ("position" %in% names(dots)) {
            pos <- dots$position
            valid <- (is.character(pos) && length(pos) == 1L) || inherits(pos, "Position")
            if (!valid)
                dots$position <- NULL
        }
        dots
    }
    .resolve_palette <- function(keys, colors_z, default_vals, brewer_fallback = "Set2") {
        # default (NULL) → use provided defaults
        if (is.null(colors_z)) {
            vals <- default_vals
            if (is.null(names(vals)))
                names(vals) <- keys
            return(vals[keys])
        }
        # single Brewer palette
        if (is.character(colors_z) && length(colors_z) == 1 && colors_z %in% rownames(RColorBrewer::brewer.pal.info)) {
            maxc <- RColorBrewer::brewer.pal.info[colors_z, "maxcolors"]
            n <- length(keys)
            if (is.na(maxc) || maxc < n) {
                stop(sprintf("Palette '%s' must support at least %d colors.", colors_z,
                  n))
            }
            vals <- RColorBrewer::brewer.pal(n, colors_z)
            names(vals) <- keys
            return(vals)
        }
        # unnamed vector → in order
        if (is.null(names(colors_z))) {
            if (length(colors_z) < length(keys)) {
                stop(sprintf("Provided %d colors but need %d for: %s", length(colors_z),
                  length(keys), paste(keys, collapse = ", ")))
            }
            vals <- colors_z[seq_along(keys)]
            names(vals) <- keys
            return(vals)
        }
        # named vector → must cover keys (with a special case for group labels
        # 'X (N=...)')
        if (!all(keys %in% names(colors_z))) {
            # try matching by base names if keys are 'Group (N=..)'
            base_keys <- sub(" \\(N=.*\\)$", "", keys)
            if (all(base_keys %in% names(colors_z))) {
                vals <- colors_z[base_keys]
                names(vals) <- keys
                return(vals)
            }
            missing <- setdiff(keys, names(colors_z))
            stop(sprintf("Missing color(s) for: %s", paste(missing, collapse = ", ")))
        }
        colors_z[keys]
    }

    # ---- Validate representation ----
    valid_reps <- c("differential", "split_by_base", "split_by_motif")
    if (!representation %in% valid_reps) {
        stop(sprintf("representation' must be one of: %s", paste(valid_reps, collapse = ", ")))
    }
    original_col_z <- col_z

    # ---- Grouping logic ----
    if (is.null(col_z) && !is.null(vals_z))
        stop("If 'col_z' is NULL, 'vals_z' must also be NULL.")
    is_grouped <- !is.null(col_z)
    if (!is_grouped) {
        df_fragments$..group.. <- "All Fragments"
        col_z <- "..group.."
    } else if (!col_z %in% names(df_fragments)) {
        stop(sprintf("Column '%s' not found in the dataframe.", col_z))
    }
    if (representation == "differential" && !is_grouped) {
        stop("Differential analysis requires a grouping column 'col_z'.")
    }

    # ---- Force UPPERCASE sequences ----
    if (motif_type %in% c("Start", "Both") && end_motif_5p %in% names(df_fragments)) {
        df_fragments[[end_motif_5p]] <- toupper(df_fragments[[end_motif_5p]])
    }
    if (motif_type %in% c("End", "Both") && end_motif_3p %in% names(df_fragments)) {
        df_fragments[[end_motif_3p]] <- toupper(df_fragments[[end_motif_3p]])
    }

    # ---- Filter rows & groups ----
    df_filtered <- df_fragments %>%
        dplyr::filter(!is.na(.data[[col_z]]))
    if (is_grouped && !is.null(vals_z)) {
        df_filtered <- df_filtered %>%
            dplyr::filter(.data[[col_z]] %in% vals_z)
    }
    if (nrow(df_filtered) == 0)
        stop("No data remains after filtering. Check 'col_z' and 'vals_z'.")
    if (is.null(vals_z))
        vals_z <- sort(unique(df_filtered[[col_z]]))
    df_filtered[[col_z]] <- factor(df_filtered[[col_z]], levels = vals_z)

    # ---- Build motif table (ACGT only; length = 3) ----
    motifs_list <- list()
    if (motif_type %in% c("Start", "Both")) {
        motifs_list$start <- df_filtered %>%
            dplyr::select(group = dplyr::all_of(col_z), motif = dplyr::all_of(end_motif_5p))
    }
    if (motif_type %in% c("End", "Both")) {
        motifs_list$end <- df_filtered %>%
            dplyr::select(group = dplyr::all_of(col_z), motif = dplyr::all_of(end_motif_3p))
    }

    analysis_df <- dplyr::bind_rows(motifs_list) %>%
        tidyr::drop_na(motif) %>%
        dplyr::mutate(motif = stringr::str_sub(motif, 1, motif_size)) %>%
        dplyr::filter(stringr::str_detect(motif, "^[ACGT]{3}$"))

    if (!is.null(motif_start)) {
        analysis_df <- analysis_df %>%
            dplyr::filter(stringr::str_sub(motif, 1, 1) %in% motif_start)
    }
    if (nrow(analysis_df) == 0)
        stop("No valid 3-base motifs found in the data after filtering.")

    # ---- Proportions per group ----
    proportions_df <- analysis_df %>%
        dplyr::group_by(group) %>%
        dplyr::count(motif, name = "count") %>%
        dplyr::mutate(proportion = count/sum(count)) %>%
        dplyr::ungroup()

    # ---- Prepare plotting data ----
    base_levels <- c("A", "C", "G", "T")
    plot_data <- proportions_df %>%
        dplyr::mutate(first_base = factor(stringr::str_sub(motif, 1, 1), levels = base_levels),
            second_base = factor(stringr::str_sub(motif, 2, 2), levels = base_levels),
            third_base = factor(stringr::str_sub(motif, 3, 3), levels = base_levels),
            motif = factor(motif, levels = sort(unique(motif))))

    # ---- Build plot per representation ----
    dots <- .clean_dots(list(...))
    user_set_width <- "width" %in% names(dots)

    plot_xlab <- ""
    if (representation == "differential") {
        # Need exactly 2 groups
        if (length(vals_z) != 2)
            stop("Differential analysis requires exactly two values in 'vals_z'.")

        pseudocount <- 1e-09
        diff_data <- plot_data %>%
            dplyr::select(motif, group, proportion, first_base, second_base, third_base) %>%
            tidyr::pivot_wider(names_from = group, values_from = proportion, values_fill = 0) %>%
            dplyr::mutate(log2FC = log2((.data[[vals_z[1]]] + pseudocount)/(.data[[vals_z[2]]] +
                pseudocount)), sign = ifelse(log2FC >= 0, "Positive", "Negative"))

        # Colors: Positive/Negative
        diff_keys <- c("Positive", "Negative")
        diff_defaults <- c(Positive = "#016FB9", Negative = "#CE2D4F")
        fill_vals <- .resolve_palette(diff_keys, colors_z, diff_defaults)

        # Layer (no bottom gap, symmetric limits with +10% headroom)
        bar_args <- c(list(mapping = ggplot2::aes(x = third_base, y = log2FC, fill = sign),
            stat = "identity", color = "black", linewidth = 0.2), dots)
        if (!user_set_width)
            bar_args$width <- 1

        final_plot <- ggplot2::ggplot(diff_data) + do.call(ggplot2::geom_bar, bar_args) +
            ggplot2::geom_hline(yintercept = 0, color = "darkgrey", linewidth = 1) +
            ggh4x::facet_nested(~first_base + second_base, scales = "free_x") + ggplot2::scale_fill_manual(values = fill_vals,
            name = "Direction")

        max_abs <- max(abs(diff_data$log2FC), na.rm = TRUE)
        lim <- c(-1.1 * max_abs, 1.1 * max_abs)
        final_plot <- final_plot + ggplot2::scale_y_continuous(limits = lim, expand = c(0,
            0))
        plot_xlab <- "Third Base"
    } else {
        # common additions for proportion plots
        total_counts <- analysis_df %>%
            dplyr::count(group, name = "total_n")
        plot_data <- plot_data %>%
            dplyr::left_join(total_counts, by = "group") %>%
            dplyr::mutate(group_label = paste0(group, " (N=", total_n, ")"))

        ordered_levels <- plot_data %>%
            dplyr::mutate(group = factor(group, levels = vals_z)) %>%
            dplyr::arrange(group) %>%
            dplyr::pull(group_label) %>%
            unique()
        plot_data$group_label <- factor(plot_data$group_label, levels = ordered_levels)

        if (representation == "split_by_base") {
            # Fill by 2nd base
            base_defaults <- c(A = "#FD96A9", C = "#E88B00", G = "#0D539E", T = "#6CAE75")
            fill_vals <- .resolve_palette(base_levels, colors_z, base_defaults)

            bar_args <- c(list(mapping = ggplot2::aes(x = third_base, y = proportion,
                fill = second_base), stat = "identity", color = "black", linewidth = 0.2),
                dots)
            if (!user_set_width)
                bar_args$width <- 1

            final_plot <- ggplot2::ggplot(plot_data) + do.call(ggplot2::geom_bar,
                bar_args) + ggh4x::facet_nested(~group_label + first_base + second_base,
                scales = "free_x") + ggplot2::scale_fill_manual(values = fill_vals,
                name = "2nd Base") + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,
                0.1)))
            plot_xlab <- "Third Base"
        } else if (representation == "split_by_motif") {
            # Colors by group_label
            group_keys <- levels(plot_data$group_label)

            # default palette for groups
            if ("Set2" %in% rownames(RColorBrewer::brewer.pal.info) && length(group_keys) <=
                RColorBrewer::brewer.pal.info["Set2", "maxcolors"]) {
                grp_defaults <- RColorBrewer::brewer.pal(length(group_keys), "Set2")
                names(grp_defaults) <- group_keys
            } else {
                grp_defaults <- grDevices::hcl.colors(length(group_keys))
                names(grp_defaults) <- group_keys
            }
            fill_vals <- .resolve_palette(group_keys, colors_z, grp_defaults)

            bar_args <- c(list(mapping = ggplot2::aes(x = motif, y = proportion,
                fill = group_label), stat = "identity", color = "black", linewidth = 0.2),
                dots)
            if (!user_set_width)
                bar_args$width <- 0.9
            if (!"position" %in% names(dots))
                bar_args$position <- ggplot2::position_dodge(preserve = "single")

            final_plot <- ggplot2::ggplot(plot_data) + do.call(ggplot2::geom_bar,
                bar_args) + ggplot2::facet_wrap(~first_base, scales = "free_x", nrow = 1) +
                ggplot2::scale_fill_manual(values = fill_vals, name = if (is.null(original_col_z))
                  NULL else original_col_z) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,
                0.1)))
            plot_xlab <- "Motif"
        }
    }

    # ---- Titles & theme ----
    auto_title <- switch(representation, differential = "Differential Motif Usage",
        split_by_base = "Motif Proportions at Fragment Ends", split_by_motif = "Motif Proportions by Group")
    plot_subtitle <- switch(representation, differential = paste("Comparison:", vals_z[1],
        "vs", vals_z[2]), split_by_base = "Hierarchical view by base position", split_by_motif = "Grouped by motif")
    plot_ylab <- if (representation == "differential") {
        paste0("Log2 Fold Change (", vals_z[1], "/", vals_z[2], ")")
    } else {
        "Proportion"
    }

    final_plot <- final_plot + ggplot2::labs(title = if (is.null(title) || is.na(title) ||
        !nzchar(title) || identical(title, "NA"))
        auto_title else title, subtitle = plot_subtitle, y = plot_ylab, x = plot_xlab) + ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = if (representation ==
            "split_by_motif")
            90 else 0, vjust = 0.5, hjust = if (representation == "split_by_motif")
            1 else 0.5), strip.background = ggplot2::element_rect(fill = "white", color = "black",
            linewidth = 0.5), strip.text = ggplot2::element_text(face = "bold"),
            panel.spacing = grid::unit(0.2, "lines"), legend.position = "right",
            panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(linetype = "dashed",
                color = "grey85"), panel.grid.minor = ggplot2::element_blank())

    # ---- Save if requested ----
    if (!is.null(output_path) && is.character(output_path) && length(output_path) ==
        1L && !is.na(output_path) && nzchar(output_path)) {
        out_dir <- dirname(output_path)
        if (!dir.exists(out_dir))
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave(filename = output_path, plot = final_plot, width = ggsave_params$width,
            height = ggsave_params$height, units = ggsave_params$units, dpi = ggsave_params$dpi,
            bg = ggsave_params$bg)
        return(invisible(NULL))
    }

    final_plot
}

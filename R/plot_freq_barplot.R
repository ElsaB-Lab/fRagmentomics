#' Plot Overall Nucleotide Frequency
#'
#' @description
#' Creates a bar plot to compare the overall proportion of each nucleotide
#' (A/C/G/T; optional 'Other') in the end motifs. Error bars show 95% confidence intervals.
#'
#' @param df_fragments The input data frame containing fragment sequence data.
#' @param end_motif_5p Character string. Column name for 5' end sequences.
#' @param end_motif_3p Character string. Column name for 3' end sequences.
#' @param motif_type Character string. Which ends to analyze: 'Start', 'End', or 'Both'.
#' @param motif_size A single integer (>= 1) specifying the k-mer length to analyze.
#' @param col_z Character string or 'NULL'. Column name for grouping. If 'NULL', all fragments are pooled.
#' @param vals_z A character vector of group names from 'col_z' to include.
#'   If 'NULL', all unique groups in 'col_z' are used.
#' @param ... Additional aesthetics/arguments passed to [ggplot2::geom_col()] and [ggplot2::geom_errorbar()].
#'   (e.g., 'alpha', 'position', or 'width').
#' @param colors_z A character vector of colors for the groups, or a single
#'   RColorBrewer palette name (e.g., 'Set2'). Named vectors are aligned to 'vals_z'.
#' @param title Character or 'NA'. Plot title. If 'NULL', 'NA', or empty, a default title is used.
#' @param output_path Character or 'NA'. If provided and non-empty, the plot is saved to this path.
#' @param ggsave_params A named list of arguments passed to [ggplot2::ggsave()].
#' @param show_pvalue Logical. If 'TRUE' and there are at least two groups, append a global
#'   Chi-squared p-value to the caption.
#' @param drop_non_acgt Logical. If 'FALSE', characters other than A/C/G/T are tallied into an
#'   'Other' category.
#'
#' @return A 'ggplot' object. If 'output_path' is provided and non-empty, the plot is saved
#'   to file and the function returns 'invisible(NULL)'.
#'
#' @importFrom dplyr filter select mutate all_of group_by summarise ungroup bind_rows distinct arrange
#' @importFrom stringr str_sub str_count
#' @importFrom tidyr pivot_longer
#' @importFrom stats prop.test chisq.test xtabs na.omit
#' @importFrom purrr map2 map_dbl
#' @importFrom scales percent_format pvalue
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' ## --- Create a dataset for demonstration ---
#' set.seed(42)
#'
#' # Helper to generate random DNA sequences with base bias
#' generate_biased_dna <- function(n_seq, len, prob) {
#'     bases <- c('A', 'C', 'G', 'T')
#'     replicate(n_seq, paste(sample(bases, len, replace = TRUE, prob = prob), collapse = ''))
#' }
#'
#' # 50 'MUT' fragments biased toward 'C'
#' df_mut <- data.frame(
#'     Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'     Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'     Fragment_Status_Simple = 'MUT'
#' )
#'
#' # 50 'WT' fragments biased toward 'G'
#' df_wt <- data.frame(
#'     Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
#'     Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
#'     Fragment_Status_Simple = 'WT'
#' )
#'
#' # Combine into a single data frame
#' example_df <- rbind(df_mut, df_wt)
#'
#' ## --- Function calls ---
#'
#' # 1) Default plot: compare MUT vs WT for 3-mers from both ends
#' p1 <- plot_freq_barplot(example_df)
#' print(p1)
#'
#' # 2) First-nucleotide only (k = 1) on 5' end, with custom colors
#' p2 <- plot_freq_barplot(
#'     df_fragments = example_df,
#'     motif_type   = 'Start',
#'     motif_size   = 1,
#'     colors_z     = c('MUT' = '#d95f02', 'WT' = '#1b9e77'),
#'     title        = "5' First Base Composition"
#' )
#' print(p2)
#'
#' # 3) Ungrouped: overall nucleotide frequencies across all fragments
#' p3 <- plot_freq_barplot(example_df, col_z = NULL, title = 'Overall Composition')
#' print(p3)
#'
#' # 4) Subset of groups (if you had >2 groups, e.g., 'MUT', 'WT', 'AMB')
#' p4 <- plot_freq_barplot(
#'     df_fragments = example_df,
#'     vals_z       = c('MUT', 'WT')
#' )
#' print(p4)
#'
#' # 5) Include non-ACGT characters tallied as 'Other'
#' example_df$Fragment_Bases_5p[1:3] <- c('NNNNNNNNNN', 'ACGTNACGTN', 'TTTNAAAAAA')
#' p5 <- plot_freq_barplot(example_df,
#'     motif_size = 2, drop_non_acgt = FALSE,
#'     title = "Including 'Other' (non-ACGT)"
#' )
#' print(p5)
#'
#' # 6) Save to file with specific dimensions
#' # plot_freq_barplot(
#' #   df_fragments = example_df,
#' #   output_path  = file.path(tempdir(), 'nucleotide_frequency.png'),
#' #   ggsave_params = list(width = 7, height = 5, units = 'in', dpi = 300, bg = 'white')
#' # )
plot_freq_barplot <- function(df_fragments, end_motif_5p = "Fragment_Bases_5p", end_motif_3p = "Fragment_Bases_3p",
    motif_type = "Both", motif_size = 3, col_z = "Fragment_Status_Simple", vals_z = NULL,
    ..., colors_z = "Dark2", title = NULL, output_path = NA_character_, ggsave_params = list(width = 14,
        height = 5, units = "in", dpi = 300, bg = "white"), show_pvalue = FALSE,
    drop_non_acgt = TRUE) {
    # ---- Capture & sanitize dots to avoid duplicate/invalid args ----
    dots <- list(...)
    if (length(dots)) {
        keep <- (lengths(dots) > 0) & (!vapply(dots, is.null, logical(1)))
        dots <- dots[keep]
    }
    if ("position" %in% names(dots)) {
        pos <- dots$position
        valid_pos <- (is.character(pos) && length(pos) == 1L) || inherits(pos, "Position")
        if (!valid_pos)
            dots$position <- NULL
    }
    user_set_width <- "width" %in% names(dots)

    # ---- Basic checks ----
    stopifnot(is.data.frame(df_fragments))
    if (!motif_type %in% c("Start", "End", "Both")) {
        stop("otif_type' must be one of 'Start', 'End', or 'Both'.")
    }
    if (!is.numeric(motif_size) || length(motif_size) != 1L || motif_size < 1L) {
        stop("otif_size' must be a single integer >= 1.")
    }

    # Sequence column presence
    if (motif_type %in% c("Start", "Both") && !end_motif_5p %in% names(df_fragments)) {
        stop(sprintf("Column '%s' not found in 'df_fragments'.", end_motif_5p))
    }
    if (motif_type %in% c("End", "Both") && !end_motif_3p %in% names(df_fragments)) {
        stop(sprintf("Column '%s' not found in 'df_fragments'.", end_motif_3p))
    }

    # ---- Grouping ----
    if (is.null(col_z) && !is.null(vals_z)) {
        stop("If 'col_z' is NULL, 'vals_z' must also be NULL.")
    }
    is_grouped <- !is.null(col_z)
    if (!is_grouped) {
        df_fragments$..group.. <- "All Fragments"
        col_z <- "..group.."
    } else if (!col_z %in% names(df_fragments)) {
        stop(sprintf("Column '%s' not found in the dataframe.", col_z))
    }

    # ---- Uppercase sequences (forced) ----
    if (motif_type %in% c("Start", "Both")) {
        df_fragments[[end_motif_5p]] <- toupper(df_fragments[[end_motif_5p]])
    }
    if (motif_type %in% c("End", "Both")) {
        df_fragments[[end_motif_3p]] <- toupper(df_fragments[[end_motif_3p]])
    }

    # Filter groups/levels
    df_filtered <- dplyr::filter(df_fragments, !is.na(.data[[col_z]]))
    if (is_grouped && !is.null(vals_z)) {
        df_filtered <- dplyr::filter(df_filtered, .data[[col_z]] %in% vals_z)
    }
    if (nrow(df_filtered) == 0)
        stop("No data remains after filtering. Check 'col_z' and 'vals_z'.")
    if (is.null(vals_z))
        vals_z <- sort(unique(df_filtered[[col_z]]))
    df_filtered[[col_z]] <- factor(df_filtered[[col_z]], levels = vals_z)

    # ---- Motif sizing ----
    seq_lens <- integer(0)
    if (motif_type %in% c("Start", "Both")) {
        seq_lens <- c(seq_lens, nchar(stats::na.omit(df_filtered[[end_motif_5p]])))
    }
    if (motif_type %in% c("End", "Both")) {
        seq_lens <- c(seq_lens, nchar(stats::na.omit(df_filtered[[end_motif_3p]])))
    }
    if (length(seq_lens) == 0)
        stop("All selected sequence columns are empty or NA.")
    max_len <- min(seq_lens)
    if (motif_size > max_len) {
        warning(sprintf("Requested 'motif_size' (%d) exceeds shortest sequence (%d). Using %d.",
            motif_size, max_len, max_len))
        motif_size <- max_len
    }

    # ---- Build motif table ----
    motifs_list <- list()
    if (motif_type %in% c("Start", "Both")) {
        motifs_list$start <- df_filtered |>
            dplyr::select(group = dplyr::all_of(col_z), motif = dplyr::all_of(end_motif_5p)) |>
            dplyr::mutate(motif = stringr::str_sub(motif, 1L, motif_size))
    }
    if (motif_type %in% c("End", "Both")) {
        motifs_list$end <- df_filtered |>
            dplyr::select(group = dplyr::all_of(col_z), motif = dplyr::all_of(end_motif_3p)) |>
            dplyr::mutate(motif = stringr::str_sub(motif, -motif_size, -1L))
    }
    dat <- dplyr::bind_rows(motifs_list) |>
        dplyr::filter(!is.na(motif))

    # ---- Count bases ----
    count_df <- dat |>
        dplyr::mutate(a_count = stringr::str_count(motif, "A"), c_count = stringr::str_count(motif,
            "C"), g_count = stringr::str_count(motif, "G"), t_count = stringr::str_count(motif,
            "T"), Other = if (!drop_non_acgt)
            nchar(motif) - (a_count + c_count + g_count + t_count) else 0L) |>
        dplyr::select(group, a_count, c_count, g_count, t_count, Other) |>
        tidyr::pivot_longer(cols = c(a_count, c_count, g_count, t_count, Other),
            names_to = "nuc_key", values_to = "count") |>
        dplyr::mutate(nucleotide = dplyr::recode(nuc_key, a_count = "A", c_count = "C",
            g_count = "G", t_count = "T", Other = "Other")) |>
        dplyr::select(group, nucleotide, count)

    if (drop_non_acgt) {
        count_df <- dplyr::filter(count_df, nucleotide %in% c("A", "C", "G", "T"))
    }

    count_df <- count_df |>
        dplyr::group_by(group, nucleotide) |>
        dplyr::summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop")

    # --- Drop 'Other' entirely if it has zero counts (prevents NA facet) ---
    if (!drop_non_acgt) {
        has_other <- any(count_df$nucleotide == "Other" & count_df$total_count >
            0)
        if (!has_other) {
            count_df <- dplyr::filter(count_df, nucleotide != "Other")
        }
    }

    # ---- Frequencies & CIs ----
    plot_data <- count_df |>
        dplyr::group_by(group) |>
        dplyr::mutate(total_bases_in_group = sum(total_count)) |>
        dplyr::ungroup() |>
        dplyr::mutate(frequency = total_count/total_bases_in_group, ci = purrr::map2(total_count,
            total_bases_in_group, ~stats::prop.test(.x, .y)$conf.int), ci_low = purrr::map_dbl(ci,
            ~.x[1]), ci_high = purrr::map_dbl(ci, ~.x[2]))

    # Factor levels (include 'Other' only if present)
    nuc_levels <- c("A", "C", "G", "T")
    if (!drop_non_acgt && any(plot_data$nucleotide == "Other")) {
        nuc_levels <- c(nuc_levels, "Other")
    }
    plot_data <- dplyr::filter(plot_data, nucleotide %in% nuc_levels)
    plot_data$nucleotide <- factor(plot_data$nucleotide, levels = nuc_levels)
    plot_data$group <- factor(plot_data$group, levels = vals_z)
    y_max <- min(1, max(plot_data$ci_high, na.rm = TRUE) * 1.1)

    # ---- Labels & global χ² ----
    group_totals <- dplyr::distinct(plot_data, group, total_bases_in_group) |>
        dplyr::arrange(factor(group, levels = vals_z))
    x_labels <- paste0(group_totals$group, " (N=", group_totals$total_bases_in_group,
        ")")
    names(x_labels) <- as.character(group_totals$group)

    plot_caption <- "Bars show overall proportion with 95% confidence intervals."
    if (is_grouped && length(vals_z) >= 2 && show_pvalue) {
        # drop rows where a nucleotide never appears across groups
        tmp <- dplyr::filter(plot_data, total_count > 0)
        if (dplyr::n_distinct(tmp$nucleotide) >= 2) {
            contingency <- stats::xtabs(total_count ~ nucleotide + group, data = tmp)
            chi <- stats::chisq.test(contingency, simulate.p.value = TRUE, B = 5000L)
            plot_caption <- paste0(plot_caption, "\nGlobal comparison (Chi-squared), final_plot = ",
                scales::pvalue(chi$p.value))
        }
    }

    # ---- Colors (by group) ----
    if (length(colors_z) == 1L && colors_z %in% rownames(RColorBrewer::brewer.pal.info)) {
        max_n <- RColorBrewer::brewer.pal.info[colors_z, "maxcolors"]
        if (length(vals_z) > max_n) {
            warning(sprintf("Palette '%s' supports up to %d colors; you have %d groups. Reusing colors.",
                colors_z, max_n, length(vals_z)))
            pal <- RColorBrewer::brewer.pal(max_n, colors_z)
            colors_vec <- rep(pal, length.out = length(vals_z))
        } else {
            colors_vec <- RColorBrewer::brewer.pal(length(vals_z), colors_z)
        }
        names(colors_vec) <- vals_z
    } else {
        colors_vec <- colors_z
        if (is.null(names(colors_vec))) {
            if (length(colors_vec) != length(vals_z)) {
                stop("If 'colors_z' is an unnamed vector, its length must match the number of groups.")
            }
            names(colors_vec) <- vals_z
        } else {
            colors_vec <- colors_vec[vals_z]
            if (any(is.na(colors_vec)))
                stop("Named 'colors_z' must provide a color for each group level.")
        }
    }

    # ---- Title ----
    auto_title <- paste0("Overall Nucleotide Frequency (", motif_size, "-mers)")
    if (is.null(title) || is.na(title) || !nzchar(title) || identical(title, "NA")) {
        title_txt <- auto_title
    } else {
        title_txt <- title
    }

    # ---- Build layers with do.call (clean 'dots') ----
    col_args <- c(list(mapping = NULL, data = NULL), dots)
    if (!user_set_width)
        col_args$width <- 0.7
    geom_col_layer <- do.call(ggplot2::geom_col, col_args)

    geom_errorbar_layer <- ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high),
        width = 0.25)

    # ---- Plot ----
    final_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = frequency,
        fill = group)) + geom_col_layer + geom_errorbar_layer + ggplot2::facet_wrap(~nucleotide,
        scales = "free_x", nrow = 1) + ggplot2::scale_x_discrete(labels = x_labels) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,
            y_max), expand = c(0, 0)) + ggplot2::scale_fill_manual(values = colors_vec,
        name = if (!is_grouped)
            NULL else col_z, labels = x_labels) + ggplot2::theme_bw(base_size = 14) + ggplot2::labs(title = title_txt,
        x = NULL, y = "Overall Frequency", caption = plot_caption) + ggplot2::theme(legend.position = "right",
        axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "italic"),
        strip.text = ggplot2::element_text(face = "bold", size = 12), strip.background = ggplot2::element_rect(fill = "white",
            color = NA), panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(size = 11), legend.text = ggplot2::element_text(size = 10))

    # ---- Save if requested ----
    if (is.null(output_path) || !is.character(output_path) || length(output_path) !=
        1L || is.na(output_path) || !nzchar(output_path)) {
        return(final_plot)
    }

    ggplot2::ggsave(filename = output_path, plot = final_plot, width = ggsave_params$width,
        height = ggsave_params$height, units = ggsave_params$units, dpi = ggsave_params$dpi,
        bg = ggsave_params$bg)
    invisible(NULL)
}

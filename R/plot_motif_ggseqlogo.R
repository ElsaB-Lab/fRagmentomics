#' Plot sequence motif composition
#'
#' @description
#' Creates a sequence logo plot showing the proportion of each nucleotide at each
#' position, with flexible grouping/faceting.
#'
#' @param df_fragments Data frame containing fragment sequence data.
#' @param end_motif_5p Character. Column name for 5' end sequences.
#' @param end_motif_3p Character. Column name for 3' end sequences.
#' @param motif_type Character. One of 'Start', 'End', or 'Both'.
#' @param motif_size Integer (>=1). Length of the motif to analyze.
#' @param col_z Character or NULL. Grouping/faceting column. If NULL, all fragments are pooled.
#' @param vals_z Character vector or NULL. Subset of groups from 'col_z' to include.
#'   If NULL, all unique groups are used.
#' @param colors_z NULL (use ggseqlogo defaults), a single RColorBrewer palette name
#'   (e.g., 'Dark2'), or a named vector for 'A/C/G/T', e.g.
#'   'c(A='#1B9E77', C='#D95F02', G='#7570B3', T='#E7298A')'.
#' @param title Character or NA. Plot title; if NULL/NA/'NA'/empty, a default title is used.
#' @param output_path Character or NA. If provided and non-empty, the plot is saved to this file.
#' @param ggsave_params Named list passed to [ggplot2::ggsave()].
#' @param ... Extra arguments forwarded to [ggseqlogo::ggseqlogo()] (e.g., 'stack_width', 'font', or 'col_scheme').
#'
#' @return A 'ggplot' object (invisibly NULL if saved).
#'
#' @importFrom dplyr filter group_split
#' @importFrom magrittr %>%
#' @importFrom purrr map map_chr map_int set_names
#' @importFrom stringr str_detect str_sub str_split
#' @importFrom tidyr crossing
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot theme_void labs aes geom_text theme element_text element_rect element_blank margin geom_vline
#' @importFrom ggseqlogo ggseqlogo make_col_scheme
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
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
#' # 50 'MUT' fragments biased toward 'C' at the ends
#' df_mut <- data.frame(
#'     Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'     Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'     Fragment_Status_Simple = 'MUT'
#' )
#'
#' # 50 'WT' fragments biased toward 'G' at the ends
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
#' # 1) Default plot: 3-mer from both 5' and 3' ends, separated by a dash,
#' #    faceted by group ('MUT' and 'WT').
#' p1 <- plot_qqseqlogo_meme(example_df)
#' print(p1)
#'
#' # 2) Single-end motif: 5-mer from the 5' end only.
#' p2 <- plot_qqseqlogo_meme(
#'     df_fragments = example_df,
#'     motif_type   = 'Start',
#'     motif_size   = 5,
#'     title        = '5' motif (k=5)'
#' )
#' print(p2)
#'
#' # 3) Custom colors using an RColorBrewer palette (first 4 colors mapped to A/C/G/T).
#' #    Note: the '-' separator in 'Both' is not colored.
#' p3 <- plot_qqseqlogo_meme(
#'     df_fragments = example_df,
#'     motif_type   = 'Both',
#'     motif_size   = 3,
#'     colors_z     = 'Dark2',
#'     title        = 'Both ends (palette = Dark2)'
#' )
#' print(p3)
#'
#' # 4) Fully custom nucleotide colors (named vector).
#' custom_cols <- c(A = '#1B9E77', C = '#D95F02', G = '#7570B3', T = '#E7298A')
#' p4 <- plot_qqseqlogo_meme(
#'     df_fragments = example_df,
#'     motif_type   = 'Start',
#'     motif_size   = 3,
#'     colors_z     = custom_cols,
#'     title        = 'Custom nucleotide colors'
#' )
#' print(p4)
#'
#' # 5) Ungrouped: analyze all fragments together (single facet).
#' p5 <- plot_qqseqlogo_meme(example_df, col_z = NULL, title = 'All fragments pooled')
#' print(p5)
#'
#' # 6) Passing extra ggseqlogo options via '...' (e.g., stack width and font)
#' p6 <- plot_qqseqlogo_meme(
#'     df_fragments = example_df,
#'     motif_type   = 'End',
#'     motif_size   = 4,
#'     stack_width  = 0.9,
#'     font         = 'helvetica_regular',
#'     title        = '3' motif (k=4, custom stack width)'
#' )
#' print(p6)
#'
#' # 7) Save to file (commented out for CRAN)
#' # out_file <- file.path(tempdir(), 'motif_logo.png')
#' # plot_qqseqlogo_meme(
#' #   df_fragments  = example_df,
#' #   motif_type    = 'Both',
#' #   motif_size    = 3,
#' #   title         = 'Saved motif logo',
#' #   output_path   = out_file,
#' #   ggsave_params = list(width = 7, height = 5, units = 'in', dpi = 300, bg = 'white')
#' # )
#'
plot_qqseqlogo_meme <- function(df_fragments, end_motif_5p = "Fragment_Bases_5p",
    end_motif_3p = "Fragment_Bases_3p", motif_type = "Both", motif_size = 3, col_z = "Fragment_Status_Simple",
    vals_z = NULL, colors_z = NULL, title = NULL, output_path = NA_character_, ggsave_params = list(width = 12,
        height = 6, units = "in", dpi = 300, bg = "white"), ...) {
    stopifnot(is.data.frame(df_fragments))
    if (!motif_type %in% c("Start", "End", "Both"))
        stop("motif_type' must be 'Start', 'End', or 'Both'.")
    if (!is.numeric(motif_size) || length(motif_size) != 1L || motif_size < 1L)
        stop("motif_size' must be a single integer >= 1.")
    if (is.null(col_z) && !is.null(vals_z))
        stop("If 'col_z' is NULL, 'vals_z' must also be NULL.")
    is_grouped <- !is.null(col_z)
    if (!is_grouped) {
        df_fragments$..group.. <- "All Fragments"
        col_z <- "..group.."
    } else if (!col_z %in% names(df_fragments)) {
        stop(sprintf("Column '%s' not found in the dataframe.", col_z))
    }

    # Uppercase sequences (forced)
    if (motif_type %in% c("Start", "Both") && end_motif_5p %in% names(df_fragments)) {
        df_fragments[[end_motif_5p]] <- toupper(df_fragments[[end_motif_5p]])
    }
    if (motif_type %in% c("End", "Both") && end_motif_3p %in% names(df_fragments)) {
        df_fragments[[end_motif_3p]] <- toupper(df_fragments[[end_motif_3p]])
    }

    # Filter rows & levels
    df_filtered <- dplyr::filter(df_fragments, !is.na(.data[[col_z]]))
    if (is_grouped && !is.null(vals_z))
        df_filtered <- dplyr::filter(df_filtered, .data[[col_z]] %in% vals_z)
    if (nrow(df_filtered) == 0)
        stop("No data remains after filtering. Check 'col_z' and 'vals_z'.")
    if (is.null(vals_z))
        vals_z <- unique(df_filtered[[col_z]])
    df_filtered[[col_z]] <- factor(df_filtered[[col_z]], levels = vals_z)

    # Cap motif_size
    min_len_available <- Inf
    if (motif_type %in% c("Start", "Both")) {
        s5 <- stats::na.omit(df_filtered[[end_motif_5p]])
        if (length(s5) > 0)
            min_len_available <- min(min_len_available, min(nchar(s5)))
    }
    if (motif_type %in% c("End", "Both")) {
        s3 <- stats::na.omit(df_filtered[[end_motif_3p]])
        if (length(s3) > 0)
            min_len_available <- min(min_len_available, min(nchar(s3)))
    }
    if (!is.finite(min_len_available))
        stop("All selected sequence columns are empty or NA.")
    if (motif_size > min_len_available) {
        warning(sprintf("Requested 'motif_size' (%d) exceeds shortest sequence (%d). Using %d.",
            motif_size, min_len_available, min_len_available), call. = FALSE)
        motif_size <- min_len_available
    }

    # Extract motifs per group (A/C/G/T only; '-' allowed as separator for
    # 'Both')
    process_group <- function(group_df) {
        group_name <- as.character(unique(group_df[[col_z]]))
        start_motifs <- character(0)
        end_motifs <- character(0)

        if (motif_type %in% c("Start", "Both")) {
            s <- stats::na.omit(group_df[[end_motif_5p]])
            if (length(s)) {
                n_short <- sum(nchar(s) < motif_size)
                if (n_short > 0) {
                  warning(sprintf("For group '%s', removed %d 5' sequences shorter than motif_size (%d).",
                    group_name, n_short, motif_size), call. = FALSE)
                }
                s <- s[nchar(s) >= motif_size]
                if (length(s))
                  start_motifs <- substr(s, 1, motif_size)
                # filtration non-ACGT
                start_motifs <- start_motifs[!stringr::str_detect(start_motifs, "[^ACGT]")]
            }
        }

        if (motif_type %in% c("End", "Both")) {
            s <- stats::na.omit(group_df[[end_motif_3p]])
            if (length(s)) {
                n_short <- sum(nchar(s) < motif_size)
                if (n_short > 0) {
                  warning(sprintf("For group '%s', removed %d 3' sequences shorter than motif_size (%d).",
                    group_name, n_short, motif_size), call. = FALSE)
                }
                s <- s[nchar(s) >= motif_size]
                if (length(s))
                  end_motifs <- stringr::str_sub(s, -motif_size, -1L)
                # filtration non-ACGT
                end_motifs <- end_motifs[!stringr::str_detect(end_motifs, "[^ACGT]")]
            }
        }

        if (motif_type == "Both") {
            n <- min(length(start_motifs), length(end_motifs))
            if (n == 0) {
                return(character(0))
            }
            paste0(start_motifs[seq_len(n)], "-", end_motifs[seq_len(n)])
        } else if (motif_type == "Start")
            start_motifs else end_motifs
    }

    groups <- df_filtered %>%
        dplyr::group_split(.data[[col_z]])
    motifs <- purrr::map(groups, process_group)
    names(motifs) <- purrr::map_chr(groups, ~as.character(unique(.x[[col_z]])))
    motifs <- motifs[purrr::map_int(motifs, length) > 0]
    if (!length(motifs)) {
        return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No data to display"))
    }
    names(motifs) <- paste0(names(motifs), " (N=", purrr::map_int(motifs, length),
        ")")

    # --- Color scheme (default ggseqlogo if NULL) ---
    dots <- list(...)
    color_scheme <- NULL

    if (!is.null(colors_z)) {
        # Single RColorBrewer palette name
        if (is.character(colors_z) && length(colors_z) == 1 && colors_z %in% rownames(RColorBrewer::brewer.pal.info)) {
            maxc <- RColorBrewer::brewer.pal.info[colors_z, "maxcolors"]
            if (is.na(maxc) || maxc < 4) {
                stop(sprintf("Palette '%s' must provide at least 4 colors.", colors_z))
            }
            pal <- RColorBrewer::brewer.pal(4, colors_z)
            names(pal) <- c("A", "C", "G", "T")
            color_scheme <- ggseqlogo::make_col_scheme(chars = names(pal), cols = unname(pal))

            # Named vector (expecting A/C/G/T; names are normalized to
            # uppercase)
        } else if (is.character(colors_z) && !is.null(names(colors_z))) {
            req <- c("A", "C", "G", "T")
            nm <- toupper(names(colors_z))
            names(colors_z) <- nm
            if (!all(req %in% nm)) {
                stop(sprintf("Custom colors are incomplete: missing %s.", paste(setdiff(req,
                  nm), collapse = ", ")))
            }
            color_scheme <- ggseqlogo::make_col_scheme(chars = req, cols = unname(colors_z[req]))

            # Unnamed vector of 4 colors → mapped in order A, C, G, T
        } else if (is.character(colors_z) && is.null(names(colors_z))) {
            if (length(colors_z) != 4L) {
                stop(sprintf("An unnamed vector must contain exactly 4 colors (received: %d).",
                  length(colors_z)))
            }
            pal <- stats::setNames(colors_z, c("A", "C", "G", "T"))
            color_scheme <- ggseqlogo::make_col_scheme(chars = names(pal), cols = unname(pal))
        } else {
            stop("'colors_z' must be NULL, a valid RColorBrewer palette name, ",
                "a named vector for A/C/G/T, or an unnamed vector of 4 colors.")
        }
    }

    # Let a user-supplied col_scheme in ... take precedence
    if (!is.null(dots$col_scheme))
        color_scheme <- NULL

    # ---- CALL ggseqlogo ----
    args <- c(list(motifs), list(method = "prob"), if (is.null(dots$stack_width)) list(stack_width = 0.95) else list(),
        if (is.null(dots$font)) list(font = "helvetica_regular") else list(), if (!is.null(color_scheme)) list(col_scheme = color_scheme) else list(),
        dots)
    base_logo <- withCallingHandlers(do.call(ggseqlogo::ggseqlogo, args), warning = function(w) {
        msg <- conditionMessage(w)  # Remove this exact warning
        if (grepl("The `<scale>` argument of `guides\\(\\)` cannot be `FALSE`", msg)) {
            invokeRestart("muffleWarning")
        }
    })

    final_plot <- base_logo + ggplot2::labs(title = if (is.null(title) || is.na(title) ||
        !nzchar(title) || identical(title, "NA")) {
        "Sequence Motif Composition"
    } else {
        title
    }, y = "Frequency", x = "Position")

    # Custom x labels & theme (white strips), as before
    motif_len <- nchar(motifs[[1]][1])
    custom_labels <- NULL
    if (motif_type == "Start") {
        custom_labels <- as.character(seq_len(motif_len))
    } else if (motif_type == "End") {
        custom_labels <- as.character((-motif_len):-1)
    } else if (motif_type == "Both")
        custom_labels <- as.character(c(seq_len(motif_size), "", (-motif_size):-1))

    facet_names <- names(motifs)
    base_annotation_df <- tibble::tibble(x_pos = seq_len(motif_len), label = custom_labels)
    annotation_df <- tidyr::crossing(group = factor(facet_names, levels = facet_names),
        base_annotation_df)

    final_plot$coordinates$clip <- "off"
    final_plot <- final_plot + ggplot2::geom_text(data = annotation_df, ggplot2::aes(x = x_pos,
        y = -0.1, label = label), size = 5, inherit.aes = FALSE) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
        face = "bold", size = 16), strip.text = ggplot2::element_text(face = "bold",
        size = 14), axis.title.x = ggplot2::element_text(size = 14, margin = ggplot2::margin(t = 15)),
        axis.title.y = ggplot2::element_text(size = 14), strip.background = ggplot2::element_rect(fill = "white",
            color = NA), panel.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(), plot.margin = ggplot2::margin(t = 5.5,
            r = 5.5, b = 20, l = 5.5, "pt"))

    if (motif_type == "Both") {
        final_plot <- final_plot + ggplot2::geom_vline(xintercept = motif_size +
            1, linetype = "dashed", color = "grey40", linewidth = 0.8)
    }

    # Save if requested
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

#' Plot sequence motif composition
#'
#' @description Creates a sequence logo plot showing the proportion of each nucleotide at each specified position, with flexible
#' grouping and faceting.
#'
#' @param df_fragments The input dataframe containing fragment sequence data.
#' @param end_motif_5p Character string. The column name for 5' end sequences.
#' @param end_motif_3p Character string. The column name for 3' end sequences.
#' @param motif_type Character string. Which ends to analyze: 'Start', 'End', or 'Both'.
#' @param motif_size A single integer specifying the length of the motif.
#' @param col_z Character string. The column name for grouping/faceting. If NULL, no grouping is applied.
#' @param vals_z A character vector of group names to include. If NULL, all groups in 'col_z' are used.
#' @param colors_z The color scheme for nucleotides. Can be NULL (default ggseqlogo colors), a character string naming
#' an RColorBrewer palette (e.g., "Dark2"), or a named character vector (e.g., c("A"="blue", "C"="red", ...)).
#'
#' @return A 'ggplot' object.
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
#'
#' @export
plot_qqseqlogo_meme <- function(df_fragments,
                                end_motif_5p = "Fragment_Bases_5p",
                                end_motif_3p = "Fragment_Bases_3p",
                                motif_type = "Both",
                                motif_size = 3,
                                col_z = "Fragment_Status_Simple",
                                vals_z = NULL,
                                colors_z = NULL) {
    # --- 1. Input Validation ---
    if (is.null(col_z) && !is.null(vals_z)) stop("If 'col_z' is NULL, 'vals_z' must also be NULL.")
    if (!is.null(col_z) && !col_z %in% names(df_fragments)) stop(paste("Column", col_z, "not found in the dataframe."))
    if (!motif_type %in% c("Start", "End", "Both")) stop("motif_type must be one of 'Start', 'End', or 'Both'.")

    # --- 2. Data Filtering and Preparation ---
    df_filtered <- df_fragments
    if (!is.null(col_z)) {
        df_filtered <- df_filtered %>%
            filter(!is.na(.data[[col_z]]))
        if (!is.null(vals_z)) {
            df_filtered <- df_filtered %>%
                filter(.data[[col_z]] %in% vals_z)
        }
        df_filtered[[col_z]] <- factor(df_filtered[[col_z]], levels = if (!is.null(vals_z)) vals_z else unique(df_filtered[[col_z]]))
    } else {
        df_filtered$placeholder_group <- "All Fragments"
        col_z <- "placeholder_group"
    }

    # --- 3. Motif Extraction to a Named List ---
    process_group <- function(group_df) {
        group_name <- as.character(unique(group_df[[col_z]]))
        start_motifs <- character(0)
        end_motifs <- character(0)

        if (motif_type %in% c("Start", "Both")) {
            sequences <- na.omit(group_df[[end_motif_5p]])
            if (length(sequences) > 0) {
                min_len_5p <- min(nchar(sequences))
                if (min_len_5p < motif_size) {
                    warning(sprintf(
                        "For group '%s', 5' motif_size (%d) truncated to %d due to short sequences.",
                        group_name, motif_size, min_len_5p
                    ), call. = FALSE)
                }

                sequences <- sequences[nchar(sequences) >= motif_size]
                if (length(sequences) > 0) {
                    start_motifs <- substr(sequences, 1, motif_size)
                }
            }
        }

        if (motif_type %in% c("End", "Both")) {
            sequences <- na.omit(group_df[[end_motif_3p]])
            if (length(sequences) > 0) {
                min_len_3p <- min(nchar(sequences))
                if (min_len_3p < motif_size) {
                    warning(sprintf(
                        "For group '%s', 3' motif_size (%d) truncated to %d due to short sequences.",
                        group_name, motif_size, min_len_3p
                    ), call. = FALSE)
                }


                sequences <- sequences[nchar(sequences) >= motif_size]
                if (length(sequences) > 0) {
                    end_motifs <- str_sub(sequences, -motif_size, -1)
                }
            }
        }

        if (motif_type == "Both") {
            common_length <- min(length(start_motifs), length(end_motifs))
            return(paste0(start_motifs[1:common_length], "-", end_motifs[1:common_length]))
        } else if (motif_type == "Start") {
            return(start_motifs)
        } else {
            return(end_motifs)
        }
    }

    list_of_groups <- df_filtered %>% group_split(.data[[col_z]])
    list_of_motifs <- map(list_of_groups, process_group) %>%
        set_names(map_chr(list_of_groups, ~ paste0(unique(.x[[col_z]]), " (N=", nrow(.x), ")")))
    list_of_motifs <- list_of_motifs[map_int(list_of_motifs, length) > 0]

    # --- 4. Filter out motifs containing 'N' ---
    original_counts <- map_int(list_of_motifs, length)
    list_of_motifs <- map(list_of_motifs, ~ .x[!stringr::str_detect(.x, "N")])
    new_counts <- map_int(list_of_motifs, length)
    total_removed <- sum(original_counts - new_counts)
    if (total_removed > 0) {
        warning(paste(total_removed, "motifs containing 'N' were found and removed before plotting."))
    }
    list_of_motifs <- list_of_motifs[map_int(list_of_motifs, length) > 0]

    # --- 5. Color Scheme Setup ---
    if (length(list_of_motifs) == 0) {
        message("No data available to plot after filtering.")
        return(ggplot() +
            theme_void() +
            labs(title = "No data to display"))
    }
    color_scheme <- NULL
    if (!is.null(colors_z)) {
        unique_nucleotides <- list_of_motifs %>%
            unlist() %>%
            str_split("") %>%
            unlist() %>%
            unique() %>%
            sort()
        n_colors_needed <- length(unique_nucleotides)

        if (is.character(colors_z) && length(colors_z) == 1 && colors_z %in% rownames(brewer.pal.info)) {
            max_palette_colors <- brewer.pal.info[colors_z, "maxcolors"]
            if (n_colors_needed > max_palette_colors) {
                stop(paste0("Palette '", colors_z, "' only has ", max_palette_colors, " colors, but ", n_colors_needed, " are needed."))
            }
            palette_cols <- brewer.pal(n_colors_needed, colors_z)
            color_scheme <- make_col_scheme(chars = unique_nucleotides, cols = palette_cols)
        } else if (is.character(colors_z) && !is.null(names(colors_z))) {
            if (!all(unique_nucleotides %in% names(colors_z))) {
                missing_nucs <- setdiff(unique_nucleotides, names(colors_z))
                stop(paste("Custom colors provided, but missing definitions for nucleotide(s):", paste(missing_nucs, collapse = ", ")))
            }
            color_scheme <- make_col_scheme(chars = names(colors_z), cols = colors_z)
        } else {
            stop("'colors_z' must be NULL, a valid RColorBrewer palette name, or a named character vector.")
        }
    }

    # --- 6. Plotting & Theming ---
    if (length(list_of_motifs) == 0) {
        message("No data available to plot after filtering.")
        return(ggplot() +
            theme_void() +
            labs(title = "No data to display"))
    }

    # Create the base ggseqlogo plot
    p <- ggseqlogo(list_of_motifs, method = "prob", col_scheme = color_scheme, stack_width = 0.95, font = "helvetica_regular") +
        labs(
            title = "Sequence Motif Composition",
            y = "Frequency",
            x = "Position"
        )

    # Prepare custom axis labels
    motif_len <- nchar(list_of_motifs[[1]][1])
    custom_labels <- NULL
    if (motif_type == "Start") {
        custom_labels <- as.character(1:motif_len)
    } else if (motif_type == "End") {
        custom_labels <- as.character((-motif_len):-1)
    } else if (motif_type == "Both") {
        start_labels <- 1:motif_size
        end_labels <- (-motif_size):-1
        custom_labels <- as.character(c(start_labels, "", end_labels))
    }

    # Create a robust annotation dataframe for custom labels
    facet_names <- names(list_of_motifs)
    base_annotation_df <- tibble(
        x_pos = 1:motif_len,
        label = custom_labels
    )
    annotation_df <- tidyr::crossing(
        group = factor(facet_names, levels = facet_names),
        base_annotation_df
    )

    # Modify the existing coordinate system to allow drawing outside the panel
    p$coordinates$clip <- "off"

    p <- p +
        # Add custom axis labels using geom_text
        geom_text(
            data = annotation_df,
            aes(x = x_pos, y = -0.1, label = label),
            size = 5,
            inherit.aes = FALSE
        ) +
        # Apply custom theme
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            strip.text = element_text(face = "bold", size = 14),
            axis.title.x = element_text(size = 14, margin = margin(t = 15)),
            axis.title.y = element_text(size = 14),
            strip.background = element_rect(fill = "grey90", color = "grey90"),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            # Add margin at the bottom to make space for the labels
            plot.margin = margin(t = 5.5, r = 5.5, b = 20, l = 5.5, "pt")
        )

    # Add separator line for "Both" motif type
    if (motif_type == "Both") {
        p <- p + geom_vline(
            xintercept = motif_size + 1,
            linetype = "dashed",
            color = "grey40",
            linewidth = 0.8
        )
    }

    return(p)
}

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
#' @param sample_id Sample identifier.
#' @param output_path Character vector for the plot output path.
#' @param ggsave_params A named list of arguments to be passed to 'ggplot2::ggsave()'. For example,
#'   'list(width = 8, height = 6, units = "in", dpi = 300, bg = "white")'. If not provided, sensible defaults will be used.
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
#'
#' @examples
#' ## --- Create a dataset for demonstration ---
#' # Set a seed for reproducibility
#' set.seed(42)
#'
#' # Helper function to generate random DNA sequences with a bias
#' generate_biased_dna <- function(n_seq, len, prob) {
#'   bases <- c("A", "C", "G", "T")
#'   replicate(n_seq, paste(sample(bases, len, replace = TRUE, prob = prob), collapse = ""))
#' }
#'
#' # Create 50 "MUT" fragments with a high proportion of 'C' at the ends
#' df_mut <- data.frame(
#'   Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'   Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.2, 0.5, 0.15, 0.15)),
#'   Fragment_Status_Simple = "MUT"
#' )
#'
#' # Create 50 "WT" fragments with a high proportion of 'G' at the ends
#' df_wt <- data.frame(
#'   Fragment_Bases_5p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
#'   Fragment_Bases_3p = generate_biased_dna(50, 10, prob = c(0.15, 0.15, 0.5, 0.2)),
#'   Fragment_Status_Simple = "WT"
#' )
#'
#' # Combine into a single dataframe
#' example_df <- rbind(df_mut, df_wt)
#'
#' ## --- Function Calls ---
#'
#' # 1. Default plot: Shows a 3-base motif from both 5' and 3' ends,
#' #    separated by a dash, for each group ("MUT" and "WT").
#' p1 <- plot_qqseqlogo_meme(example_df)
#' print(p1)
#'
#' # 2. Analyze a longer, single-end motif: Shows a 5-base motif
#' #    from only the 5' end ('motif_type = "Start"').
#' p2 <- plot_qqseqlogo_meme(
#'   df_fragments = example_df,
#'   motif_type = "Start",
#'   motif_size = 5
#' )
#' print(p2)
#'
#' # 3. Customizing colors: Use a named RColorBrewer palette.
#' #    Note the separator "-" is not a nucleotide and won't be colored.
#' p3 <- plot_qqseqlogo_meme(
#'   df_fragments = example_df,
#'   colors_z = "Set1"
#' )
#' print(p3)
#'
#' # You can also provide a named vector for full control over colors.
#' custom_cols <- c("A" = "#1B9E77", "C" = "#D95F02", "G" = "#7570B3", "T" = "#E7298A")
#' p4 <- plot_qqseqlogo_meme(
#'   df_fragments = example_df,
#'   motif_type = "Start",
#'   colors_z = custom_cols
#' )
#' print(p4)
#'
#' # 4. Ungrouped plot: Analyzes all fragments together as a single group.
#' p5 <- plot_qqseqlogo_meme(example_df, col_z = NULL)
#' print(p5)
#'
#' # 5. Save plot with default settings.
#' # plot_qqseqlogo_meme(
#' #   df_fragments = example_df,
#' #   sample_id = "test01_motif",
#' #   output_path = "test01_motif_plot.png"
#' # )
#'
#' # 6. Save plot with custom dimensions.
#' # plot_qqseqlogo_meme(
#' #   df_fragments = example_df,
#' #   sample_id = "test02_motif_custom",
#' #   output_path = "test02_motif_custom_plot.png",
#' #   ggsave_params = list(width = 15, height = 10, units = "cm")
#' # )
#'
plot_qqseqlogo_meme <- function(df_fragments,
                                end_motif_5p = "Fragment_Bases_5p",
                                end_motif_3p = "Fragment_Bases_3p",
                                motif_type = "Both",
                                motif_size = 3,
                                col_z = "Fragment_Status_Simple",
                                vals_z = NULL,
                                colors_z = NULL,
                                sample_id = NA_character_,
                                output_path = NA_character_,
                                ggsave_params = list()) {
  # --- 1. Input Validation ---
  if (is.null(col_z) && !is.null(vals_z)) stop("If 'col_z' is NULL, 'vals_z' must also be NULL.")
  if (!is.null(col_z) && !col_z %in% names(df_fragments)) {
    stop(sprintf("Column '%s' not found in the dataframe.", col_z))
  }
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

  # Find the length of the shortest sequence in the relevant columns.
  min_len_available <- Inf

  # Check sequences if they are part of the analysis.
  if (motif_type %in% c("Start", "Both") && end_motif_5p %in% names(df_filtered)) {
    sequences_5p <- na.omit(df_filtered[[end_motif_5p]])
    if (length(sequences_5p) > 0) {
      min_len_available <- min(min_len_available, min(nchar(sequences_5p)))
    }
  }

  if (motif_type %in% c("End", "Both") && end_motif_3p %in% names(df_filtered)) {
    sequences_3p <- na.omit(df_filtered[[end_motif_3p]])
    if (length(sequences_3p) > 0) {
      min_len_available <- min(min_len_available, min(nchar(sequences_3p)))
    }
  }

  # If the requested motif_size is larger than the data supports, adjust it.
  if (motif_size > min_len_available) {
    # Warn the user that motif_size has been automatically reduced.
    warning(sprintf(
      "Requested 'motif_size' (%d) is larger than the shortest available sequence (%d). Size has been adjusted to %d.",
      motif_size, min_len_available, min_len_available
    ), call. = FALSE)

    # Cap the motif_size at the maximum possible value.
    motif_size <- min_len_available
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
      return(paste0(start_motifs[seq_len(common_length)], "-", end_motifs[seq_len(common_length)]))
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
    warning(sprintf(
      "%d motifs containing 'N' were found and removed before plotting.",
      total_removed
    ))
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
        stop(sprintf(
          "Palette '%s' only has %d colors, but %d are needed.",
          colors_z, max_palette_colors, n_colors_needed
        ))
      }
      palette_cols <- brewer.pal(n_colors_needed, colors_z)
      color_scheme <- make_col_scheme(chars = unique_nucleotides, cols = palette_cols)
    } else if (is.character(colors_z) && !is.null(names(colors_z))) {
      if (!all(unique_nucleotides %in% names(colors_z))) {
        missing_nucs <- setdiff(unique_nucleotides, names(colors_z))
        stop(sprintf(
          "Custom colors provided, but missing definitions for nucleotide(s): %s",
          paste(missing_nucs, collapse = ", ")
        ))
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
  final_plot <- ggseqlogo(list_of_motifs, method = "prob", col_scheme = color_scheme, stack_width = 0.95, font = "helvetica_regular") +
    labs(
      title = "Sequence Motif Composition",
      y = "Frequency",
      x = "Position"
    )

  # Prepare custom axis labels
  motif_len <- nchar(list_of_motifs[[1]][1])
  custom_labels <- NULL
  if (motif_type == "Start") {
    custom_labels <- as.character(seq_len(motif_len))
  } else if (motif_type == "End") {
    custom_labels <- as.character((-motif_len):-1)
  } else if (motif_type == "Both") {
    start_labels <- seq_len(motif_size)
    end_labels <- (-motif_size):-1
    custom_labels <- as.character(c(start_labels, "", end_labels))
  }

  # Create a robust annotation dataframe for custom labels
  facet_names <- names(list_of_motifs)
  base_annotation_df <- tibble(
    x_pos = seq_len(motif_len),
    label = custom_labels
  )
  annotation_df <- tidyr::crossing(
    group = factor(facet_names, levels = facet_names),
    base_annotation_df
  )

  # Modify the existing coordinate system to allow drawing outside the panel
  final_plot$coordinates$clip <- "off"

  final_plot <- final_plot +
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
    final_plot <- final_plot + geom_vline(
      xintercept = motif_size + 1,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.8
    )
  }

  # --- 7. Save the plot to a file if an output_path is provided ---
  if (is.null(output_path) ||
      (is.character(output_path) && length(output_path) == 1 &&
      (is.na(output_path) || !nzchar(output_path)))) {
    return(final_plot)
  }

  if (!(is.character(output_path) && length(output_path) == 1)) {
    stop("'output_path' must be a single character string.")
  }

  save_plot_if_needed(final_plot, output_path, ggsave_params)
  return(invisible(NULL))
}

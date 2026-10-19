#' Read genomic regions from a BED file
#'
#' @description
#' Reads a BED file containing genomic regions and converts its coordinates to
#' the 1-based, fully-closed convention used by BAM files.
#'
#' @details
#' BED files use **0-based, half-open** coordinates: the first base of a region
#' is `start`, while `end` is the position just after the last base. This
#' function converts them to **1-based, inclusive** coordinates with
#' `start = start + 1`, `end = end`, and `width = end - start + 1`.
#'
#' The first three columns must be `chrom`, `start`, and `end`. An optional
#' fourth column provides the region name; when it is absent or empty,
#' `region_id` is built as `chrom:start-end` using 1-based coordinates. Blank
#' lines, lines starting with `#`, and `track`/`browser` lines are ignored.
#'
#' @param region_path Path to a `.bed` or `.bed.gz` file.
#'
#' @return A `data.frame` with columns `region_id`, `chrom`, `start`, `end`, and
#'   `width`, using 1-based inclusive coordinates.
#'
#' @keywords internal
read_bed <- function(region_path) {
  if (length(region_path) != 1 || !is.character(region_path)) {
    stop("region_path must be a single character string.")
  }
  if (!grepl("\\.bed(\\.gz)?$", region_path, ignore.case = TRUE)) {
    stop("The region file must have a .bed or .bed.gz extension: ", region_path)
  }
  if (!file.exists(region_path)) {
    stop("The region file does not exist: ", region_path)
  }

  con <- if (grepl("\\.gz$", region_path, ignore.case = TRUE)) {
    gzfile(region_path, "r")
  } else {
    file(region_path, "r")
  }
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)

  skip <- grepl("^\\s*$", lines) |
    grepl("^#", lines) |
    grepl("^track(\\s|$)", lines, ignore.case = TRUE) |
    grepl("^browser(\\s|$)", lines, ignore.case = TRUE)

  valid_lines <- which(!skip)
  region_id <- character(length(valid_lines))
  chrom <- character(length(valid_lines))
  start <- integer(length(valid_lines))
  end <- integer(length(valid_lines))
  n_valid <- 0L

  for (i in valid_lines) {
    fields <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3) {
      stop(sprintf(
        "Line %d: a BED region must have at least 3 columns (chrom, start, end).",
        i
      ))
    }

    start_0 <- suppressWarnings(as.integer(fields[2]))
    end_0 <- suppressWarnings(as.integer(fields[3]))
    if (is.na(start_0)) {
      stop(sprintf("Line %d: 'start' must be an integer.", i))
    }
    if (is.na(end_0)) {
      stop(sprintf("Line %d: 'end' must be an integer.", i))
    }
    if (start_0 < 0) {
      stop(sprintf("Line %d: 'start' must be non-negative.", i))
    }
    if (end_0 <= start_0) {
      stop(sprintf("Line %d: 'end' must be greater than 'start'.", i))
    }

    region_start <- start_0 + 1L
    region_end <- end_0
    region_name <- if (length(fields) >= 4) trimws(fields[4]) else ""
    if (!nzchar(region_name)) {
      region_name <- paste0(fields[1], ":", region_start, "-", region_end)
    }

    n_valid <- n_valid + 1L
    region_id[n_valid] <- region_name
    chrom[n_valid] <- fields[1]
    start[n_valid] <- region_start
    end[n_valid] <- region_end
  }

  if (n_valid == 0L) {
    stop("No valid region found in the BED file: ", region_path)
  }

  region_id <- region_id[seq_len(n_valid)]
  chrom <- chrom[seq_len(n_valid)]
  start <- start[seq_len(n_valid)]
  end <- end[seq_len(n_valid)]

  data.frame(
    region_id = region_id,
    chrom = chrom,
    start = start,
    end = end,
    width = end - start + 1L,
    stringsAsFactors = FALSE
  )
}

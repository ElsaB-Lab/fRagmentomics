setup_region_test_files <- function() {
  list(
    region = system.file(
      "extdata/bed", "exemple-regions.bed",
      package = "fRagmentomics"
    ),
    bam = system.file(
      "extdata/bam", "cfdna-egfr-del_chr7_55241864_55243064_10k.bam",
      package = "fRagmentomics"
    )
  )
}

test_that("Region run returns the region schema without mutation columns", {
  files <- setup_region_test_files()

  res <- run_fRagmentomics_region(
    region = files$region,
    bam = files$bam,
    sample_id = "region_sample",
    report_bam_info = TRUE,
    report_softclip = TRUE,
    n_cores = 1L
  )

  expect_s4_class(res, "DataFrame")
  expect_gt(nrow(res), 0)

  expect_identical(
    names(res),
    c(
      "Sample_Id", "Chromosome", "Region_Id", "Region_Start", "Region_End",
      "Region_Width", "Fragment_Id", "Fragment_QC", "Fragment_Size",
      "Position_5p", "Position_3p", "POS_5p", "POS_3p", "FLAG_5p", "FLAG_3p",
      "MAPQ_5p", "MAPQ_3p", "CIGAR_5p", "CIGAR_3p", "TLEN",
      "Fragment_Bases_5p", "Fragment_Bases_3p", "Fragment_Basqs_5p",
      "Fragment_Basqs_3p", "Nb_Fragment_Bases_Softclip_5p",
      "Nb_Fragment_Bases_Softclip_3p"
    )
  )

  mutation_columns <- c(
    "Position", "Ref", "Alt", "Input_Mutation", "VAF",
    "Fragment_Status_Simple", "Fragment_Status_Detail", "Read_5p_Status",
    "Read_3p_Status", "BASE_5p", "BASE_3p", "BASQ_5p", "BASQ_3p"
  )
  expect_false(any(mutation_columns %in% names(res)))
})

test_that("Region run writes output_path, creates its folder, and returns NULL", {
  files <- setup_region_test_files()

  nested_dir <- file.path(tempdir(), "frag_test_region", "subdir")
  out_file <- file.path(nested_dir, "output.tsv")
  if (dir.exists(nested_dir)) unlink(nested_dir, recursive = TRUE, force = TRUE)

  res <- run_fRagmentomics_region(
    region = files$region,
    bam = files$bam,
    output_path = out_file,
    verbose = TRUE,
    n_cores = 1L
  )

  expect_null(res)
  expect_true(file.exists(out_file))

  res_df <- utils::read.table(
    out_file,
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )
  expect_gt(nrow(res_df), 0)
  expect_true(all(c("Region_Id", "Region_Start", "Region_End", "Region_Width") %in%
    names(res_df)))

  unlink(file.path(tempdir(), "frag_test_region"), recursive = TRUE, force = TRUE)
})

test_that("Region run warns once for regions wider than 10000 bp without erroring", {
  files <- setup_region_test_files()

  wide_bed <- tempfile(fileext = ".bed")
  writeLines(
    c(
      "chr7\t10600\t10900\tnormal-region",
      "chr7\t0\t20100\twide-region"
    ),
    wide_bed
  )
  on.exit(unlink(wide_bed), add = TRUE)

  expect_warning(
    res <- run_fRagmentomics_region(
      region = wide_bed,
      bam = files$bam,
      n_cores = 1L
    ),
    "wider than 10000"
  )
  expect_gt(nrow(res), 0)
})

test_that("Region run reports counters in verbose mode", {
  files <- setup_region_test_files()

  msgs <- testthat::capture_messages(
    res <- run_fRagmentomics_region(
      region = files$region,
      bam = files$bam,
      verbose = TRUE,
      n_cores = 1L
    )
  )

  expect_s4_class(res, "DataFrame")
  expect_true(any(grepl("Processing [0-9]+ region", msgs)))
  expect_true(any(grepl("Processed [0-9]+ fragment", msgs)))
})

test_that("Region run warns and skips a region without reads", {
  files <- setup_region_test_files()

  bed_path <- tempfile(fileext = ".bed")
  writeLines(
    c(
      "chr7\t10600\t10900\tnormal-region",
      "chr7\t100000\t100100\tempty-region"
    ),
    bed_path
  )
  on.exit(unlink(bed_path), add = TRUE)

  expect_warning(
    res <- run_fRagmentomics_region(
      region = bed_path,
      bam = files$bam,
      n_cores = 1L
    ),
    "No read covers the region"
  )
  expect_gt(nrow(res), 0)
})

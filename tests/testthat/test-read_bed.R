test_that("read_bed converts BED 0-based half-open coordinates to 1-based inclusive", {
  bed_path <- tempfile(fileext = ".bed")
  writeLines("chr1\t100\t200\tregion_a", bed_path)
  on.exit(unlink(bed_path), add = TRUE)

  regions <- read_bed(bed_path)

  expect_equal(nrow(regions), 1)
  expect_equal(regions$region_id, "region_a")
  expect_equal(regions$chrom, "chr1")
  expect_equal(regions$start, 101L)
  expect_equal(regions$end, 200L)
  expect_equal(regions$width, 100L)
})

test_that("read_bed falls back to chrom:start-end when the name is absent or empty", {
  bed_path <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t100\t200", "chr2\t10\t20\t"), bed_path)
  on.exit(unlink(bed_path), add = TRUE)

  regions <- read_bed(bed_path)

  expect_equal(regions$region_id, c("chr1:101-200", "chr2:11-20"))
})

test_that("read_bed ignores comments, track/browser lines and blank lines", {
  bed_path <- tempfile(fileext = ".bed")
  writeLines(
    c(
      "# a comment",
      "track name=example",
      "browser position chr1:1-10",
      "",
      "chr1\t10\t20\tfirst",
      "chr2\t30\t40\tsecond"
    ),
    bed_path
  )
  on.exit(unlink(bed_path), add = TRUE)

  regions <- read_bed(bed_path)

  expect_equal(nrow(regions), 2)
  expect_equal(regions$region_id, c("first", "second"))
})

test_that("read_bed reads .bed.gz files", {
  bed_path <- tempfile(fileext = ".bed.gz")
  con <- gzfile(bed_path, "w")
  writeLines("chr1\t100\t200\tcompressed", con)
  close(con)
  on.exit(unlink(bed_path), add = TRUE)

  regions <- read_bed(bed_path)

  expect_equal(regions$region_id, "compressed")
  expect_equal(regions$start, 101L)
  expect_equal(regions$end, 200L)
})

test_that("read_bed reports invalid inputs with explicit errors", {
  expect_error(
    read_bed(tempfile(fileext = ".bed")),
    "does not exist"
  )

  txt_path <- tempfile(fileext = ".txt")
  writeLines("chr1\t100\t200", txt_path)
  on.exit(unlink(txt_path), add = TRUE)
  expect_error(read_bed(txt_path), "extension")

  bed_path <- tempfile(fileext = ".bed")
  writeLines("chr1\t100", bed_path)
  on.exit(unlink(bed_path), add = TRUE)
  expect_error(read_bed(bed_path), "at least 3 columns")

  writeLines("chr1\t-1\t10", bed_path)
  expect_error(read_bed(bed_path), "start")

  writeLines("chr1\t100\t100", bed_path)
  expect_error(read_bed(bed_path), "end")

  writeLines(c("# only comments", "track name=x"), bed_path)
  expect_error(read_bed(bed_path), "No valid region")
})

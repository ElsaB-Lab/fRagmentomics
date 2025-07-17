test_that("Test function remove softclip", {
    # Test 1: Correctly removes soft clips from both ends
    read_stats_input <- list(
        CIGAR = "10S80M10S",
        SEQ = paste0(
            paste(rep("A", 10), collapse = ""),
            paste(rep("G", 80), collapse = ""),
            paste(rep("C", 10), collapse = "")
        ),
        QUAL = paste(rep("F", 100), collapse = "")
    )
    result <- remove_softclip(read_stats_input)

    # Assert: Define expected output and check if the result matches
    expected_output <- list(
        SEQ   = paste(rep("G", 80), collapse = ""),
        QUAL  = substr(read_stats_input$QUAL, 11, 90),
        CIGAR = "80M"
    )

    expect_equal(result, expected_output)

    # Test 2: Correctly removes 5p soft clips only
    read_stats_input <- list(
        CIGAR = "5S95M",
        SEQ = paste0(
            paste(rep("T", 5), collapse = ""),
            paste(rep("G", 95), collapse = "")
        ),
        QUAL = paste(rep("F", 100), collapse = "")
    )
    result <- remove_softclip(read_stats_input)

    # Assert
    expected_output <- list(
        SEQ   = paste(rep("G", 95), collapse = ""),
        QUAL  = substr(read_stats_input$QUAL, 6, 100),
        CIGAR = "95M"
    )

    expect_equal(result, expected_output)

    # Test 3: Correctly removes 3p soft clips only
    read_stats_input <- list(
        CIGAR = "92M8S",
        SEQ = paste0(
            paste(rep("G", 92), collapse = ""),
            paste(rep("C", 8), collapse = "")
        ),
        QUAL = paste(rep("F", 100), collapse = "")
    )
    result <- remove_softclip(read_stats_input)

    # Assert
    expected_output <- list(
        SEQ   = paste(rep("G", 92), collapse = ""),
        QUAL  = substr(read_stats_input$QUAL, 1, 92),
        CIGAR = "92M"
    )

    expect_equal(result, expected_output)

    # Test 4: Returns the original data when there are no soft clips
    read_stats_input <- list(
        SEQ   = paste(rep("G", 100), collapse = ""),
        QUAL  = paste(rep("F", 100), collapse = ""),
        CIGAR = "100M"
    )
    result <- remove_softclip(read_stats_input)

    # Assert: The result should be identical to the input
    expect_equal(result, read_stats_input)

    # Test 5: Returns empty strings when the entire read is soft-clipped
    read_stats_input <- list(
        CIGAR = "150S",
        SEQ   = paste(rep("N", 150), collapse = ""),
        QUAL  = paste(rep("#", 150), collapse = "")
    )
    result <- remove_softclip(read_stats_input)

    # Assert: All fields should be empty strings
    expected_output <- list(
        SEQ   = "",
        QUAL  = "",
        CIGAR = ""
    )

    expect_equal(result, expected_output)


    # Test 6: Correctly handles complex CIGARs with insertions/deletions
    read_stats_input <- list(
        CIGAR = "10S40M5I45M5S",
        SEQ   = paste(rep("A", 100), collapse = ""),
        QUAL  = paste(rep("F", 100), collapse = "")
    )
    result <- remove_softclip(read_stats_input)

    # Assert
    expected_output <- list(
        SEQ   = substr(read_stats_input$SEQ, 11, 95), # 100 - 10 - 5 = 85 chars long
        QUAL  = substr(read_stats_input$QUAL, 11, 95),
        CIGAR = "40M5I45M"
    )

    expect_equal(result, expected_output)

    # Test 7: Correctly ignores hard clips and only removes soft clips
    read_stats_input <- list(
        CIGAR = "5H10S80M5S5H",
        SEQ   = paste(rep("A", 95), collapse = ""),
        QUAL  = paste(rep("F", 95), collapse = "")
    )
    result <- remove_softclip(read_stats_input)

    # Assert: The 'S' parts should be gone, but 'H' should remain
    expected_output <- list(
        SEQ   = substr(read_stats_input$SEQ, 11, 90), # 95 - 10 - 5 = 80 chars
        QUAL  = substr(read_stats_input$QUAL, 11, 90),
        CIGAR = "5H80M5H"
    )

    expect_equal(result, expected_output)
})

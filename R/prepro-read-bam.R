preprocess_bam <- function(
    bam_file,
    chr,
    pos,
    neg_offset,
    pos_offset,
    flag_keep,
    flag_remove) 
    {
        # Convert into hexadecimal number
        flag_keep_str <- sprintf("0x%X", flag_keep)
        flag_remove_str <- sprintf("0x%X", flag_remove)

        # Extract reads around the interested position 
        # sprintf formats a string (%s replaced a string and %d an int)
        cmd_pos <- sprintf(
            "samtools view -b %s %s:%d-%d | samtools view -f %s -F %s - | cut -f1-11",
            bam_file, chr, pos, pos, flag_keep_str, flag_remove_str
        )
        
        # System to do a command in the terminal. Intern to have the result in R
        # The output is a vector of strings
        sam_pos_lines <- system(cmd_pos, intern = TRUE)
        
        # Conversion into a df
        sam_pos_df <- read.table(textConnection(sam_pos_lines),
                                sep = "\t", 
                                header = FALSE,
                                stringsAsFactors = FALSE)
        

        # Extract wider reads around the interested position 
        # Important to have the read who doesn't carry the mutation
        cmd_ext <- sprintf(
            "samtools view -b %s %s:%d-%d | samtools view -f %s -F %s - | cut -f1-11",
            bam_file, chr, neg_offset, pos_offset, flag_keep_str, flag_remove_str
        )
        
        sam_ext_lines <- system(cmd_ext, intern = TRUE)
        
        sam_ext_df <- read.table(textConnection(sam_ext_lines),
                                sep = "\t", 
                                header = FALSE,
                                stringsAsFactors = FALSE)
        
        # Select only the lines where QNAME (col 1) is in sam_pos_df
        # fragments_of_interest is a list of the QNAMEs in sam_pos_df
        fragments_of_interest <- unique(sam_pos_df[[1]])
        
        # We filter sam_ext_df 
        sam_final_df <- sam_ext_df[sam_ext_df[[1]] %in% fragments_of_interest, ]
        
        # Return df with sam data
        return(sam_final_df)
}

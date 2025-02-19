# @created: 11 Oct 24
# @modified: 11 Oct 24
# @authors: Juliette Samaniego
#
# Run R package

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

# set command line arguments =============================================================================================================
args <- commandArgs(trailingOnly = TRUE)

# stop the script if no command line argument
if (length(args)==0) {
  stop("Requires command line argument.")
} else if (length(args)==1) {
  args[2] = "{wildcards.sample}"
  args[3] = "{wildcards.alt}"
  args[4] = "{output.tsv}"
}


# function ===============================================================================================================================
# function takes a sample tsv file, creates new columns in the dataframe, removes duplicates and saves a new tsv file with all changes
RemoveDuplicates <- function(inputfile, sample, alternatebase,  outputfile) {
    Base_read_mut_mut <- paste(alternatebase, alternatebase, sep = "_")
    Base_read_mut_NA <- paste(alternatebase, "NA", sep = "_")
    Base_read_NA_mut <- paste("NA", alternatebase, sep = "_")
        
    df_file <- lapply(inputfile, function(inputfile) {
    df <- read_tsv(inputfile, show_col_types=F, progress=F, col_types = cols(Base_Read_1 = "c", Base_Read_2 = "c"), quote = "")
    df[, "Samples_ID"] <- sample
    

    # Creating a column which gives the identity of the bases of each read linked by "_" from the reading and position of interest
    df$Base_Read <- paste(df$Base_Read_1, df$Base_Read_2, sep="_")

    # Creation of a "Mutation_Status" column to define whether the fragment is mutated "MUT" or wild-type "WT"
    df$Mutation_Status <- sapply(1:nrow(df), function(i) {
    Status <- "WT"
    x <- df[i,]
    if (!is.na(x$Base_Read)) {
        if (x$Base_Read == Base_read_mut_mut) {
           Status <- "MUT"
            }
        }
    if (!is.na(x$Base_Read)) {
        if (x$Base_Read == Base_read_mut_NA) {
           Status <- "MUT"
            }
        }
    if (!is.na(x$Base_Read)) {
        if (x$Base_Read == Base_read_NA_mut) {
           Status <- "MUT"
            }
        }
     return(Status) } )

    # Creation of a column "last_4_characters" which retrieves the last 4 characters of the name given to the fragment
    # corresponding to a part of the UMI
    n_last <- 4
    last_4_characters <- substr(df$Fragment, nchar(df$Fragment) - n_last + 1, nchar(df$Fragment))
    df$last_4_characters <- substr(df$Fragment, nchar(df$Fragment) - n_last + 1, nchar(df$Fragment))

    # Creation of a "frag_id" column which gathers the last 4 characters of the name given to the fragment
    # and the starting positions of each of the reads (read 1 & 2)
    df$frag_id <- paste(df$last_4_characters,df$Startpos_read_1,df$Startpos_read_2, sep="_")

    # Create a column "last_4_characters_inverted" which inverts the last 4 chercaters of the name given to the fragment to obtain XXYY -> YYXX
    df$last_4_characters_inverted <- paste0(substr(x = df$last_4_characters, start = 3, stop = 4),
                                            substr(x = df$last_4_characters, start = 1, stop = 2))

    # Creation of column "UMI_characters" grouping the last 4 charcaters and the last 4 inverted separated by "_" definig the UMI given by FMI
    df$UMI_characters <- paste0(pmin(df$last_4_characters, df$last_4_characters_inverted), '_' ,
                                pmax(df$last_4_characters, df$last_4_characters_inverted))
    # Create two column linking the starting positions of the reads in one direction and the other
    df$Startpos1_Startpos2 <- paste(df$Startpos_read_1, df$Startpos_read_2, sep="_")
    df$Startpos2_Startpos1 <- paste(df$Startpos_read_2, df$Startpos_read_1, sep="_")

    # Creation of "Startpos_read" column taking the larger value between Startpos1_Startpos2 and Startpos2_Startpos1,
    df$Startpos_read <- pmax(df$Startpos1_Startpos2, df$Startpos2_Startpos1)

    # Creation of a "frag_id_2" column linking UMI_characters and Startpos_read
    df$frag_id_2 <- paste(df$UMI_characters, df$Startpos_read, sep="_")

    # Creation of a "Start_test" column which will give the starting position by removing the number of soft-clipped bases on the left
    df$Start_test <- sapply(1:nrow(df), function(i) {
    Status <- 0
    x <- df[i,]
    if (x$Startpos_read_1 <= x$Startpos_read_2){
      Status <- x$Startpos_read_1 - x$Bases_soft_clip_left
    } else {
      Status <- x$Startpos_read_2 - x$Bases_soft_clip_left
    }
    return(Status) } )

    # Creation of a "Start_test_read_1" column which will give the starting position by removing the number of soft-clipped bases
    # on the left if Startpos_read_1 <= Startpos_read_2, otherwise it will be equal to 0
    df$Start_test_read_1 <- sapply(1:nrow(df), function(i) {
    Status <- 0
    x <- df[i,]
    if (x$Startpos_read_1 <= x$Startpos_read_2){
      Status <- x$Startpos_read_1 - x$Bases_soft_clip_left
    }
    return(Status) } )
    # Creation of a "Start_test_read_2" column which will give the starting position by removing the number of soft-clipped bases
    # on the left if Startpos_read_1 >= Startpos_read_2, otherwise it will be equal to 0
    df$Start_test_read_2 <- sapply(1:nrow(df), function(i) {
    Status <- 0
    x <- df[i,]
    if (x$Startpos_read_1 >= x$Startpos_read_2){
      Status <- x$Startpos_read_2 - x$Bases_soft_clip_left
    }
    return(Status) } )

    # Create "Startpos_read_1_cor" column to determine correct starting position
    df$Startpos_read_1_cor <- sapply(1:nrow(df), function(i) {
    Status <- 0
    x <- df[i,]
    if (x$Start_test_read_1 == 0){
      Status <- x$Startpos_read_1
    }
    if (x$Start_test_read_1 != 0){
      Status <- x$Start_test_read_1
    }
    return(Status) } )
    #Create "Startpos_read_2_cor" column to determine correct starting position
    df$Startpos_read_2_cor <- sapply(1:nrow(df), function(i) {
    Status <- 0
    x <- df[i,]
    if (x$Start_test_read_2 == 0){
      Status <- x$Startpos_read_2
    }
    if (x$Start_test_read_2 != 0){
      Status <- x$Start_test_read_2
    }
    return(Status) } )

    # Create two column linking the correct starting positions of the reads in one direction and the other
    df$Startpos1_Startpos2_cor <- paste(df$Startpos_read_1_cor, df$Startpos_read_2_cor, sep="_")
    df$Startpos2_Startpos1_cor <- paste(df$Startpos_read_2_cor, df$Startpos_read_1_cor, sep="_")
    # Creation of "Startpos_read_cor" column taking the larger value between Startpos1_Startpos2_cor and Startpos2_Startpos1_cor,
    # this column will then be used to delete duplicates
    df$Startpos_read_cor <- pmax(df$Startpos1_Startpos2_cor, df$Startpos2_Startpos1_cor)

    # Use of duplicated function to delete duplicate lines based on the UMI_characters, Startpos_read_cor, Mutation_Status, Absolute_size columns 
    df <- df %>%
    filter(!duplicated(cbind(UMI_characters, Startpos_read_cor, Mutation_Status, Absolute_size)))

return(df)})
    
    df_file <- do.call(rbind, df_file)
    write_tsv(df_file, file=outputfile)

}

# run ==============================================================================================================================================

RemoveDuplicates(inputfile = args[1], sample = args[2], alternatebase = args[3], outputfile = args[4])



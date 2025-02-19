# Project : ElsaBLab_fRagmentomics

#' Prepare R sockets for running in parallel.
#'
#' @inheritParams build_fragments_info_table
#' @return an object of class "cluster"
#'
#' @importFrom parallel makePSOCKcluster clusterExport clusterEvalQ
#' @importFrom doParallel registerDoParallel
#'
#' @noRd
setup_parallel_computations <- function(n_cores){
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  # tmp <- clusterEvalQ(cl, library(fRagmentomics))
  cl
}


#' Create a table with statistics about fragments covering a position of interest (e.g the position of a mutation).
#'
#' @param df_sam dataframe containing reads from SAM format
#' @param chr chromosome for position of interest
#' @param pos chromosome for position of interest
#' @param ref reference sequence for position of interest
#' @param alt alternative sequence for position of interest
#' @param del_info character vector in the format "$pos_final,$rep_del"
#' @param mutation_type type of the mutation in insertion, deletion and mutation
#' @param n_cores the number of cores available
#' @return a dataframe
#'
#' @importFrom readr read_tsv
#' @importFrom stringr str_extract_all
#' @importFrom parallel stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom utils txtProgressBar setTxtProgressBar flush.console
#' @export
build_fragments_info_table <- function(df_sam, sample_id, chr, pos, ref, alt, del_info, mutation_type, n_cores=1){

  # list of all fragments names
  fragments_names <- unique(df_sam[,1,drop=T])
  n_fragments <- length(fragments_names)

  # prepare cluster
  cat(paste("-exporting functions to cluster..."))
  cl <- setup_parallel_computations(n_cores)
  cat("done!\n")

  # function to rbind with progressbar
  rbindprogress <- function(iterator){
    pb <- txtProgressBar(min = 1, max = iterator - 1, style = 3)
    count <- 0
    function(...) {
      count <<- count + length(list(...)) - 1
      setTxtProgressBar(pb, count)
      flush.console()
      rbind(...)
    }
  }

  # run a loop on fragments names
  i <- NULL
  df_fragments_info <- foreach(i=1:n_fragments,
                               .combine = rbindprogress(n_fragments),
                               .inorder = F,
                               .multicombine = F) %dopar% {
    fragment_name <- fragments_names[[i]]

    # select reads corresponding to fragment
    df_fragment_reads <- df_sam[df_sam[,1,drop=T]==fragment_name,,drop=F]

    # sanity checks
    fragment_qc <- ""

    if (!all(df_fragment_reads[,3]==chr)){
      fragment_qc <- paste(fragment_qc, "&", paste("read(s) from other chromosome than", chr))
    }

    if (nrow(df_fragment_reads)!=2){
      fragment_qc <- paste(fragment_qc, "&", paste(nrow(df_fragment_reads), "read(s)"))
    }

    if (fragment_qc != ""){
      # fragments does not pass checks
      df_fragment <- data.frame(Chromosome=chr, Position=pos, Fragment=fragment_name, Fragment_Check=fragment_qc)
    } else {
      ####
      #### load reads 1 and 2 into one-line dataframes
      ####

      # read 1 is defined as the read having 0x40 in its flag while read 2 has 0x80 in its flag value
      flag_a <- df_fragment_reads[1,"flag",drop=T]
      flag_b <- df_fragment_reads[2,"flag",drop=T]

      flag_a <- bitvalues_from_bam_flag(as.integer(flag_a), bitnames=c("isFirstMateRead","isSecondMateRead"))
      flag_b <- bitvalues_from_bam_flag(as.integer(flag_b), bitnames=c("isFirstMateRead","isSecondMateRead"))

      if ( (flag_a[,"isFirstMateRead"] == 1) & (flag_b[,"isSecondMateRead"] == 1) ){
        df_read1 <- df_fragment_reads[1, ]
        df_read2 <- df_fragment_reads[2, ]
      } else if ( (flag_b[,"isFirstMateRead"] == 1) & (flag_a[,"isSecondMateRead"] == 1) ){
        df_read1 <- df_fragment_reads[2, ]
        df_read2 <- df_fragment_reads[1, ]
      } else {
        stop("invalid read flags. One read should be first mate, the other second mate:\n", flag_a, "\n", flag_b)
      }

      ####
      #### retrieve read-level stats
      ####

      # mapping quality
      mapq1 <- df_read1$mapq
      mapq2 <- df_read2$mapq

      # tlen
      tlen1 <- df_read1$tlen
      tlen2 <- df_read2$tlen

      # reads cigar
      r_cigar1 <- df_read1$cigar
      r_cigar2 <- df_read2$cigar

      # reads mapped positions
      r_pos1 <- df_read1$pos
      r_pos2 <- df_read2$pos

      # reads queries
      r_query1 <- df_read1$seq
      r_query2 <- df_read2$seq

      # reads qualities
      r_qual1 <- df_read1$qual
      r_qual2 <- df_read2$qual

      # reads length
      read_length1 <- nchar(r_query1)
      read_length2 <- nchar(r_query2)

      # identity variable to identify the position of read 1 and 2 depending on whether they are on the left or right
      if (r_pos1 <= r_pos2){
        identity_RIGHT <- 2
        identity_LEFT <- 1
      } else {
        identity_RIGHT <- 1
        identity_LEFT <- 2
      }

      ####
      #### retrieve base and qualities
      ####

      if (mutation_type == "mutation") {
        r_info1 <- get_base_qual_from_read(pos=pos, r_pos=r_pos1, r_cigar=r_cigar1, r_query=r_query1, r_qual=r_qual1)
        r_info2 <- get_base_qual_from_read(pos=pos, r_pos=r_pos2, r_cigar=r_cigar2, r_query=r_query2, r_qual=r_qual2)

      } else if (mutation_type == "deletion") {
        r_info1 <- get_deletion(pos=pos, ref=ref, del_info=del_info, r_pos=r_pos1, r_cigar=r_cigar1)
        r_info2 <- get_deletion(pos=pos, ref=ref, del_info=del_info, r_pos=r_pos2, r_cigar=r_cigar2)

      } else if (mutation_type == "insertion") {
        r_info1 <- get_insertion(pos=pos, alt=alt, r_pos=r_pos1, r_cigar=r_cigar1, r_query=r_query1, r_qual=r_qual1)
        r_info2 <- get_insertion(pos=pos, alt=alt, r_pos=r_pos2, r_cigar=r_cigar2, r_query=r_query2, r_qual=r_qual2)

      } else {
        r_info1 <- "Mutation_type is not defined. You should define it as 'mutation', 'deletion' or 'insertion'"
        r_info2 <- "Mutation_type is not defined. You should define it as 'mutation', 'deletion' or 'insertion'"
      }


      ####
      #### retrieve bases : inserted, deleted, soft-clipped left and soft-clipped right
      ####

      # bases soft-clipped left read 1/2
      bases_softcl_left1 <- as.numeric(str_extract(r_cigar1, "^[[:digit:]]+(?=S)"))
      bases_softcl_left1 <- ifelse(is.na(bases_softcl_left1), 0, bases_softcl_left1)
      bases_softcl_left2 <- as.numeric(str_extract(r_cigar2, "^[[:digit:]]+(?=S)"))
      bases_softcl_left2 <- ifelse(is.na(bases_softcl_left2), 0, bases_softcl_left2)

      # bases soft-clipped right read 1/2
      bases_softcl_right1 <- as.numeric(str_extract(r_cigar1, "[[:digit:]]+(?=S$)"))
      bases_softcl_right1 <- ifelse(is.na(bases_softcl_right1), 0, bases_softcl_right1)
      bases_softcl_right2 <- as.numeric(str_extract(r_cigar2, "[[:digit:]]+(?=S$)"))
      bases_softcl_right2 <- ifelse(is.na(bases_softcl_right2), 0, bases_softcl_right2)

      # number of soft clipped bases genomically left
      if (r_pos1 <= r_pos2){
        bases_softcl_left <- bases_softcl_left1
        bases_softcl_right <- bases_softcl_right2
      } else {
        bases_softcl_left <- bases_softcl_left2
        bases_softcl_right <- bases_softcl_right1
      }

      # counts of unique deletions and insertions
      n_unique_del_ins <- get_nb_unique_del_ins(r_pos1, r_cigar1, r_pos2, r_cigar2)
      n_unique_del <- n_unique_del_ins$unique_del
      n_unique_ins <- n_unique_del_ins$unique_ins

      ####
      #### retrieve inner distance and absolute_size
      ####
      inner_distances <- get_fragment_inner_distance(r_pos1, r_cigar1, read_length1, r_pos2, r_cigar2, read_length2)
      inner_distance_a <- inner_distances$a
      inner_distance_b <- inner_distances$b

      # # Check that the two formulas for calculating inner distance return the same values
      # if (inner_distance_a == inner_distance_b){
      #   inner_distance <- inner_distance_a
      # } else {
      #   inner_distance <- "error"
      # }

      # absolute size
      absolute_size <- read_length1 + read_length2 + inner_distance_b

      # result dataframe
      df_fragment <- data.frame(Sample_Id=sample_id, Chromosome=chr, Position=pos, Ref=ref, Alt=alt, Fragment=fragment_name, Fragment_Check="OK",
                                Mutation_type=mutation_type, Mapq_Read_1=mapq1, Mapq_Read_2=mapq2, Base_Read_1=r_info1$base,
                                Base_Read_2=r_info2$base, Qual_Read_1=r_info1$qual, Qual_Read_2=r_info2$qual,
                                Indel_1=r_info1$indel, Indel_2=r_info2$indel,
                                TLEN=abs(tlen1), Read_length_1=read_length1, Read_length_2=read_length2,
                                Startpos_read_1=r_pos1, Startpos_read_2=r_pos2,
                                Identity_RIGHT=identity_RIGHT, Identity_LEFT=identity_LEFT,
                                Bases_del=n_unique_del, Bases_ins=n_unique_ins,
                                Bases_soft_clip_left=bases_softcl_left, Bases_soft_clip_right=bases_softcl_right,
                                Inner_distance_a=inner_distance_a, Inner_distance_b=inner_distance_b,
                                Absolute_size=absolute_size)
    }
    df_fragment
  }

  # terminate cluster
  stopCluster(cl)

  df_fragments_info
}

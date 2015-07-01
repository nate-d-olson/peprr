## Code to generate figures and tables for ROA
#### Summary of sequencing datasets --------------------------------------------
#' create df for dataset summary table
#' @param db_con peprDB connection
#' @return NULL
seq_summary_table <- function(db_con){
    # table with number of reads, read length, coverage, ect.
    seq_metrics <- dplyr::tbl(src = db_con, from="align")  %>%
        dplyr::select(accession, CATEGORY, TOTAL_READS, MEAN_READ_LENGTH)  %>%
        dplyr::filter(CATEGORY %in% c("UNPAIRED","PAIR"))  %>%
        dplyr::collect() %>%
        dplyr::mutate(accession = gsub(pattern = ".+_",
                                       replacement = "",
                                       x = accession))


    insert_tbl <- dplyr::tbl(src = db_con, from="insert_hist")  %>%
        dplyr::group_by(accession)  %>%
        dplyr::summarize(mean_insert =
            sum(insert_size *All_Reads.fr_count)/sum(All_Reads.fr_count))  %>%
        dplyr::collect()  %>%
        dplyr::mutate(accession = gsub(pattern = ".+_",
                                       replacement = "",
                                       x = accession))

    coverage_tbl <- dplyr::tbl(src = db_con, from="coverage") %>%
        dplyr::rename("accession"=SAMPLE) %>%
        dplyr::collect()

    dplyr::tbl(src = db_con, from ="exp_design")  %>%
        dplyr::collect() %>% dplyr::full_join(seq_metrics) %>%
        dplyr::full_join(insert_tbl) %>%
        dplyr::full_join(coverage_tbl) %>%
        dplyr::select(-CATEGORY) %>%
        dplyr::rename("Accession" = accession, "Platform" = plat,
                      "Vial" = vial, "Replicate" = rep,
                      "Reads"= TOTAL_READS, "Read Length" = MEAN_READ_LENGTH,
                      "Insert Size"=mean_insert, "Coverage"= COV)
}

#### Pilon Changes -------------------------------------------------------------
#' create df for pilon changes
#' @param db_con peprDB connection
#' @return NULL
pilon_changes_table <- function(db_con){
    dplyr::tbl(src = db_con, from="pilon_changes")  %>%
        dplyr::collect()
}


## Code to generate figures and tables for ROA
#### Summary of sequencing datasets --------------------------------------------
#' create df for dataset summary table
#' @param db_con peprDB connection
#' @return NULL
seq_summary_table <- function(db_con){
    # table with number of reads, read length, coverage, ect.
    seq_metrics <- dplyr::tbl(src = peprDB, from="align")  %>%
        dplyr::select(accession, CATEGORY, TOTAL_READS, MEAN_READ_LENGTH)  %>%
        dplyr::filter(CATEGORY %in% c("UNPAIRED","PAIR"))  %>%
        dplyr::collect()  %>%
        tidyr::separate(accession, c("ref","pilon", "accession"), sep = "_")

    insert_tbl <- dplyr::tbl(src = peprDB, from="insert_hist")  %>%
        dplyr::group_by(accession)  %>%
        dplyr::summarize(mean_insert =
            sum(insert_size *All_Reads.fr_count)/sum(All_Reads.fr_count))  %>%
        dplyr::collect()  %>%
        tidyr::separate(accession, c("ref","pilon", "accession"), sep = "_")

    dplyr::tbl(src = peprDB, from ="exp_design")  %>%
        dplyr::collect() %>% dplyr::full_join(seq_metrics) %>%
        dplyr::full_join(insert_tbl) %>%
        dplyr::select(-ref, -pilon,-CATEGORY)
}

#### Pilon Changes -------------------------------------------------------------
#' create df for pilon changes
#' @param db_con peprDB connection
#' @return NULL
pilon_changes_table <- function(db_con){
    dplyr::tbl(src = peprDB, from="pilon_changes")  %>%
        dplyr::collect()
}
# table of pilon results

##### Base level analysis ####
## purity plots - script to generate purity plots


## Analysis of low purity positions
# figures correlating different parameters with low purity positions - e.g. bias and coverage

# table of low purity positions

## indel analysis
# correlate with homopolymer annotation results

#### Homogeneity ####
# table, figure, depending on results


### Genomic Purity -------------------------------------------------------------
.genomic_purity_df <- function(db_con){
    pathoDF <- read.csv("~/Desktop/micro_rm/micro_rm_dev/analysis/stats/genomic_purity/pathoDF.csv",
                        stringsAsFactors = FALSE) %>%
        select(-X)
    sampleDF <- read.csv("~/Desktop/micro_rm/micro_rm_dev/analysis/stats/genomic_purity/sampleDF.csv",
                         stringsAsFactors = FALSE) %>%
        select(-X)

    df <- left_join(pathoDF, sampleDF)
    df <- tbl(src = peprDB, from = "exp_design") %>% collect()  %>%
        left_join(df, by = c("accession" = "sampleID")) %>%
        select(accession, plat.y, vial, rep, Genome, Final_Guess,
               Final_Best_Hit_Read_Numbers) %>%
        rename(plat = plat.y)

    df$Contam <- !(grepl("Salmonella", df$Genome))
    df_salmonella <- df %>% filter(Contam == FALSE) %>%
        group_by(accession, plat, vial) %>% summarize(total_prop = sum(Final_Guess),
                                                      prop_read = 1000000 *(1-total_prop))

    non_salmonella <- df %>% filter(Contam == TRUE)
}



### Figures -----------------------------------------------------
contam_counts_figure <- function(db_con){
    df_salmonella <- .genomic_purity_df(db_con)
ggplot(df_salmonella) + geom_boxplot(aes(x = plat,
                                         y = prop_read,
                                         color = plat)) +
    theme_bw() +
    labs(x = "Sequencing Platform",y = "Contaminants/Million Reads") +
    theme(legend.position = "none")
}
ggsave(filename = "contam_prop.png",width = 4, height = 3.5, dpi = 450)
### Distribution of contaminant by tax id
ggplot(non_salmonella) +
    geom_density(aes(x = Final_Best_Hit_Read_Numbers, fill = plat),
                 alpha = 0.5) +
    scale_x_log10() + theme_bw() +
    labs(x = "Number of Reads",y = "Density", fill= "Platform") +
    theme(legend.position = c(0.65,0.9), legend.direction = "horizontal")
ggsave(filename = "contam_count.png",width = 4, height = 3.5, dpi = 450)

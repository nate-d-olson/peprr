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

##### Consensus base -----------------------------------------------------------
## purity plots - script to generate purity plots


## Analysis of low purity positions
# figures correlating different parameters with low purity positions
# e.g. bias and coverage

# table of low purity positions

## indel analysis
# correlate with homopolymer annotation results

### Homogeneity ----------------------------------------------------------------
# table and figures in present_homogeneity.R



### Genomic Purity -------------------------------------------------------------
.genomic_purity_df <- function(db_con,genus){
    purity_df <- tbl(src = db_con, from = "purity") %>%
        dplyr::collect()

    df <- dplyr::tbl(src = peprDB, from = "exp_design") %>%
            dplyr::collect() %>%
            dplyr::left_join(purity_df) %>%
            dplyr::select(accession, plat, vial, rep, Genome, Final.Guess,
               Final.Best.Hit.Read.Numbers)

    df$Contam <- !(grepl(genus, df$Genome))
    df
}

#' create contamination count figure
#' @param db_con peprDB connection
#' @param genus rm genus
#' @return NULL
contam_counts_figure <- function(db_con, genus){
    df <- .genomic_purity_df(db_con, genus) %>%
                dplyr::filter(Contam == FALSE) %>%
                dplyr::group_by(accession, plat, vial) %>%
                dplyr::summarize(total_prop = sum(Final.Guess),
                                prop_read = 1000000 *(1-total_prop))
    ggplot2::ggplot(df) + ggplot2::geom_boxplot(ggplot2::aes(x = plat,
                                           y = prop_read,
                                           color = plat)) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Sequencing Platform",y = "Contaminants/Million Reads") +
        ggplot2::theme(legend.position = "none")
}

#' create contamination distribution figure
#' @param db_con peprDB connection
#' @param genus rm genus
#' @return NULL
contam_distribution_figure <- function(db_con,genus){
    df <-.genomic_purity_df(db_con, genus) %>% dplyr::filter(Contam == TRUE)
    ggplot2::ggplot(df) +
        ggplot2::geom_density(ggplot2::aes(x = Final.Best.Hit.Read.Numbers, fill = plat),
                     alpha = 0.5) +
        ggplot2::scale_x_log10() + ggplot2::theme_bw() +
        ggplot2::labs(x = "Number of Reads",y = "Density", fill= "Platform") +
        ggplot2::theme(legend.position = c(0.65,0.9), legend.direction = "horizontal")
}

#' create contaminant heatmap
#' @param db_con peprDB connection
#' @param genus rm genus
#' @return NULL
contam_heatmap_figure <- function(db_con,genus){
    df <-.genomic_purity_df(db_con, genus) %>% dplyr::filter(Contam == TRUE) %>%
        dplyr::filter(Final.Best.Hit.Read.Numbers > 1)
    genome_count_df <- df %>% dplyr::group_by(Genome)  %>%
        dplyr::summarize(count = n())  %>%
        dplyr::filter(count > 2)
    filt_df <- df %>% dplyr::filter(Genome %in% genome_count_df$Genome) %>%
        tidyr::separate(Genome, c("X1","ti", "X2", "org"),
                        sep = "\\|",remove = FALSE) %>%
        tidyr::separate(org, into = c("Genus", "species"),
                        sep = "_",extra = "drop") %>%
        dplyr::mutate(org_name = paste(Genus, species, sep = " ")) %>%
        dplyr::group_by(accession, org_name, plat, vial, rep) %>%
        dplyr::summarise(count = sum(Final.Best.Hit.Read.Numbers)) %>%
        dplyr::filter(org_name != "Unknown. NA") %>%
        dplyr::arrange(count)
    ggplot2::ggplot(filt_df) +
        ggplot2::geom_raster(ggplot2::aes(y = org_name, x = accession, fill = count)) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(~plat, scale = "free_x") +
        ggplot2::labs(x = "Datasets",y = "Organism", fill= "Hit Count") +
        ggplot2::theme(legend.position = "bottom",
                       legend.direction = "horizontal",
                       axis.text.x = ggplot2::element_text(angle = 90))
}
# need to work out how to best cleanup taxa information .....
# using taxize package
# tax_id <- unique(filt_df$org)  %>% str_replace_all("-", " ")  %>% map(.f = get_uid)  %>% map(.f = classification, db = "ncbi")
# for( i in tax_id){id <-names(i);id_df <- i[[id]]; id_df$id <- id;}

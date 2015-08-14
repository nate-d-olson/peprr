### Genomic Purity -------------------------------------------------------------
.genomic_purity_df <- function(db_con,genus, db_ref_root = "micro_rm_patho_db_"){
    purity_df <- dplyr::tbl(src = db_con, from = "purity") %>%
        dplyr::collect() %>%
        dplyr::mutate(accession = stringr::str_replace(pattern = db_ref_root,
                                                       replacement = "",accession))

    purity_df <- dplyr::tbl(src = db_con, from = "exp_design") %>%
        dplyr::collect() %>%
        dplyr::left_join(purity_df) %>%
        dplyr::select(accession, plat, vial, rep, Genome, Final.Guess,
                      Final.Best.Hit.Read.Numbers)

    purity_df$Contam <- !(grepl(genus, purity_df$Genome))
    purity_df
}

#' create contamination count figure
#' @param db_con peprDB connection
#' @param genus rm genus
#' @return NULL
contam_counts_figure <- function(db_con, genus){
    df <- .genomic_purity_df(db_con, genus) %>%
        dplyr::filter(Contam == TRUE) %>%
        dplyr::group_by(accession, plat, vial) %>%
        dplyr::summarize(contam_prop = sum(Final.Guess))
    ggplot2::ggplot(df) + ggplot2::geom_boxplot(ggplot2::aes(x = plat,
                                                             y = contam_prop,
                                                             color = plat)) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Sequencing Platform",y = "Proportion Total Contaminants") +
        ggplot2::theme(legend.position = "none")
}

#' create contamination distribution figure
#' @param db_con peprDB connection
#' @param genus rm genus
#' @return NULL
contam_distribution_figure <- function(db_con,genus){
    .genomic_purity_df(db_con, genus) %>%
        dplyr::filter(Contam == TRUE) %>%
        ggplot2::ggplot() +
        ggplot2::geom_density(ggplot2::aes(x = Final.Best.Hit.Read.Numbers,
                                           fill = plat),
                              alpha = 0.5) +
        ggplot2::scale_x_log10() + ggplot2::theme_bw() +
        ggplot2::labs(x = "Number of Reads",y = "Density",
                      fill= "Platform") +
        ggplot2::theme(legend.position = c(0.65,0.9),
                       legend.direction = "horizontal")
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

contam_point_line_figure <- function(db_con, genus){
    purity_df <- .genomic_purity_df(db_con, genus) %>%
        dplyr::filter(Contam == TRUE) %>%
        dplyr::filter(Final.Best.Hit.Read.Numbers > 1)

    genome_count_df <- purity_df %>%
        dplyr::group_by(Genome)  %>%
        dplyr::summarize(count = n())  %>%
        dplyr::filter(count > 2)

    filt_df <- purity_df %>% dplyr::filter(Genome %in% genome_count_df$Genome) %>%
        tidyr::separate(Genome, c("X1","ti", "X2", "org"),
                        sep = "\\|",remove = FALSE) %>%
        tidyr::separate(org, into = c("Genus", "species"),
                        sep = "_",extra = "drop") %>%
        dplyr::mutate(org_name = paste(Genus, species, sep = " ")) %>%
        dplyr::group_by(accession, org_name, plat, vial, rep) %>%
        dplyr::summarise(count = sum(Final.Best.Hit.Read.Numbers)) %>%
        dplyr::filter(org_name != "uncultured bacterium") %>%
        dplyr::arrange(count)

    org_order <- filt_df  %>% dplyr::group_by(org_name) %>%
        dplyr::summarise(med_count = median(count))  %>%
        dplyr::arrange(med_count)  %>% .$org_name

    filt_df$org_name <- factor(filt_df$org_name, levels = org_order)
    ggplot2::ggplot(filt_df) +
        ggplot2::geom_point(ggplot2::aes(y = count, x = org_name, color = plat)) +
        ggplot2::geom_line(ggplot2::aes(y = count, x = org_name, color = plat), alpha = 0.5) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Organism", y = "Reads", color = "Platform") +
        ggplot2::theme(legend.position = c(0.85, 0.15))
}

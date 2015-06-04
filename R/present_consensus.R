## code for generating base level purity results
# db_con <- dplyr::src_sqlite(db_path)
# library(dplyr)
# library(ggplot)
# cb tables  - raw vcf input
# pur tables - per sample, position purity summary - includes forward and reverse ref and alt counts
# pur_miseq_pooled - position purity summary by platform
# pur_join - pairwise table for purity scatter plot

## need to rerun make peprDB - wtih new data
## excluding novel seqeunces added by pilon

purity_scatter_plot <- function (db_con, plat1_name = "Miseq", plat2_name = "PGM") {
    dplyr::tbl(db_con, "pur_join") %>%
        dplyr::filter(plat1 < 0.98 | plat2 < 0.98, CHROM != "pilon_novel_001" ) %>%
        dplyr::collect() %>%
        dplyr::filter(CHROM %in% grep("novel", CHROM,value = TRUE,invert = TRUE)) %>%
        ggplot2::ggplot() +
            ggplot2::geom_point(ggplot2::aes(x = plat1, y = plat2),
                                alpha = 0.5) +
            ggplot2::labs(x = "MiSeq", y = "PGM") + ggplot2::theme_bw()
}

low_purity_table <- function(db_con){
    low_pur <- dplyr::tbl(db_con, "pur_join") %>%
        dplyr::filter(plat1 < 0.98, plat2 < 0.98, CHROM != "pilon_novel_001" ) %>% dplyr::collect() %>%
        dplyr::select(-plat1, -plat2)
    low_miseq <- dplyr::tbl(db_con, "pur_miseq_pooled") %>%
                    dplyr::filter(CHROM %in% low_pur$CHROM,
                                  POS %in% low_pur$POS) %>%
                    dplyr::collect()  %>%
                    dplyr::right_join(low_pur)
    low_pgm <- dplyr::tbl(db_con, "pur_pgm_pooled") %>%
        dplyr::filter(CHROM %in% low_pur$CHROM,
                      POS %in% low_pur$POS) %>%
        dplyr::collect()  %>%
        dplyr::right_join(low_pur)
    low_join <- dplyr::left_join(low_miseq, low_pgm, by = c("CHROM", "POS"))
    names(low_join) <- names(low_join) %>%
        stringr::str_replace_all("x","miseq") %>%
        stringr::str_replace_all("y","pgm")
    low_join
}

get_low_pur_metrics <- function(df, low_pur){
    df  %>% dplyr::filter(CHROM %in% low_pur$CHROM,
                  POS %in% low_pur$POS) %>%
        dplyr::collect()  %>%
        dplyr::right_join(low_pur) %>%
        dplyr::select(CHROM, POS, RPB, MQB, BQB, MQSB, MQ0F) %>%
        dplyr::group_by(CHROM, POS) %>%
        dplyr::summarise_each(funs(mean))
}
low_pur_metrics <- function(db_con){
    low_pur <- dplyr::tbl(db_con, "pur_join") %>%
        dplyr::filter(plat1 < 0.98, plat2 < 0.98) %>% dplyr::collect() %>%
        dplyr::select(-plat1, -plat2)
    low_miseq <- dplyr::tbl(db_con, "cb_miseq") %>%
        get_low_pur_metrics(low_pur) %>%
        dplyr::mutate(plat = "miseq")
    low_pgm <- dplyr::tbl(db_con, "cb_pgm") %>% get_low_pur_metrics(low_pur) %>%
        dplyr::mutate(plat = "pgm")
    dplyr::bind_rows(low_miseq, low_pgm)
}

low_pur_metrics_figure <- function(db_con){
    low_pur_metrics(db_con) %>% tidyr::gather("metric","value", 3:7) %>%
        ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = metric, y = value, color = plat)) +
            ggplot2::facet_wrap(~POS) + ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
}


gel_indel_metrics <- function(df){
     df %>% dplyr::filter(INDEL != 0) %>%
        select(CHROM, POS, REF,ALT, QUAL, IDV, IMF, INDEL) %>%
        group_by(CHROM,POS,REF, ALT) %>%
        summarise_each(funs(mean)) %>%
        collect()
}

indel_metrics <- function(db_con){
    miseq_indel <- dplyr::tbl(db_con, "cb_miseq") %>%
        gel_indel_metrics() %>%
        dplyr::mutate(plat = "miseq")
    dplyr::tbl(db_con, "cb_pgm") %>%
        gel_indel_metrics() %>%
        dplyr::mutate(plat = "pgm") %>%
        bind_rows(miseq_indel)
}

# ## indel scatter plot
# indel_df <- indel_metrics(db_con)
# indel_pair <- indel_df  %>% group_by(CHROM, POS)  %>%
#     select(-QUAL, -IDV)  %>%
#     tidyr::spread(key = plat, value = IMF,drop = TRUE) %>%
#     filter(!is.na(miseq), !is.na(pgm))
# ggplot(indel_pair) + geom_point(aes(x = miseq, y = pgm)) + ylim(0,1) + xlim(0,1)
# ggplot(indel_pair) + geom_point(aes(x = miseq, y = pgm)) + ylim(0.2,1) + xlim(0.02,1)

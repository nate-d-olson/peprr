### Homogeneity Data Frames ----------------------------------------------------
.convert_percents <- function(x){
    stringr::str_replace(x, pattern = "%", replacement = "") %>%
        as.numeric()
}

#' data frame for homogeneity snp results
#' @param db_con peprDB connection
#' @return data_frame
create_homogeneity_df <- function(db_con){
    dplyr::tbl(src = db_con, from = "varscan_snp")  %>%
        dplyr::select(chrom, position, ref, var, normal, tumor,
               normal_var_freq, tumor_var_freq,
               variant_p_value, somatic_p_value)  %>%
        dplyr::collect()  %>%
        dplyr::mutate(norm_freq = .convert_percents(normal_var_freq),
               tumor_freq = .convert_percents(tumor_var_freq)) %>%
        dplyr::select(-normal_var_freq, -tumor_var_freq) %>%
        tidyr::gather("type","accession", normal, tumor) %>%
        tidyr::gather("freq_type", "freq", norm_freq, tumor_freq)
}

#' data frame for pairwise homogeneity snp results
#' @param db_con peprDB connection
#' @return data_frame
create_pairwise_homogeneity_df <- function(db_con){
    dplyr::tbl(src = db_con, from = "varscan_snp")  %>%
        dplyr::select(chrom, position, ref, var, normal, tumor,
               normal_var_freq, tumor_var_freq,
               variant_p_value, somatic_p_value)  %>%
        dplyr::collect()  %>%
        dplyr::mutate(norm_freq = .convert_percents(normal_var_freq),
               tumor_freq = .convert_percents(tumor_var_freq)) %>%
        dplyr::rename(ds1=normal, ds2=tumor)
}

### Tables ---------------------------------------------------------------------
## Questions 1. what about correction for multiple comparisons, should it be
## 0.05/total pairs or number of pair where the variant was called?

#' data frame for pairwise homogeneity snp results
#' @param db_con peprDB connection
#' @param rename_col boolean if TRUE rename columns for use in ROA table
#' @return data_frame
homogeneity_sig_table <- function(db_con, rename_cols = FALSE){
    df <- create_pairwise_homogeneity_df(db_con)  %>%
        dplyr::mutate(sig_p = ifelse(somatic_p_value < 0.05, 1, 0)) %>%
        dplyr::group_by(position)  %>%
        # note the number of datasets used to calculate number of pairwise
        # comparisons is hard coded
        dplyr::summarize(prop_pairs = n()/sum(c(15:1)),
                         med_freq = median(norm_freq),
                         min_p = min(somatic_p_value),
                         sig_n = sum(sig_p))
    if(rename_cols){
        df <- df %>% dplyr::rename("Position" = position,
                                   "Proportion of Pairs" = prop_pairs,
                                   "Median Frequency" = med_freq,
                                   "Minimum P-value" = min_p,
                                   "N Significant" = sig_n)
    }
    df
}

### Figures --------------------------------------------------------------------

## Distribution of somatic variant call p-values filtering on positions with
## variants called for 1/4 of pairs and minimum p value less than 0.8
#' Figure with distribution of p values for pair-wise variants
#' @param db_con peprDB connection
#' @param prop_pairs_cutoff cutoff values for positions with proportions of pairs with variant call default 0.25
#' @param min_p_value cutoff values for positions with minimum variant call p-value default 0.8
#' @return figure
homogeneity_pvalue_figure <- function(db_con, prop_pairs_cutoff = 0.25, min_p_value = 0.8){
    sig_var <- homogeneity_sig_table(db_con)
    create_homogeneity_df(db_con) %>%
        dplyr::filter(position %in%
            sig_var$position[sig_var$prop_pairs > prop_pairs_cutoff &
                                    sig_var$min_p < min_p_value])  %>%
        ggplot2::ggplot() +
            ggplot2::geom_histogram(ggplot2::aes(x = somatic_p_value,
                                                    fill = as.factor(position))) +
            ggplot2::theme_bw() +
            ggplot2::facet_wrap(~position)+
            ggplot2::labs(x = "Distibution of P-values for Pairwise Comparisons",
                          y = "Count") +
            ggplot2::theme(legend.position = "none")
}

#' Figure with distribution of variant frequency for pair-wise variants
#' @param db_con peprDB connection
#' @param prop_pairs_cutoff cutoff values for positions with proportions of pairs with variant call default 0.25
#' @param freq_cutoff cutoff value for positions variant frequency default 0.98
#' @return figure
homogeneity_freq_figure <- function(db_con, prop_pairs_cutoff = 0.25,
                                    freq_cutoff = 99){
    sig_var <- homogeneity_sig_table(db_con)
    create_homogeneity_df(db_con) %>% dplyr::filter(position %in%
           sig_var$position[sig_var$prop_pairs > prop_pairs_cutoff &
                                sig_var$med_freq < freq_cutoff]) %>%
        ggplot2::ggplot() +
            ggplot2::geom_histogram(ggplot2::aes(x = freq,
                                                 fill = as.factor(position))) +
            ggplot2::theme_bw() +
            ggplot2::xlim(0,100) +
            ggplot2::labs(x = "Variant Frequency", y = "Count") +
            ggplot2::facet_wrap(~position) +
            ggplot2::theme(legend.position = "none")
}

## Pairwise plot of somatic p-values for positions with p value < 0.05
## Note position hard coded.
## not necessary at this point
# ggplot(pair_wise_homogeneity[pair_wise_homogeneity$position == 1103278,]) +
#     geom_raster(aes(x = ds1, y = ds2, fill = -10*log(somatic_p_value))) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90), legend.position = c(0.95,0.35)) +
#     labs(x = "Dataset", y = "Dataset", fill = "-10logP")


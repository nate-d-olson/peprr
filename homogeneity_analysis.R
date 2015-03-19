library(peprr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

peprDB_path = "~/Desktop/micro_rm/pepr-data/peprDB.sqlite"
## Load database
peprDB <- init_peprDB(db_path = peprDB_path,create = FALSE)

### Homogeneity Data Frames ----------------------------------------------------
.convert_percents <- function(x){
    str_replace(x, pattern = "%", replacement = "") %>%
        as.numeric()
}

homogeneity <- tbl(src = peprDB, from = "varscan_snp")  %>%
    select(chrom, position, ref, var, normal, tumor,
           normal_var_freq, tumor_var_freq,
           variant_p_value, somatic_p_value)  %>%
    collect()  %>%
    mutate(norm_freq = .convert_percents(normal_var_freq),
           tumor_freq = .convert_percents(tumor_var_freq)) %>%
    select(-normal_var_freq, -tumor_var_freq) %>%
    gather("type","accession", normal, tumor) %>%
    gather("freq_type", "freq", norm_freq, tumor_freq)

pair_wise_homogeneity <- tbl(src = peprDB, from = "varscan_snp")  %>%
    select(chrom, position, ref, var, normal, tumor,
           normal_var_freq, tumor_var_freq,
           variant_p_value, somatic_p_value)  %>%
    collect()  %>%
    mutate(norm_freq = .convert_percents(normal_var_freq),
           tumor_freq = .convert_percents(tumor_var_freq)) %>%
    rename(ds1=normal, ds2=tumor)

### Tables ---------------------------------------------------------------------
## Questions 1. what about correction for multiple comparisons, should it be
## 0.05/total pairs or number of pair where the variant was called?

sig_var<- pair_wise_homogeneity  %>% mutate(sig_p =
                                        ifelse(somatic_p_value < 0.05, 1, 0)) %>%
    group_by(position)  %>%
    summarize(prop_pairs = n()/sum(c(15:1)), med_freq = median(norm_freq),
              min_p = min(somatic_p_value), sig_n = sum(sig_p))
sig_table <- sig_var %>% rename("Position" = position,
                                "Proportion of Pairs" = prop_pairs,
                                "Median Frequency" = med_freq,
                                "Minimum P-value" = min_p,
                                "N Significant" = sig_n)

### Figures --------------------------------------------------------------------

## Distribution of somatic variant call p-values filtering on positions with
## variants called for 1/4 of pairs and minimum p value less than 0.8
filter(homogeneity, position %in%
                    sig_var$position[sig_var$prop_pairs > 0.25 &
                                         sig_var$min_p < 0.8])  %>%
    ggplot() + geom_histogram(aes(x = somatic_p_value, fill = as.factor(position))) +
    theme_bw() +
    facet_wrap(~position)+ labs(x = "Distibution of P-values for Pairwise Comparisons", y = "Count") + facet_wrap(~position) +
    theme(legend.position = "none")
ggsave(filename = "somatic_p_dist.png",width = 4, height = 3.5, dpi = 450)

## Distribution of variant frequency filtering on positions with
## variants called for 1/4 of pairs and freq < 0.05
filter(homogeneity, position %in%
           sig_var$position[sig_var$prop_pairs > 0.25 & sig_var$med_freq < 98]) %>%
    ggplot() + geom_histogram(aes(x = freq,
                                       fill = as.factor(position))) +
    theme_bw() + xlim(0,100) + labs(x = "Variant Frequency", y = "Count") + facet_wrap(~position) +
    theme(legend.position = "none")
ggsave(filename = "freq_dist.png",width = 4, height = 3.5, dpi = 450)

## Pairwise plot of somatic p-values for positions with p value < 0.05
## Note position hard coded.
ggplot(pair_wise_homogeneity[pair_wise_homogeneity$position == 1103278,]) +
    geom_raster(aes(x = ds1, y = ds2, fill = -10*log(somatic_p_value))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90), legend.position = c(0.95,0.35)) +
    labs(x = "Dataset", y = "Dataset", fill = "-10logP")


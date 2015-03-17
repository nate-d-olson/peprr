## Create db -------------------------------------------------------------------
library(peprr)

peprDB <- init_peprDB("~/Desktop/MG001",create = TRUE)

load_peprMeta("~/Desktop/micro_rm/pepr/RM8375_pipeline.yaml",db_con = peprDB)

load_metrics("~/Desktop/MG001/Re.Order.RC.CFSAN008157.HGAP_miseq_qc_stats/",
             db_con = peprDB)

load_fastqc("~/Desktop/MG001/Re.Order.RC.CFSAN008157.HGAP_miseq_qc_stats/",
            db_con = peprDB)

load_varscan("~/Desktop/MG001/Re.Order.RC.CFSAN008157.HGAP_miseq_homogeneity/",
             db_con = peprDB)
## error with real data replacement has 1 row data has 0

## =============================================================================
##
## Tables ----------------------------------------------------------------------
##
## =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(knitr)

### Dataset summary ------------------------------------------------------------
seq_metrics <- tbl(src = peprDB, from="align")  %>%
                select(accession, CATEGORY, TOTAL_READS, MEAN_READ_LENGTH)  %>%
                filter(CATEGORY %in% c("UNPAIRED","PAIR"))  %>%
                collect()  %>%
                separate(accession, c("ref","pilon", "accession"), sep = "_")

insert_tbl <- tbl(src = peprDB, from="insert_hist")  %>%
                group_by(accession)  %>%
                summarize(mean_insert =
                    sum(insert_size *All_Reads.fr_count)/sum(All_Reads.fr_count))  %>%
                collect()  %>% separate(accession, c("ref","pilon", "accession"), sep = "_")

seq_summary <- tbl(src = peprDB, from ="exp_design")  %>%
                    collect() %>% full_join(seq_metrics) %>%
                    full_join(insert_tbl) %>%
                    select(-ref)

# need to add coverage information
kable(seq_summary)

### Homogeneity summary --------------------------------------------------------
.convert_percents <- function(x){
    str_replace(x, pattern = "%", replacement = "") %>%
        as.numeric()
}

homogeneity <- tbl(src = peprDB, from = "varscan_indel")  %>%
                    select(chrom, position, ref, var, normal, tumor,
                           normal_var_freq, tumor_var_freq,
                           variant_p_value, somatic_p_value)  %>%
                    collect()  %>%
                    mutate(norm_freq = .convert_percents(normal_var_freq),
                           tumor_freq = .convert_percents(tumor_var_freq)) %>%
                    select(-normal_var_freq, -tumor_var_freq) %>%
                    gather("type","accession", normal, tumor) %>%
                    gather("freq_type", "freq", norm_freq, tumor_freq)

### Genomic Purity summary --------------------------------------------------------
pathoDF <- read.csv("~/Desktop/micro_rm/micro_rm_dev/analysis/stats/genomic_purity/pathoDF.csv",
                    stringsAsFactors = FALSE)
df <- join_all(list(pathoDF, sampleDF, meta))
#using pathoscope_analysis_v2.Rmd as guide

## =============================================================================
##
## Figures ---------------------------------------------------------------------
##
## =============================================================================

library(ggplot2)
ggplot(homogeneity) + geom_bar(aes(x = "accession", y = "freq"))

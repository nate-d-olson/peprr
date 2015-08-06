library(peprr)
library(dplyr)
source("rm_metadata.R")
peprDB <- dplyr::src_sqlite(db_path)



### Seq summary values
seq_summary <- seq_summary_table(peprDB)

pgm_med_read_length <- seq_summary %>% filter(Platform == "pgm") %>%
                            .[["Read Length"]] %>% median() %>% round(digits = 0)
miseq_med_read_length <- seq_summary %>% filter(Platform == "miseq") %>%
                            .[["Read Length"]] %>% median() %>% round(digits = 0)
pacbio_med_read_length <- seq_summary %>% filter(Platform == "pacbio") %>%
                            .[["Read Length"]] %>% median() %>% round(digits = 0)

mean_miseq_library_read_count <- seq_summary %>% filter(Platform == "miseq") %>%
                                    .[["Reads"]] %>% mean() %>% round(digits = 0)
mean_pgm_library_read_count <- seq_summary %>% filter(Platform == "pgm") %>%
                                    .[["Reads"]] %>% mean() %>% round(digits = 0)

mean_pgm_library_coverage <- seq_summary %>% filter(Platform == "pgm") %>%
    .[["Coverage"]] %>% mean() %>% round(digits = 0)
mean_miseq_library_coverage <- seq_summary %>% filter(Platform == "miseq") %>%
    .[["Coverage"]] %>% mean() %>% round(digits = 0)

pgm_total_coverage <- seq_summary %>% filter(Platform == "pgm") %>%
                        .[["Coverage"]] %>% sum() %>% round(digits = 0)
miseq_total_coverage <- seq_summary %>% filter(Platform == "miseq") %>%
                            .[["Coverage"]] %>% sum() %>% round(digits = 0)

pacbio_total_coverage <- seq_summary %>% filter(Platform == "pacbio") %>%
                            .[["Coverage"]] %>% sum() %>% round(digits = 0)

total_coverage <- miseq_total_coverage + pgm_total_coverage + pacbio_total_coverage


## base level purity
# in base_purity_analysis.R

## genomic contaminants
max_lib_contam <- peprr:::.genomic_purity_df(peprDB, rm_genus) %>%
                            filter(Contam == TRUE) %>%
                            group_by(accession) %>%
                            summarize(prop_contam = sum(Final.Guess)) %>%
                            .$prop_contam %>% max(na.rm = TRUE)
max_contam <- paste0(as.character(round(max_lib_contam,6) * 100),"%") # maximum contamination per dataset

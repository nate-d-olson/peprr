library(dplyr)
library(tidyr)
library(stringr)
library(knitr)

peprDB_path = "~/Desktop/micro_rm/pepr-data/peprDB.sqlite"
## Load database
peprDB <- init_peprDB(db_path = peprDB_path,create = FALSE)

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
    select(-ref, -pilon,-CATEGORY)

# need to add coverage information
kable(seq_summary,digits = 0)

library(dplyr)
source("gather_vcf_functions.R")
db_con <- src_sqlite("/Users/nolson/Desktop/micro_rm/pepr-data/peprDB.sqlite", create =TRUE)
i <- "pgm"
consensus_dir <- "/Users/nolson/Desktop/micro_rm/pepr-data/MG001/consensus_tsv"

load_consensus <- function(consensus_dir, db_con){
    for(i in c("miseq", "pgm")){
        tsv_file <- list.files(consensus_dir,
                               pattern = paste0("*",i,".tsv"),
                               full.names = TRUE)

        cb_table <- paste0("cb_",i)
        cbtsv_tosql(tsv_file, db_con,cb_table)

        pur_table <- paste0("pur_",i)
        pur_tbl(cb_table, db_con = db_con, pur_table)

        pur_plat_table <- paste0("pur_",i, "_pooled")
        pur_plat(pur_table, db_con = db_con, pur_plat_table)
    }
    pur_plat_join(pur_plat_tbl1 = "pur_miseq_pooled",
                  pur_plat_tbl2 = "pur_pgm_pooled",
                  db_con,
                  tbl_name = "pur_pooled_join",
                  plat1_name = "miseq", plat2_name = "pgm")

}

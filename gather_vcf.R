## Loads vcf files into sqlite database and creates new table with purity values
source("gather_vcf_functions.R")

# database
vcf_db <- src_sqlite("../../../data/RM8375/RM8375-v01.sqlite", create = TRUE)

#loading and calculating purity PGM
pgm_tsv <- "../../../analysis/stats/sequence_purity/RM8375-PGM.tsv"
cbtsv_tosql(pgm_tsv, vcf_db,"cb_pgm")
pur_tbl("cb_pgm", db_conn = vcf_db, "pur_pgm")
pur_plat("pur_pgm", vcf_db, "pur_pgm_pooled")

#loading and calculating purity MiSeq
miseq_tsv <- "../../../analysis/stats/sequence_purity/RM8375-MiSeq.tsv"
cbtsv_tosql(miseq_tsv, vcf_db,"cb_miseq")
pur_tbl("cb_miseq", db_conn = vcf_db, "pur_miseq")
pur_plat("pur_miseq", vcf_db, "pur_miseq_pooled")

# joining pooled platform data tables
pur_plat_join(pur_plat_tbl1 = "pur_miseq_pooled", pur_plat_tbl2 = "pur_pgm_pooled", 
              vcf_db, tbl_name = "pur_pooled_join", plat1_name = "MiSeq", plat2_name = "PGM")
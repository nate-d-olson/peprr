# loading database for use in developing code for tables and figures
library(peprr)
library(dplyr)
db_path <- "~/Desktop/micro_rm/pepr-data/MG001/peprDB.sqlite"
param_yaml <- "~/Desktop/micro_rm/pepr/RM8375_pipeline_reorder_rc_pilon.yaml"
qc_stats_dir <- "~/Desktop/micro_rm/pepr-data/MG001/Re.Order.RC.CFSAN008157.HGAP_miseq_qc_stats/"
homogeneity_dir <- "~/Desktop/micro_rm/pepr-data/MG001/Re.Order.RC.CFSAN008157.HGAP_homogeneity/"
consensus_dir <- "~/Desktop/micro_rm/pepr-data/MG001/Re.Order.RC.CFSAN008157.HGAP_consensus_base/"
pilon_dir <- "~/Desktop/micro_rm/pepr-data/MG001/Re.Order.RC.CFSAN008157.HGAP_pilon/"
purity_dir <- "~/Desktop/micro_rm/pepr-data/MG001/genomic_purity/"

createPeprDB(db_path, param_yaml, qc_stats_dir, homogeneity_dir,
             consensus_dir, purity_dir, pilon_dir)

## Functions for loading data into sqlit db
library(data.table)
library(dplyr)
library(tidyr)

cbtsv_tosql <- function (tsv_file,db_conn, tbl_name) {
  # tsv_file generated using vcf2tsv in vcflib
  # db_conn dplyr sqlite db connection
  # tbl_name of table creating in sqlite db
  vcf <- fread(tsv_file,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  vcf <- rename(vcf, PLATDP=DP) # added to fix issue with two DP columns
  copy_to(vcf_db, vcf, name = tbl_name, temporary = FALSE, indexes = list("CHROM","POS","SAMPLE")) 
}

# Calculating purity and filtering indels
pur_tbl <- function (vcf_tbl, db_conn, tbl_name) {
  # vcf_tbl sqlite table generated from a vcf file using cbtsv, note requires DP4 info for each sample
  # db_conn dplyr sqlite db connection
  # tbl_name of table creating in sqlite db
  
  # Select desired columns and filtering indels
  vcf_dp4 <- tbl(src = vcf_db, from = vcf_tbl)  %>% 
    select(CHROM, POS, SAMPLE, INDEL, DP, DP4, SP) %>% 
    filter(INDEL==0) %>% collect()
  
  vcf_dp4 <- vcf_dp4  %>% 
    separate(DP4, c("Ref_For","Ref_Rev","Alt_For","Alt_Rev")) %>%
    mutate(Ref_For = as.numeric(Ref_For),
           Ref_Rev = as.numeric(Ref_Rev),
           Alt_For = as.numeric(Alt_For),
           Alt_Rev = as.numeric(Alt_Rev),
           Ref=Ref_For + Ref_Rev, 
           Alt = Alt_For + Alt_Rev, 
           Pur = Ref/(Ref+ Alt))
  copy_to(db_conn, vcf_dp4, name=tbl_name,temporary = FALSE, indexes = list("CHROM","POS","SAMPLE"))
}

# Generating purity by platform summary
pur_plat <-function(pur_plat_tbl, db_conn, tbl_name){
  tbl(src = vcf_db, from=pur_plat_tbl) %>% 
    group_by(CHROM, POS) %>% 
    summarize(Ref = sum(Ref), Alt = sum(Alt)) %>%
    mutate(Pur = Ref/(Ref + Alt)) %>%
    compute(name = tbl_name, temporary = FALSE)
}

# Joining two purity platform tables
pur_plat_join <- function(pur_plat_tbl1, pur_plat_tbl2, vcf_db, tbl_name, plat1_name = "plat1", plat2_name = "plat2"){
  plat1_tbl <- tbl(src = vcf_db, from = pur_plat_tbl1) %>%
    group_by(CHROM,POS) %>%
    rename(plat1=Pur) %>%
    select(CHROM, POS, plat1)
  plat2_tbl <- tbl(src = vcf_db, from = pur_plat_tbl2) %>%
    group_by(CHROM,POS) %>% 
    rename(plat2=Pur) %>%
    select(CHROM, POS, plat2)
  # may want to change to outer_join later
  inner_join(plat1_tbl, plat2_tbl) %>% compute(name="pur_join", temporary = FALSE) 
}
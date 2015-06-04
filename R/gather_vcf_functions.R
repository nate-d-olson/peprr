# ## Functions for loading data into sqlite db
.cbtsv_tosql <- function (tsv_file,db_con,tbl_name) {
  # tsv_file generated using vcf2tsv in vcflib
  # db_conn dplyr sqlite db connection
  #tbl_name of table creating in sqlite db

  vcf <- data.table::fread(tsv_file,sep = "\t",header = TRUE,stringsAsFactors = FALSE) %>%
            dplyr::rename(PLATDP=DP) # added to fix issue with two DP columns
  dplyr::copy_to(db_con, vcf, name = tbl_name, temporary = FALSE,
                 indexes = list("CHROM","POS","SAMPLE"))
  ## removing from workspace
  rm(vcf)
}

# Calculating purity and filtering indels
.pur_tbl <- function (vcf_tbl, db_con,tbl_name) {
  # vcf_tbl sqlite table generated from a vcf file using cbtsv, note requires
  # DP4 info for each sample db_conn dplyr sqlite db connection tbl_name of
  # table creating in sqlite db select desired columns and filtering indels
  vcf_dp4 <-dplyr::tbl(src = db_con, from = vcf_tbl)  %>%
    dplyr::select(CHROM, POS, SAMPLE, INDEL, DP, DP4, SP) %>%
    filter(INDEL==0) %>% dplyr::collect()

  vcf_dp4 <- vcf_dp4  %>%
    tidyr::separate(DP4, c("Ref_For","Ref_Rev","Alt_For","Alt_Rev")) %>%
    dplyr::mutate(Ref_For = as.numeric(Ref_For),
           Ref_Rev = as.numeric(Ref_Rev),
           Alt_For = as.numeric(Alt_For),
           Alt_Rev = as.numeric(Alt_Rev),
           Ref=Ref_For + Ref_Rev,
           Alt = Alt_For + Alt_Rev,
           Pur = Ref/(Ref+ Alt))
  dplyr::copy_to(db_con, vcf_dp4, name=tbl_name,temporary = FALSE, indexes = list("CHROM","POS","SAMPLE"))
  rm(vcf_dp4)
}

# Generating purity by platform summary
.pur_plat <-function(pur_plat_tbl, db_con,tbl_name){
 dplyr::tbl(src = db_con, from=pur_plat_tbl) %>%
    dplyr::group_by(CHROM, POS) %>%
    dplyr::summarize(Ref = sum(Ref), Alt = sum(Alt)) %>%
    dplyr::mutate(Pur = Ref/(Ref + Alt)) %>%
    dplyr::compute(name =tbl_name, temporary = FALSE)
}

# Joining two purity platform tables
.pur_plat_join <- function(pur_plat_tbl1, pur_plat_tbl2, db_con,tbl_name,
                           plat1_name = "plat1", plat2_name = "plat2"){
  plat1_tbl <-dplyr::tbl(src = db_con, from = pur_plat_tbl1) %>%
    dplyr::group_by(CHROM,POS) %>%
    dplyr::rename(plat1=Pur) %>%
    dplyr::select(CHROM, POS, plat1)
  plat2_tbl <-dplyr::tbl(src = db_con, from = pur_plat_tbl2) %>%
    dplyr::group_by(CHROM,POS) %>%
    dplyr::rename(plat2=Pur) %>%
    dplyr::select(CHROM, POS, plat2)
  # may want to change to outer_join later
  dplyr::inner_join(plat1_tbl, plat2_tbl) %>%
      dplyr::compute(name="pur_join", temporary = FALSE)
  rm(plat1_tbl, plat2_tbl)
}

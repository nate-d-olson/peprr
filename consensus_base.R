#per base purity
calc_purity <- function(I16){
    return( sum(I16[1:2]) / sum(I16[1:4]) )
}

#per base purity probabilities
calc_pure_prob <- function(I16,p){
    return(pbinom(q=sum(I16[1:2]),size = sum(I16[1:4]),prob = p,))
}

#per base purity quantiles
calc_pure_prob <- function(I16,p){
    return(pbinom(q=sum(I16[1:2]),size = sum(I16[1:4]),prob = p,))
}

#calculate purity quantiles
calc_quantile <- function(I16,p){
    return(qbinom(p,size = sum(I16[1:4]), sum(I16[1:2]) / sum(I16[1:4]))/sum(I16[1:04]))
}

## calculating position probabilites and purity
process_vcf_purity <- function (vcf_file, ref, db_conn){
    #read vcf
    vcf <- VariantAnnotation::readVcf(vcf_file, geno=ref)

    dataset_name <- stringr::str_split(vcf_file,pattern = "_")  %>%
                        unlist() %>% dplyr::last() %>%
                        stringr::str_replace(".vcf", "")

    #   #get I16 info into a data.table
    I16_names <- c("R_Q13_F","R_Q13_R","NR_Q13_F","NR_Q13_R","RS_BQ",
                   "R_BQ_SSq","NR_BQ_S","NR_BQ_SSq","R_MQ_S","R_MQ_SSq",
                   "NR_MQ_S","NR_MQ_SSq","R_TD_S","R_TD_SSq","NR_TD_S","NR_TD_SSq")
    I16 <- plyr::ldply(VariantAnnotation::info(vcf)$I16) %>% dplyr::tbl_df() %>% data.table::setnames(I16_names)
    #
    #   # calculate purity
    PUR <- sapply(VariantAnnotation::info(vcf)$I16,FUN = calc_purity)
    PUR_prob97 <- sapply(VariantAnnotation::info(vcf)$I16,FUN=calc_pure_prob, p = 0.97)
    PUR_Q2.5 <- sapply(VariantAnnotation::info(vcf)$I16,FUN=calc_quantile, p = 0.025)
    PUR_Q50 <- sapply(VariantAnnotation::info(vcf)$I16,FUN=calc_quantile, p = 0.5)
    PUR_Q97.5 <- sapply(VariantAnnotation::info(vcf)$I16,FUN=calc_quantile, p = 0.975)
    #
    #generate datatable
    vcf_tbl <- data.table::data.table(CHROM = stringr::str_sub(string = rownames(VariantAnnotation::info(vcf)),
                                          start = 1,end = 8),
                          POS = IRanges::ranges(vcf)@start, WIDTH = IRanges::ranges(vcf)@width,
                          INDEL=VariantAnnotation::info(vcf)$INDEL, DP = VariantAnnotation::info(vcf)$DP,
                          QUAL = vcf@fixed$QUAL, PUR, PUR_prob97,
                          PUR_Q2.5, PUR_Q50, PUR_Q97.5)

    dplyr::copy_to(dest = db_con,df = vcf_tbl, name = stringr::str_c(dataset_name, "vcf", sep = "_"))

    I16$POS <- vcf_tbl$POS
    I16$WIDTH <- vcf_tbl$WIDTH
    I16$CHROM <- vcf_tbl$CHROM
    dplyr::copy_to(dest = db_con,df = I16,name = stringr::str_c(dataset_name, "I16", sep = "_"))

    rm(vcf,vcf_tbl,I16,PUR,PUR_prob97)
}

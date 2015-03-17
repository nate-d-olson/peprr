# ## Gather output from pathoscope and create a single data table for each loading into R.
# library(stringr)
# library(plyr)
# library(dplyr)
# source("pathoscopeR.R")
# genomic_purity_bioinf_dir = "../../bioinf/genome_purity/"
# pathoscope_params = "../../bioinf/genome_purity/RM8375_pathoscope_pipeline_params.txt"
#
# sampleDF <- parse_pathoparams(pathoscope_params)
# sampleDF$filename <- str_join(genomic_purity_bioinf_dir,sampleDF$sample,"/",
#                                 sampleDF$sample,"-sam-report.tsv",sep = "")
# write.csv(sampleDF, "sampleDF.csv")
#
# pathoDF <- ldply(list(sampleDF$filename)[[1]], parse_sam_report)
# write.csv(pathoDF, "pathoDF.csv")
#
#

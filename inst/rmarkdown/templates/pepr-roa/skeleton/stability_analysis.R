## gel image analysis
library(peprrDnaStability)
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
data_dir <- "~/Desktop/MG002-gels/images/cropped/"

batch_data <- batch_process_gels(data_dir)

## combine data frames
intensity_df <- dplyr::data_frame()
markers_df <- dplyr::data_frame()
binned_df <- dplyr::data_frame()
for(x in batch_data){
    gel_name <- x$gel
    intensity_df <- x$intensity_dat %>%
        dplyr::mutate(gel = gel_name) %>%
        dplyr::bind_rows(intensity_df)
    markers_df <- x$marker_dat %>%
        dplyr::mutate(gel = gel_name) %>%
        dplyr::bind_rows(markers_df)
    binned_df <- x$bin_dat %>%
        dplyr::mutate(gel = gel_name) %>%
        dplyr::bind_rows(binned_df)
}

## adding metadata
meta <- read.csv("~/Desktop/MG002-gels/metadata/MG002_stability_metadata.csv",
                 stringsAsFactors = F) %>% filter(condition != "Blank")
intensity <- left_join(intensity_df, meta) %>% filter(condition != "Blank")
markers <- left_join(markers_df, meta) %>% filter(condition != "Blank")
binned <- left_join(binned_df, meta) %>% filter(condition != "Blank")


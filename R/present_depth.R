## Depth Figures
make_bin_cov_df <- function (db_con, platforms = c("miseq","pgm", "pacbio")){
    cov_df <- dplyr::data_frame()
    for(plat in platforms){
        tbl_name <- paste0("depth_", plat)
        cov_df <- dplyr::tbl(src = db_con, from =tbl_name) %>%
            dplyr::group_by(SAMPLE, REF) %>%
            #dplyr::mutate(POS_BIN = round(POS/1000,0)) %>%
            dplyr::collect() %>%
            dplyr::mutate(POS_BIN = cut(POS,1000)) %>%
            dplyr::group_by(SAMPLE, REF,POS_BIN) %>%
            dplyr::summarise(POS = mean(POS),
                             MEAN_COV = mean(COV),
                             MIN_COV = min(COV),
                             MAX_COV = max(COV)) %>%
            dplyr::bind_rows(cov_df)
    }
    cov_df
}

## Coverage distribution
# ggplot(cov_df) +
#     geom_histogram(aes(x = MEAN_COV)) +
#     facet_grid(SAMPLE~REF, scales = "free")

## Coverage by genome position
# ggplot(cov_df) +
#     geom_smooth(aes(x = MEAN_POS, y = MEAN_COV, color = SAMPLE)) +
#     facet_wrap(~REF, scales = "free", ncol =1) + theme_bw()

## Coverage by GC
## Need to add ref to depth table and calculate GC

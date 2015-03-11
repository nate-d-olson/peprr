## supplemental
.parse_varscan_filename <- function(filename){
    dat_name <- stringr::str_split(filename,pattern = "_") %>% unlist()
    dat_name[length(dat_name)-1] %>%
        stringr::str_split(pattern = "-") %>% unlist()
}

.read_varscan <- function(file){
    df <- read.table(file, sep ="\t",
                     header = T, stringsAsFactors = F)
    ds_names <- .parse_varscan_filename(file)
    df$normal <- ds_names[1]; df$tumor <- ds_names[2]
    df
}



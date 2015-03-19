library(dplyr)
library(tidyr)
library(stringr)
library(knitr)



### Genomic Purity Data Frames -------------------------------------------------
pathoDF <- read.csv("~/Desktop/micro_rm/micro_rm_dev/analysis/stats/genomic_purity/pathoDF.csv",
                    stringsAsFactors = FALSE) %>%
                select(-X)
sampleDF <- read.csv("~/Desktop/micro_rm/micro_rm_dev/analysis/stats/genomic_purity/sampleDF.csv",
                    stringsAsFactors = FALSE) %>%
                select(-X)

df <- left_join(pathoDF, sampleDF)
df <- tbl(src = peprDB, from = "exp_design") %>% collect()  %>%
    left_join(df, by = c("accession" = "sampleID")) %>%
    select(accession, plat.y, vial, rep, Genome, Final_Guess,
           Final_Best_Hit_Read_Numbers) %>%
    rename(plat = plat.y)

df$Contam <- !(grepl("Salmonella", df$Genome))
df_salmonella <- df %>% filter(Contam == FALSE) %>%
    group_by(accession, plat, vial) %>% summarize(total_prop = sum(Final_Guess),
                                         prop_read = 1000000 *(1-total_prop))

non_salmonella <- df %>% filter(Contam == TRUE)
# %>% mutat
#     group_by(plat) %>% mutate(total_contam = sum(Final_Guess), prop_contam =)

### Figures -----------------------------------------------------
## Contaminant Read counts
ggplot(df_salmonella) + geom_boxplot(aes(x = plat,
                                         y = prop_read,
                                         color = plat)) +
    theme_bw() +
    labs(x = "Sequencing Platform",y = "Contaminants/Million Reads") +
    theme(legend.position = "none")
ggsave(filename = "contam_prop.png",width = 4, height = 3.5, dpi = 450)
### Distribution of contaminant by tax id
ggplot(non_salmonella) +
    geom_density(aes(x = Final_Best_Hit_Read_Numbers, fill = plat),
                 alpha = 0.5) +
    scale_x_log10() + theme_bw() +
    labs(x = "Number of Reads",y = "Density", fill= "Platform") +
    theme(legend.position = c(0.65,0.9), legend.direction = "horizontal")
ggsave(filename = "contam_count.png",width = 4, height = 3.5, dpi = 450)



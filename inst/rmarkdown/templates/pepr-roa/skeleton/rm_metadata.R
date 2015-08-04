# When using this template change for RM specific information
rm_number <- "RM8376"
rm_strain <- "Staphlyococcus aureus NGS100"
rm_genus <- "Staphylococcus"
rm_concentration <- 50
rm_volume <- 60
rm_mass <- round(rm_concentration*rm_volume/1000,0)
strain_source <- "Children's National Hostipal."
rm_genome_size <- 2.8*10^6

# Database location
## Requires a peprDB generated from the pepr pipeline
db_path <- "/Volumes/External-SSD/micro_rm/pepr-data/MG002/MG002.sqlite"

## Shipment info
box_count <- 19
box_size <- 81
exclusion_box <- 13
exclusion_size <- 38
vial_count <- (box_count-length(exclusion_box))*box_size+exclusion_size

## sequencing datasets
## number of pacbio replicate datasets
pb_replicate_number <- "three"

## optical mapping
## restriction enzyme used in optical mapping
restriction_enzyme <- "NcoI"

## Chrom rename
rename_chroms <- c("unitig_0|quiver|quiver|quiver" = "Chromosome",
                 "unitig_1|quiver|quiver|quiver" = "Plasmid")

## stability metadata
gel_dir <- "/Volumes/External-SSD/micro_rm/RM8376/stability_study/"
gel_metadata <- paste0(gel_dir,"MG002_stability_metadata.csv")
gel_data_dir <- paste0(gel_dir,"/images/cropped/")

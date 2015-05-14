# When using this template change for RM specific information
rm_number <- "RM8375"
rm_strain <- "Salmonella enterica subspecies enterica serovar Typhimurim strain LT2"
rm_genus <- "Salmonella"
rm_concentration <- 50
rm_volume <- 60
rm_mass <- round(rm_concentration*rm_volume/1000,0)
rm_vial_number <- 1500
strain_source <- "the Food and Drug Administration's Center for Food Safety, who originally received from The Institute for Genomic Research."
rm_genome_size <- 4.8*10^6

## DNA stability analysis
source("~/Desktop/dnaStability/stability_roa.Rdata")

# Database location
db_path <- "~/Desktop/micro_rm/pepr-data/MG001/peprDB.sqlite"

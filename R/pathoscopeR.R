#
# #################################################################################################################
# ##
# ## Summary: Function for processing pathoscope results, functions modified from
# ##            https://github.com/nate-d-olson/genomic_purity
# ## Date: 12/16/2014
# ## Author: Nate Olson
# ## Affiliation: National Institute for Standards and Technology
# ## Dependencies: reshape2, dplyr, stringr
# ## Output files:
# ##
# #################################################################################################################
# # require(reshape2)
# # require(stringr)
# # require(plyr)
# # require(dplyr)
#
# ## Example input file
# #====================
# #Total Number of Aligned Reads:  623  Total Number of Mapped Genomes:	9
# #Genome	MapQ Guess	MapQ Best Hit	MapQ Best Hit Read Numbers	MapQ High Confidence Hits	MapQ Low Confidence Hits	Alignment Guess	Alignment Best Hit	Alignment Best Hit Read Numbers	Alignment High Confidence Hits	Alignment Low Confidence Hits
# #ti|1133852|org|Escherichia_coli_O104:H4_str._2011C-3493	0.9999766754464884	0.9678972712680578	603.0	0.9678972712680578	0.0	0.9999766754464884	0.9678972712680578	603.0	0.9678972712680578	0.0
# #ti|585056|org|Escherichia_coli_UMN026	2.3324078726666834e-05	0.012841091492776886	8.0	0.012841091492776886	0.0	2.3324078726666834e-05	0.012841091492776886	8.0	0.012841091492776886	0.0
# #ti|585397|org|Escherichia_coli_ED1a	3.8955107904062245e-10	0.0032102728731942215	2.0	0.0032102728731942215	0.0	3.8955107904062245e-10	0.0032102728731942215	2.0	0.0032102728731942215	0.0
# #ti|562|org|Escherichia_coli	6.583931938866915e-11	0.009630818619582664	6.0	0.009630818619582664	0.0	6.583931938866915e-11	0.009630818619582664	6.0	0.009630818619582664	0.0
# #ti|216592|org|Escherichia_coli_042	1.939460620497064e-11	0.0016051364365971107	1.0	0.0016051364365971107	0.0	1.939460620497064e-11	0.0016051364365971107	1.0	0.0016051364365971107	0.0
# #ti|585055|org|Escherichia_coli_55989	2.5328608787508853e-24	0.0016051364365971107	1.0	0.0016051364365971107	0.0	2.5328608787508853e-24	0.0016051364365971107	1.0	0.0016051364365971107	0.0
# #ti|566546|org|Escherichia_coli_W	7.817298332120808e-26	0.0016051364365971107	1.0	0.0016051364365971107	0.0	6.305185886832915e-26	0.0008025682182985554	0.5	0.0016051364365971107	0.0Example input file format
#
# parse_sam_report <- function(inputfile,sampleID){#, sampleID){
#   # 'Parse pathoscope sam-report.tsv output files from pathoscope'
#   ## Example input file
#   #====================
#   #Total Number of Aligned Reads:  623  Total Number of Mapped Genomes:  9
#   #Genome	MapQ Guess	MapQ Best Hit	MapQ Best Hit Read Numbers	MapQ High Confidence Hits	MapQ Low Confidence Hits	Alignment Guess	Alignment Best Hit	Alignment Best Hit Read Numbers	Alignment High Confidence Hits	Alignment Low Confidence Hits
#   #ti|1133852|org|Escherichia_coli_O104:H4_str._2011C-3493	0.9999766754464884	0.9678972712680578	603.0	0.9678972712680578	0.0	0.9999766754464884	0.9678972712680578	603.0	0.9678972712680578	0.0
#   #ti|585056|org|Escherichia_coli_UMN026	2.3324078726666834e-05	0.012841091492776886	8.0	0.012841091492776886	0.0	2.3324078726666834e-05	0.012841091492776886	8.0	0.012841091492776886	0.0
#   #ti|585397|org|Escherichia_coli_ED1a	3.8955107904062245e-10	0.0032102728731942215	2.0	0.0032102728731942215	0.0	3.8955107904062245e-10	0.0032102728731942215	2.0	0.0032102728731942215	0.0
#   #ti|562|org|Escherichia_coli	6.583931938866915e-11	0.009630818619582664	6.0	0.009630818619582664	0.0	6.583931938866915e-11	0.009630818619582664	6.0	0.009630818619582664	0.0
#   #ti|216592|org|Escherichia_coli_042	1.939460620497064e-11	0.0016051364365971107	1.0	0.0016051364365971107	0.0	1.939460620497064e-11	0.0016051364365971107	1.0	0.0016051364365971107	0.0
#   #ti|585055|org|Escherichia_coli_55989	2.5328608787508853e-24	0.0016051364365971107	1.0	0.0016051364365971107	0.0	2.5328608787508853e-24	0.0016051364365971107	1.0	0.0016051364365971107	0.0
#   #ti|566546|org|Escherichia_coli_W	7.817298332120808e-26	0.0016051364365971107	1.0	0.0016051364365971107	0.0	6.305185886832915e-26	0.0008025682182985554	0.5	0.0016051364365971107	0.0Example input file format
#   if(file.info(inputfile)$size == 0){
#     return(data.frame())# what to change to include warning message
#   }
#
#   input_lines <- readLines(inputfile)
#
#   #Converting input to dataframe
#   col_names <- str_split(input_lines[2], pattern = "\t")
#   col_names <- str_replace_all(col_names[[1]]," ","_") #better way to pass list to str_replace
#   report <- colsplit(input_lines[c(-1,-2)], pattern = "\t", names = col_names)
#
#   #adding meta-data
#   meta <- str_split(input_lines[1], pattern = "\t")[[1]]
#   #report$sampleID <-sampleID
#   report$filename <- inputfile
#   report$Aligned_Reads <- meta[2]
#   report$Mapped_Genome <- meta[4]
#
#   return(report)
# }
#
# parse_pathoparams <- function(inputfile){
#   # function returns a data frame with the sampleID and platfor for the datasets analyzed
#   samples <-readLines(inputfile) %>% grep(pattern = "datasets",value = TRUE) %>% str_replace(pattern = "datasets=","")
#   samplesDF <- samples %>% str_split(pattern = ",") %>% ldply(colsplit,pattern = ":",names = c("sampleID","plat"))
#   return(samplesDF)
# }

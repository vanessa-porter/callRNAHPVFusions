#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla
.libPaths("/projects/vporter_prj/R/x86_64-centos7-linux-gnu-library/4.0")

## ---------------------------------------------------------------------------
## HPV Integration in Illumina Data
## Vanessa Porter, June 2022
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))

## ---------------------------------------------------------------------------
## LOAD INPUT 
## ---------------------------------------------------------------------------

# Make help options
option_list = list(
  make_option(c("-b", "--bed"), type="character", default=NULL,
              help="bed file with integration sites", metavar="character"),
  make_option(c("-o", "--out"), type="character", default = "il_integration_sites.bed",
              help="Output file name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$out

bed <- read.delim(opt$bed, header = F)
bed$name <- paste0("hpv.site", 1:nrow(bed))

write.table(bed, file = out, quote = F, sep = "\t", col.names = F, row.names = F)

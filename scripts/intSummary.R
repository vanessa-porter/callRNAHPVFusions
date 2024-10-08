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
  make_option(c("-v", "--vcf"), type="character", default=NULL,
              help="Filtered SV VCF file", metavar="character"),
  make_option(c("-d", "--depth"), type="character", default=NULL,
              help="Samtools depth file", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default=NULL,
              help="events bed file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = "integration_sites.txt",
              help="Output file name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir

## read in the files
#vcf <- read.delim("/projects/hpv_nanopore_prj/htmcp/illumina_call_integration/output/HTMCP-03-06-02006/int/hpvIntFilt.vcf", comment.char = "#", header = FALSE)
#depth <- read.delim("/projects/hpv_nanopore_prj/htmcp/illumina_call_integration/output/HTMCP-03-06-02006/int/hpvSitesDepth.txt", header = FALSE)
#bed <- read.delim("/projects/hpv_nanopore_prj/htmcp/illumina_call_integration/output/HTMCP-03-06-02006/int/il_integration_events.bed", header = FALSE)
vcf <- read.delim(opt$vcf, comment.char = "#", header = FALSE)
depth <- read.delim(opt$depth, header = FALSE)
bed <- read.delim(opt$bed, header = FALSE)

## ---------------------------------------------------------------------------
## Organize data 
## ---------------------------------------------------------------------------

# extract # of supporting reads from vcf info
info <- strsplit(vcf$V8, ";")
count <- lapply(info, function(x){grep("PAIR_COUNT", x, value = T)[2]})
count <- unlist(count)
count <- gsub("PAIR_COUNT=", "", count)
count <- as.numeric(count)

# clean up the HPV sites
hpv <- vcf$V5
hpv <- gsub("A|C|G|T", "", hpv)
hpv <- gsub('\\[|\\]', "", hpv)
hpvSites <- as.data.frame(matrix(data = unlist(strsplit(hpv, ":")), ncol = 2, byrow = T))

# get the hpv sites per event
bed$event <- paste0("event", 1:nrow(bed))
events <- strsplit(bed$V4, ",")
names(events) <- bed$event
sites <- unique(unlist(events))

# bring together in a dataframe
df <- data.frame(chr=vcf$V1, pos=vcf$V2, hpv.chr=hpvSites$V1, hpv.pos=hpvSites$V2, HPV.site=sites, event=NA,
                 nreads=count, depth=depth$V3)
for (i in 1:nrow(df)) {
  site <- df$HPV.site[i]
  df$event[i] <- names(events)[grep(paste0("\\b", site, "\\b"), events)]
}

df$VAF <- df$nreads / df$depth
eventbed <- bed[,c(1:3,5)]

## ---------------------------------------------------------------------------
## Save files
## ---------------------------------------------------------------------------

write.table(df, file = paste0(out, "/il_summary.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(eventbed, file = paste0(out, "/il_integration_events.bed"), quote = F, sep = "\t", row.names = F, col.names = F)

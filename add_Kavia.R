#!/usr/bin/Rscript

##################################################################################
############################### add_Kavia.R ###################################### 
##  PURPOSE: This script is used to add Kaviar annotations to VCF files:
##            - assumes the VCF file has already been functionally annotated
##            - adds Kaviar annotation to existing Info fields
##
##  INPUT: CGI or illumina VCF file
##  OUTPUT: VCF file contanining any previous annotations and Kaviar annotations
##        - tsv file of existing + kaviar annotations to add to impala database 
###################################################################################
############
## Setup  ##
############
source("user_config.R")
source("path_config.R")

##check if output directory exists, if not, create
dir.create(file.path(out_dir))

##create kaviar directory
dir.create(file.path(paste(out_dir, "kaviar", sep="/")))
kaviar_out <- paste(out_dir, "kaviar", sep="/")

##set working directory to output directory
setwd(out_dir)

##list required packages
##installation issues = install through biocLite()
packages <- c("VariantAnnotation", "IRanges")
##create function to install and/or load packages 
load_reqs <- function(lib_arg){
  if (require(lib_arg, character.only=TRUE)== TRUE){
    library(lib_arg, character.only=TRUE)
  } else {
    source("http://www.bioconductor.org/biocLite.R") 
    biocLite(lib_arg) 
  }
}
##apply function to list of required packages
lapply(packages, load_reqs)
  
##TODO 
## pipeline for complete genomics
## pipeline for illumina
##add param to read in only required fields
##add option for single vcf instead of directory 

#######################
## Prep Reference(s) ##
#######################
##RUNONCE only use when updating kaviar db
##convert commas in kaviar REF to periods to it doesn't break R
#system(paste("zcat <", kav_db, "| awk -F $'\t' 'BEGIN {OFS = FS} {gsub(\",\", \".\", $5);print}' | bgzip > Kav_DB_edit.vcf.gz"), wait=TRUE)
#paste("zcat <", test_vcf, "| awk -F $'\t' 'BEGIN {OFS = FS} {gsub(\",\", \".\", $5);print}' | bgzip > test_edit.vcf.gz")

##read in edited kaviar vcf and human ref
kaviar <- readVcf("/Documents/snp_annotato.R/ref_db/kav_edit.vcf.gz", humie_ref)

##rename chromosomes to match with vcf files
kaviar <- renameSeqlevels(kaviar, c("1"="chr1"))

##make data frame of needed kaviar fields
kav.df = cbind(chr= seqlevels(kaviar), as.data.frame(ranges(kaviar), stringsToFactors=FALSE), 
               ref=as.character(values(rowData(kaviar))[["REF"]]), alt=values(rowData(kaviar))[["ALT"]], info(kaviar))
##convert alts col to character list for matching
kav.df$alt.value = as.character(kav.df$alt.value)

#alt <CN0> what to do with CNV? 

##################################
## Gather VCF files to process  ## 
##################################
##data frame *.vcf.gz files in directory path
vcf_path <- data.frame(path=list.files(vcf_dir, pattern="*.vcf.gz$", full=TRUE))


##read in everything but sample data for speediness
vcf_param = ScanVcfParam(samples=NA)
vcf <- readVcf("./ref_db/vcf_edit.vcf.gz", humie_ref, param=vcf_param)

##data frame vcf info 
vcf.df = cbind(chr= seqlevels(vcf), as.data.frame(ranges(vcf), stringsToFactors=FALSE), fixed(vcf), info(vcf))

vcf.df = cbind(chr= seqlevels(vcf), as.data.frame(ranges(vcf), stringsToFactors=FALSE), 
               ref=as.character(values(rowData(vcf))[["REF"]]), alt=values(rowData(vcf))[["ALT"]], info(vcf))

##convert alt alleles to character for matching
vcf.df$alt.value = as.character(vcf.df$alt.value)


#################
## Match SNP's ##
#################
match1 = data.frame(merge(vcf.df, kav.df, by=c("chr", "start", "end", "ref")))

testing = head(data.frame(vcf_alt =match1$alt.value.x, kav_alt = match1$alt.value.y), 15)

testing[pmatch(testing$vcf_alt, testing$kav_alt),]
testing

?pmatch


##keep matches with same alt allele
matches = unique(match1[unlist(sapply(match1$ALT.x, grep, match1$ALT.y, fixed=TRUE)),])
nots = unique(match1[unlist(sapply(match1$ALT.x, grep, match1$ALT.y, fixed=TRUE, invert=TRUE)),])

data.frame(matches$ALT.x, matches$ALT.y)
data.frame(nots$ALT.x, nots$ALT.y)



##add column names
colnames(matches) = c("chr", "start", "end", "ref", "vcf_name", "vcf_alt", "kav_name", "kav_alt", "AF", "AC", "AN", "END", "DS")

?grep


##add the kaviar fields for all matches
vcf.df$AF = matches[match(vcf.df$name, matches$vcf_name),]$kav_AF
vcf.df$AC = matches[match(vcf.df$name, matches$vcf_name),]$kav_AC
vcf.df$AN = matches[match(vcf.df$name, matches$vcf_name),]$kav_AN
vcf.df$END = matches[match(vcf.df$name, matches$vcf_name),]$kav_END
vcf.df$DS = matches[match(vcf.df$name, matches$vcf_name),]$kav_DS

vcf_annots = cbind(vcf.df$AF, vcf.df$AC)
vcf_annots



Add new variables with the info<- setter.

> newInfo <- info(vcf)
> newInfo$foo <- 1:5
> info(vcf) <- newInfo
> names(info(vcf))
[1] "NS"  "DP"  "AF"  "AA"  "DB"  "H2"  "foo"
> info(vcf)$foo
[1] 1 2 3 4 5














######################
### unused snippets ##
######################

##add "chr" to chromosome name for matching with vcf files
system(paste("zcat <", kav_path, 
             "| awk 'BEGIN{OFS=\"\t\"} {if ($1 !~ /^#/) {printf \"chr\"} print $0}' | bgzip >", paste(kaviar_out, "kav_edit.vcf.gz", sep="/"), 
             sep=" "), wait=TRUE)

##index the kav file
system(paste("tabix", paste(kaviar_out, "kav_edit.vcf.gz", sep="/"), sep=" "), wait=TRUE)
##unzip vcf file and create copy to convert to txt
system(paste("gzcat < ", vcf, " > ", 
      paste(kaviar_out, lapply(strsplit(basename(vcf), "[.]"), 
                               function(x) paste(x[1], x[2], sep=".")), sep="/"), sep=""), wait=TRUE)
##read in uncompressed vcf file
vcf.df = read.csv(paste(kaviar_out, lapply(strsplit(basename(vcf), "[.]"), 
                                            function(x) paste(x[1], x[2], sep=".")), sep="/"), sep="\t")
dim(vcf.df)


##building of ranged kaviar file
idx = indexTabix("/Users/selasady/snp_annotato.R/ref_db/Kaviar.vcf.gz", format="vcf")
tab = TabixFile("/Users/selasady/snp_annotato.R/ref_db/Kaviar.vcf.gz", idx)
rng <- GRanges(seqnames="1", ranges=IRanges(start=c(39998057, 39998057), end=c(39999366, 39999366)))
kav_ranged = readVcf(tab, "hg19", param=rng)
ranges(kav_ranged)

ranges(kaviar)
alt(kaviar)[1] = "A,TT"
head(alt(kaviar))

writeVcf(kaviar, filename="/Documents/snp_annotato.R/ref_db/kaviar_comma.vcf.gz")







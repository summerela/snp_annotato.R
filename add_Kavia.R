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

##check if output directory exists, if not, create
dir.create(file.path(out_dir))

##create kaviar directory
dir.create(file.path(paste(out_dir, "kaviar", sep="/")))
kaviar_out <- paste(out_dir, "kaviar", sep="/")

##set working directory to output directory
setwd(out_dir)

##list required packages
##installation issues = install through biocLite()
packages <- c("VariantAnnotation", "parallel")
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
##add option for list of reference db
##add option for single vcf instead of directory 

#######################
## Prep Reference(s) ##
#######################
##add "chr" to chromosome name for matching with vcf files
system(paste("zcat <", kav_path, 
             "| awk 'BEGIN{OFS=\"\t\"} {if ($1 !~ /^#/) {printf \"chr\"} print $0}' | bgzip >", paste(kaviar_out, "kav_edit.vcf.gz", sep="/"), 
             sep=" "), wait=TRUE)

##index the kav file
system(paste("tabix", paste(kaviar_out, "kav_edit.vcf.gz", sep="/"), sep=" "), wait=TRUE)











##unzip vcf file and create copy to convert to txt
system(paste("gzcat < ", ref_db, " > ", 
             paste(kaviar_out, lapply(strsplit(basename(ref_db), "[.]"), 
                                      function(x) paste(x[1], x[2], sep=".")), sep="/"), sep=""), wait=TRUE)
##read in uncompressed vcf file
ref.df = read.table(paste(kaviar_out, lapply(strsplit(basename(ref_db), "[.]"), 
                                           function(x) paste(x[1], x[2], sep=".")), sep="/"), sep="\t")
colnames(ref.df) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
head(ref.df, 20)

kav = readVcf("./Kaviar_test.vcf.gz", genome="hg19")

info(header(kav))



##################################
## Gather VCF files to process  ## 
##################################
##data frame *.vcf.gz in directory path

vcf_dir.df <- data.frame(path=list.files(vcf_dir, pattern="*.vcf.gz$", full=TRUE))


vcf = "../../cg_test.vcf.gz"

##read in the vcf file header info
vcf_hdr = info(scanVcfHeader(vcf))
vcf_hdr






######################
### unused snippets ##
######################
##unzip vcf file and create copy to convert to txt
system(paste("gzcat < ", vcf, " > ", 
      paste(kaviar_out, lapply(strsplit(basename(vcf), "[.]"), 
                               function(x) paste(x[1], x[2], sep=".")), sep="/"), sep=""), wait=TRUE)
##read in uncompressed vcf file
vcf.df = read.csv(paste(kaviar_out, lapply(strsplit(basename(vcf), "[.]"), 
                                            function(x) paste(x[1], x[2], sep=".")), sep="/"), sep="\t")
dim(vcf.df)

##convert commas in kaviar to periods to it doesn't break R
zcat < Kaviar_test.vcf.gz | awk -F $'\t' 'BEGIN {OFS = FS} {gsub(",", "_", $5);print}' > test.vcf
bgzip test.vcf

## KAVIAR ##
##requires ref genome version/chr/coords
kaviar = readVcf("/Users/selasady/Documents/Variants/test.vcf.gz", "hg19")

##rename seqlevels so they can be compared
kaviar <- renameSeqlevels(kaviar, c("1"="chr1"))

##########################
## classes and accessor ##
##########################
##Chr, range, REF, ALT, QUAL, FILTER  
pos_info = rowData(vcf)
pos_info

vcf.df = data.frame(chr = head(seqnames(vcf)), start=head(start(vcf)), 
                  end= head(end(vcf)), ref=as.character(refs), alt=as.character(unlist(alts)))

kav.df = data.frame(chr = head(seqnames(kaviar)), start=head(start(kaviar)), 
                    end= head(end(kaviar)))
kav.df                    
                    , ref=as.character(k_refs), alt=as.character(unlist(k_alts)))

rowData(vcf)
rowData(kaviar)
unlist(k_alts)
k_alts

k_refs = head(ref(kaviar))
k_alts = head(alt(kaviar))

as.character(unlist(k_alts))
k_alts

names(k_alts) <- k_alts

refs = head(ref(vcf))
alts = head(alt(vcf))

anno.df = as(anno, "DataFrame")
elementMetadata(vcf) = cbind(vcf.df, anno.df)

##fxn annots
annot_info = info(vcf)

##header descriptions
hdr = info(header(vcf))

##genotype descrips
gt_info = geno(header(vcf))

##samples
samples = samples(header(vcf))







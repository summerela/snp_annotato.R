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
packages <- c("VariantAnnotation", "stringi")
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
## pipeline for illumina
##add option for single vcf instead of directory 

#########################
## Read in Kaviar File ##
#########################
##RUNONCE only use when updated kaviar db
##replace <CN0> with a period so we dont' break R
#zcat < Kav_ranged.vcf.gz| awk -F $'\t' 'BEGIN {OFS = FS} {gsub("\<CN0\>", ".", $5);print}' | bgzip > kav_edit.vcf.gz

##set reference 
humie_ref = "hg19"

##grab kaviar allele freq
kav_params = VRangesScanVcfParam(info="AF")

##read in kaviar vcf as VRange object
kaviar <- readVcfAsVRanges("/Documents/snp_annotato.R/ref_db/Kav_ranged.vcf.gz", humie_ref, param=kav_params)

##rename chromosomes to match with vcf files
kaviar <- renameSeqlevels(kaviar, c("1"="chr1"))

##data frame of kaviar info
kav.df = data.frame(chr =as.character(seqnames(kaviar)), start = start(kaviar), end = end(kaviar), 
                    ref = as.character(ref(kaviar)), alt=alt(kaviar), mcols(kaviar), stringsAsFactors=FALSE)
kav.df$id=sprintf("kav_%s", seq(1:length(rownames(kav.df)))) 

##################################
## Gather VCF files to process  ## 
##################################
##data frame *.vcf.gz files in directory path
#vcf_path <- data.frame(path=list.files(vcf_dir, pattern="*.vcf.gz$", full=TRUE))

##grab info from vcf header
hdr = scanVcfHeader("/Documents/snp_annotato.R/ref_db/cg_test.vcf.gz")

##read in everything but sample data for speediness
vcf_param = VRangesScanVcfParam(info=rownames(info(hdr)))
vcf <- readVcfAsVRanges("/Documents/snp_annotato.R/ref_db/cg_test.vcf.gz", humie_ref, param=vcf_param)

#################
## Match SNP's ##
#################
##create data frames of info to match on
vcf.df = data.frame(chr =as.character(seqnames(vcf)), start = start(vcf), end = end(vcf), ref = as.character(ref(vcf)), 
                    alt=alt(vcf), mcols(vcf), stringsAsFactors=FALSE)
##create row id for matching
vcf.df$id= sprintf("vcf_%s", seq(1:length(rownames(vcf.df))))

##merge based on all positional fields except vcf
col_match = data.frame(merge(vcf.df, kav.df, by=c("chr", "start", "end", "ref")))

##split each alt column by comma and bind together
M1 <- stri_list2matrix(sapply(col_match$alt.x,strsplit,','))
M2 <- stri_list2matrix(sapply(col_match$alt.y,strsplit,','))
M <- rbind(M1,M2)

##compare results
result <- apply(M,2,function(z) unique(na.omit(z[duplicated(z)])))

##add results column to col_match df for checking/subsetting
col_match$match = as.character(result)

####################
## Add annotation ##
####################
##subset only rows with matches
annots = col_match[which(col_match$match != "character(0)"),]

##for every vcf.df$id in annots$id.x, add the kaviar annotation
vcf.df$kav_AF = as.character(annots[match(vcf.df$id, annots$id.x),]$AF)

##merge information to add to existing, and make format pretty
new_annots = mcols(vcf)
new_annots$kav_AF = as.numeric(vcf.df$kav_AF)
mcols(vcf) = new_annots

##convert back to collapsed vcf
multi_vcf = asVCF(vcf, info=info(hdr), geno=geno(hdr), 

asVCF(x, info = character(), filter = character(), meta = character()
      
      
      vcf.df = data.frame(chr =as.character(seqnames(vcf)), start = start(vcf), end = end(vcf), ref = as.character(ref(vcf)), 
                          alt=alt(vcf), mcols(vcf), stringsAsFactors=FALSE)

##write to file
writeVcf(vcf, filename="./out_test.vcf")

?writeVcf

test


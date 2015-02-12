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
packages <- c("VariantAnnotation")
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
#paste("zcat <", kav_db, "| awk -F $'\t' 'BEGIN {OFS = FS} {gsub(\",\", \".\", $5);print}' | bgzip > Kav_DB.vcf.gz")



?indexTabix
idx = indexTabix("/Users/selasady/snp_annotato.R/ref_db/Kaviar.vcf.gz", format="vcf")
tab = TabixFile("/Users/selasady/snp_annotato.R/ref_db/Kaviar.vcf.gz", idx)

rng <- GRanges(seqnames="1", ranges=IRanges(start=c(39998057, 39998057), end=c(39999366, 39999366)))

kav_ranged = readVcf(tab, "hg19", param=rng)
ranges(kav_ranged)

writeVcf(kav_ranged, filename="/Users/selasady/snp_annotato.R/Kav_ranged.vcf")

?writeVcf


##read in edited kaviar vcf and human ref
kaviar <- readVcf("./kaviar/Kav_DB.vcf.gz", humie_ref)

##rename seqlevels so they can be compared
kaviar <- renameSeqlevels(kaviar, c("1"="chr1"))

##make data frame of kav info
kav.df = data.frame(chr = as.character(head(seqnames(kaviar))), start = head(start(ranges(kaviar))), end=as.integer(end(head(ranges(kaviar)))),
                    ref = as.character(head(ref(kaviar))), alt = as.character(CharacterList(head(alt(kaviar)))))

##add column of annotations
annot = head(info(kaviar))
annots = paste(as.character(annot$AF), as.character(annot$AC), as.character(annot$AN), as.character(annot$END), as.character(annot$DS), sep="|")
kav.df$annots = annots

##replace "NA" in annotation with blank
kav.df$annots = gsub("NA", "", as.character(kav.df$annots))

##################################
## Gather VCF files to process  ## 
##################################
##data frame *.vcf.gz files in directory path
vcf.df <- data.frame(path=list.files(vcf_dir, pattern="*.vcf.gz$", full=TRUE))

test_vcf = vcf.df$path[1]

ranges(vcf)

##read in vcf file(s)
vcf <- readVcf(as.character(test_vcf), humie_ref)



##data frame vcf info
head(samples(vcf))

?subsetByOverlaps
colnames(vcf)
rownames(vcf)

test =subsetByOverlaps(ranges(vcf), ranges(kaviar))

test

head(ranges(vcf))
head(ranges(kaviar))


######################



Two questions therefore:
  > 1) how do you add an INFO field in a way that doesn?t generate error
> from within R

> fl <- system.file("extdata", "ex2.vcf",
                    package="VariantAnnotation")
> vcf <- readVcf(fl, "hg19")

Original info variables.

> vcf <- readVcf(fl, "hg19")
> names(info(vcf))
[1] "NS" "DP" "AF" "AA" "DB" "H2"
> dim(info(vcf))
[1] 5 6

Add new variables with the info<- setter.

> newInfo <- info(vcf)
> newInfo$foo <- 1:5
> info(vcf) <- newInfo
> names(info(vcf))
[1] "NS"  "DP"  "AF"  "AA"  "DB"  "H2"  "foo"
> info(vcf)$foo
[1] 1 2 3 4 5


Query VCF class in R:
  
  To match on position I'd use subsetByOverlaps() on the VCF GRanges if
you want a GRanges back. If you only want a count of the hits (i.e.,
location if present) then use countOverlaps(). This is essentially
what
you're doing now but if chromosome and/or strand are important then
this
is the best way to go. Others free to chime in.

You could gain a little performance if the data are large and
chromosome
and strand are not important. Using countOverlaps() on the ranges()
extracted from the GRanges will save a some time on big data.

To match on variant name I would match on the colnames(vcf) as you are
doing.

Are any of these operations taking prohibitively long? Please let me
know if they are.
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










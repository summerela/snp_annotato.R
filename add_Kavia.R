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
system(paste("zcat <", kav_db, "| awk -F $'\t' 'BEGIN {OFS = FS} {gsub(\",\", \".\", $5);print}' | bgzip > Kav_DB_edit.vcf.gz"), wait=TRUE)

##read in edited kaviar vcf and human ref
kaviar <- readVcf("Kav_DB_edit.vcf.gz", humie_ref)

##rename chromosomes to match with vcf files
kaviar <- renameSeqlevels(kaviar, c("1"="chr1"))

##make data frame of kav info
kav.df = data.frame(name= as.character(head(names(rowData(kaviar)))),chr = as.character(head(seqnames(kaviar))), start = head(start(ranges(kaviar))), end=as.integer(end(head(ranges(kaviar)))),
                    ref = as.character(head(ref(kaviar))), alt = as.character(CharacterList(head(alt(kaviar)))), AF=unlist(head(info(kaviar)$AF)), AC=unlist(head(info(kaviar)$AC)),
                    AN=head(unlist(info(kaviar)$AN)), END = head(unlist(info(kaviar)$END)), DS = head(unlist(info(kaviar)$DS)))

# ##add column of annotations
# annot = head(info(kaviar))
# annots = paste(as.character(annot$AF), as.character(annot$AC), as.character(annot$AN), as.character(annot$END), as.character(annot$DS), sep="|")
# kav.df$annots = annots

##replace "NA" in annotation with blank
#kav.df$annots = gsub("NA", "", as.character(kav.df$annots))

##replace periods in Ref with "|"
#kav.df$alt = gsub("[.]", "|", as.character(kav.df$alt))

##################################
## Gather VCF files to process  ## 
##################################
##data frame *.vcf.gz files in directory path
vcf.df <- data.frame(path=list.files(vcf_dir, pattern="*.vcf.gz$", full=TRUE))

test_vcf = as.character(vcf.df$path[1])

##read in vcf file(s)
vcf <- readVcf(test_vcf, humie_ref)

##data frame vcf info
vcf.df = data.frame(name= as.character(head(names(rowData(vcf)))), chr = as.character(head(seqnames(vcf))), start = head(start(ranges(vcf))), end=as.integer(end(head(ranges(vcf)))),
                    ref = as.character(head(ref(vcf))), alt = as.character(CharacterList(head(alt(vcf)))), kav_annots="NA")

#################
## Match SNP's ##
#################
##generate list of identical matches, not including "alt"
match1 = merge(kav.df, vcf.df, by=c("chr", "start", "end", "ref"))
match1
##compare alts (separate to account for lists)
matches = unique(match1[unlist(sapply(match1$alt.x, grep, match1$alt.y)),])
colnames(matches) = c("chr", "start", "end", "ref", "kav_name", "kav_alt", "vcf_name", "vcf_alt", "kav_AF", "kav_AC", "kav_AN", "kav_END", "kav_DS")
matches

##for these guys
vcf.df

vcf.df[match(matches$vcf_name, vcf.df$name),]
  
kav.df[match(matches$kav_name, kav.df$name),]

##add this

colnames(kav.df)


##colnames info
names(info(vcf))
dim(info(vcf)) #170 26

Add new variables with the info<- setter.

newInfo <- info(vcf)
newInfo$foo <- 1:5
info(vcf) <- newInfo
names(info(vcf))
info(vcf)$foo




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
writeVcf(kav_ranged, filename="/Users/selasady/snp_annotato.R/Kav_ranged.vcf")







options(scipen = 999)
library(data.table)
library(ggplot2)
library(mgcv)
library(grid)
source("general_functions.R")
OUT_DIR = Sys.getenv('OUT_DIR')
source("./pheno_adjust_function.R")
save_loc = paste0(OUT_DIR, "/smooth_adjust")

if (1 == 1) {
RAW_TRAITS_ADD=Sys.getenv('RAW_TRAITS_ADD')
traits = read.csv(paste0(Sys.getenv('OUT_DIR'), "/pheno_names_new.tsv"), header=F, stringsAsFactor=F)$V1
traits = traits[!grepl("adj_", traits)]
pheno = fread(RAW_TRAITS_ADD, colClasses="str")
##traits = c(c("neut", "mono", "lymph", "baso", "eo", "gran"), traits)
print(traits)
traits = c(c("neut", "mono", "lymph", "baso", "eo", "gran"), traits)
for (trait in traits) {
    print(trait)
    adjusted_trait = adjust_variable_for_gwas(pheno, paste0(trait_gwas_to_pipeline(trait), "_tech_adj"), "INTERVAL", pipeline="additional", get_trait_list("sysmex0100_traits"), get_trait_list("sysmex01_traits"), get_trait_list("sysmex_positive_traits"))$adjusted_trait
    pheno[, eval(paste0(trait, "_age_adj"))] = adjusted_trait 
}
## get rid of duplicates
dup = pheno$genotyping_id[duplicated(pheno$genotyping_id)]
subset = !is.element(pheno$genotyping_id, dup) | (is.element(pheno$genotyping_id,dup) & pheno$interval == "baseline")
pheno = pheno[subset,]

pheno = lapply(pheno,  function(x) {ifelse(x == "", NA, x)})
pheno$clinic = pheno$centre
fwrite(pheno, paste0(OUT_DIR, "/pheno_internal_age_adj.csv"), sep="\t", na="NA", quote=FALSE)
}

if (1 == 1) {
RAW_TRAITS=Sys.getenv('RAW_TRAITS')
traits = read.csv(paste0(Sys.getenv('OUT_DIR'), "/pheno_names_orig.tsv"), header=F, stringsAsFactor=F)$V1
traits = traits[!grepl("adj_", traits)]
pheno = fread(RAW_TRAITS, colClasses="str")

print(traits)
for (trait in traits) {
    print(trait)
    adjusted_trait = adjust_variable_for_gwas(pheno, paste0(trait, "_tech_adj"), "INTERVAL", pipeline="", trait_pipeline_to_gwas(get_trait_list("sysmex0100_traits")), trait_pipeline_to_gwas(get_trait_list("sysmex01_traits")), trait_pipeline_to_gwas(get_trait_list("sysmex_positive_traits")))$adjusted_trait
    pheno[, eval(paste0(trait, "_age_adj"))] = adjusted_trait 
}
## get rid of duplicates
print(nrow(pheno))
dup = duplicated(pheno$genotyping_id)
if (sum(dup) == 0) {
    stop("there are duplicates")
}
print(nrow(pheno))
#pheno = lapply(pheno,  function(x) {ifelse(x == "", NA, x)})
pheno$clinic = pheno$centre
pheno[, (names(pheno)) := lapply(.SD, function(x) {ifelse(x == "", NA, x)}), .SDcols=names(pheno)]
pheno = cbind(data.table("FID"=pheno$genotyping_id, "IID"=pheno$genotyping_id), pheno)
fwrite(pheno, paste0(OUT_DIR, "/pheno_orig_age_adj.csv"), sep=" ", na="NA", quote=FALSE)
}

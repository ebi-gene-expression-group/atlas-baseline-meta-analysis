#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-i", "--input"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to rdata."
  ),
  make_option(
    c("-c", "--covariate"),
    action = "store",
    default = "~organism_part",
    type = 'character',
    help = "Covariate for batch correction, Default: organism_part"
  ),
  make_option(
    c("--output"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to output"
  ),
  make_option(
    c("--tsv_corrected_counts"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to tsv counts output'
  )
)

#batch effect correction packages
suppressPackageStartupMessages(require(sva)) #for ComBat
suppressPackageStartupMessages(require(SummarizedExperiment))
#library(RUVSeq) #for RUVs
#library(batchelor) #for mnnCorrect
suppressPackageStartupMessages(require(magrittr))
#library(stringr)
suppressPackageStartupMessages(require(purrr))

opt <- parse_args(OptionParser(option_list=option_list))

correct_batch_effect<-function(experiment, covariate, method=c('ComBat','RUV','MNN'), k){
  log<-experiment@assays@data %>% names %>% switch(log_counts=TRUE, counts=FALSE)
  model.data<-model.frame(covariate, experiment@colData[all.vars(covariate)])
  assays<-list()
  if(method == "ComBat") {
    print("Running ComBat...")
    #assays$corrected_counts <- ComBat(experiment@assays$data[[1]], experiment$batch, mod=model.matrix(model, data=model.data))
    assays$corrected_counts <- ComBat_seq(experiment@assays@data[[1]], experiment$batch, covar_mod=model.matrix(covariate, data=model.data))
  } else if(method == "RUV") {
    print("Running RUV...")
    assays$corrected_counts <- RUVs(experiment@assays$data[[1]], cIdx=seq_len(nrow(experiment@assays$data[[1]])), k=k,
                                    scIdx=model.data %>% expand.grid %>% apply(1,paste) %>% makeGroups, isLog=log)$normalizedCounts
  } else if(method == "MNN") {
    print("Running MNN...")
    assays$corrected_counts <- mnnCorrect(experiment@assays$data[[1]], batch=experiment$batch, k=k)@assays$data$corrected
  }
  return(SummarizedExperiment(
    assays = assays,
    colData = experiment@colData,
    metadata = experiment@metadata
  ))
}

get(load(opt$input))$rnaseq->experimentSummary

experimentSummary$batch <- droplevels(experimentSummary$batch)


batch_corrected<-correct_batch_effect(experiment = experimentSummary, covariate= as.formula(opt$covariate), method='ComBat')

if( !is.na(opt$tsv_corrected_counts) ) {
  write.table(cbind(`Gene ID`=rownames(assay(batch_corrected)),assay(batch_corrected)), file = opt$tsv_corrected_counts, sep = "\t", quote = FALSE, row.names = FALSE)
}

save( batch_corrected, file = opt$output )

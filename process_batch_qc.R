#!/usr/bin/env Rscript

##
## QC metrics, RUV and comparison with MBC cohort 
## Usage
## Usage:
##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gridExtra))


## ARGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
version="1.0"


parser <- ArgumentParser()
parser$add_argument("-i", help="3_IGG_SUBTRACTED.tsv", dest="input_file")
parser$add_argument("--validation_file", type="character", 
                    default="/Volumes/Histopathology\ Shared\ Resource/CLINICAL/Nanostring/REFERENCE_FILES/validation_samples_iggsub_normalized.txt",
                    dest="validation_file", help="validation file for controls comparison")
parser$add_argument("--mbc_md_file", type="character", default= "/Volumes/Histopathology\ Shared\ Resource/CLINICAL/Nanostring/REFERENCE_FILES/validation_mbc_metadata.txt",
                    dest="mbc_md_file", help="metastatic breast cancer metadata file")
parser$add_argument("--ab_ref_file", type="character", default= "/Volumes/Histopathology\ Shared\ Resource/CLINICAL/Nanostring/REFERENCE_FILES/ANTIBODY_REFERENCE.csv",
                    dest="ab_ref_file", help="ANTIBODY_REFERENCE.csv")
parser$add_argument("--include_bad_Ab", action="store_true", default=FALSE,
                    dest="include_bad_Ab", help="include all antibodies")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()
includeBCCL = args$includeBCCL

input_file = args$input_file
validation_file = args$validation_file
mbc_md_file = args$mbc_md_file
ab_ref_file = args$ab_ref_file
include_bad_Ab = args$include_bad_Ab


## TEST ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!grep("IGG_SUBTRACTED", input_file)){
  print("Error: The input file is not #3: IGG_SUBTRACTED.tsv")
  stop()
} else if (!grep("iggsub", basename(validation_file))){
  print("Error: The validation file is not #3: IGG_SUBTRACTED.tsv")
  stop()
}

    

#Things to include and exclue
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468+EGF")
ctrlregex = gsub("\\+", "\\\\+", paste0(paste0(controls, collapse = "|"),"|", "MDA468"))
## FUNCTIONS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## INPUT DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (include_bad_Ab){
  print("keeping all antibodies: IgG and antibodies that did not perform well or did not have dynamic range")
  #validation file
  validation = read.csv(validation_file, sep= "\t", row.names = 1, check.names = F)
  colnames(validation)=make.names(gsub("\\.1$", "", colnames(validation)))
  validation=log2(validation+1)
  #AB FILE
  ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)
  #NEW BATCH
  new_batch = data.matrix(read.table(file = input_file, sep="\t", row.names=1, stringsAsFactors=F, header=T))
  #MBC 
  mbc_md= read.table(file = mbc_md_file, sep="\t", row.names=1, stringsAsFactors=F, header=T)
  
} else {
  omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2")
  omitregex = paste0(paste0("^", omit), collapse = "|")
  #validation file
  validation = read.csv(validation_file, sep= "\t", row.names = 1, check.names = F)
  validation = validation[!grepl(omitregex, rownames(validation)),]
  
  lvalidation=log2(validation+1)
  colnames(lvalidation)=make.names(gsub("\\.1$", "", colnames(lvalidation)))
  
  #AB FILE
  ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)
  ab_ref = ab_ref[!grepl(omitregex, ab_ref$X.AbID),]
  #NEW BATCH
  new_batch = data.matrix(read.table(file = input_file, sep="\t", row.names=1, stringsAsFactors=F, header=T))
  new_batch = new_batch[!grepl(omitregex, rownames(new_batch)),]
  #MBC 
  mbc_md= read.table(file = mbc_md_file, sep="\t", row.names=1, stringsAsFactors=F, header=T)
}

#AB_ORDER
ab_order = ab_ref$X.AbID[order(ab_ref$Target)]
ab_order = ab_order[!grepl(omitregex, ab_order)]

#VALIDATION QC MATRIX #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lctl_norm = as.matrix(lvalidation)[,grepl(ctrlregex, colnames(lvalidation))]
lctl_norm.m = melt(lctl_norm)
lctl_norm.m$batch = sapply(strsplit(as.character(lctl_norm.m$Var2), "__"), `[`, 1)
lctl_norm.m$ctlname = sapply(strsplit(as.character(lctl_norm.m$Var2), "__"), `[`, 2)
lctl_norm.m$ctl_probe = factor(paste0(lctl_norm.m$ctlname,"_", lctl_norm.m$Var1))


validation_stats = data.frame(
  #do.call(rbind, (tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, summary))),  #quartiles
  "min" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, min),
  "mean" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean),
  "max" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, max),
  "stddev" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd),
  "coeff_var" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd)/tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean),
  "mean.minus.2sd" =tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean) - 
    2*(tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd)),
  "mean.plus.2sd" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean) + 
    2*(tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd)))

#only use controls 
ctl_newbatch = melt(log2(new_batch[,make.names(controls)]+1))
ctl_newbatch$ctl_probe = paste0(gsub("\\.1$", "", sapply(strsplit(as.character(ctl_newbatch$Var2), split = "_"), tail,1)),"_", ctl_newbatch$Var1)
ctl_newbatch$cellline = factor(gsub("\\.1$", "", sapply(strsplit(as.character(ctl_newbatch$Var2), split = "_"), tail,1)))

ctl_validation =cbind(validation_stats, "new_batch"=ctl_newbatch[match(rownames(validation_stats),ctl_newbatch$ctl_probe), "value"])
ctl_validation$status = ifelse(ctl_validation$new_batch < ctl_validation$mean.minus.2sd, 
                               (ctl_validation$mean - ctl_validation$new_batch)/ctl_validation$stddev, 
                               ifelse(ctl_validation$new_batch >= ctl_validation$mean.plus.2sd, 
                                      (ctl_validation$new_batch - ctl_validation$mean)/ctl_validation$stddev, "PASS"))

write.table(ctl_validation, file = "qc_controls.tsv", 
            sep="\t", quote = F, row.names = T, col.names = NA)

## QC LINEAR MODEL PLOTS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Comparing newQC with old QC 
#simple linear regression analysis

lmfitqc = lm(new_batch ~ mean, data=ctl_validation)

#~ plot linear model
ggReg = ggplot(lmfitqc$model, aes_string(x = names(lmfitqc$model)[2], y = names(lmfitqc$model)[1])) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Antibody Agreement: point per cell line per antibody",
                     "\nCorrelation= ",signif(cor((lmfitqc$model)[1],(lmfitqc$model)[2]), 5),
                     "\nAdj R2= ",signif(summary(lmfitqc)$adj.r.squared, 5),
                     "\nIntercept=",signif(lmfitqc$coef[[1]],5),
                     " p=",signif(summary(lmfitqc)$coef[1,4], 5),
                     "\nSlope=",signif(lmfitqc$coef[[2]], 5),
                     " p=",signif(summary(lmfitqc)$coef[2,4], 5)), 
       x="mean of validation samples")

twosd = 2*sd(ctl_validation$new_batch-ctl_validation$mean)
#outdat =  ctl_validation[which((ctl_validation$new_batch-ctl_validation$mean) < -twosd|(ctl_validation$new_batch-ctl_validation$mean) > twosd),]
outdat=ctl_validation
outdat$cellline = sapply(strsplit(as.character(rownames(outdat)), "_"), `[`, 1)
outdat$probe = sapply(strsplit(as.character(rownames(outdat)), "_"), `[`, 2)

baplot_ctl = ggplot(outdat) +
  geom_point(aes(x=mean, y=new_batch-mean, color=cellline)) +
  geom_hline(yintercept=c(0, twosd, -twosd), color="red", linetype = 2) +
  geom_text(data=data.frame(x=0,y=c(-twosd, twosd)), aes(x, y), label = c("2SD", "-2SD"), color="red") +
  #geom_point(data = outdat, aes(x=mean,  y=outdat$new_batch-outdat$mean,color=outdat$cellline)) +
  #geom_text(data = outdat, aes(x=mean,  y=outdat$new_batch-outdat$mean,label=rownames(outdat), hjust=-0.1), size=3) +
  labs(title = "Mean-Difference Plot of Celllines", x="Mean of validation batches", y="New Batch - Mean of validation batches")

baplot_ab = ggplot(outdat) +
  geom_point(aes(x=mean, y=new_batch-mean, color=probe)) +
  geom_hline(yintercept=c(0, twosd, -twosd), color="red", linetype = 2) +
  geom_text(data=data.frame(x=0,y=c(-twosd, twosd)), aes(x, y), label = c("2SD", "-2SD"), color="red") +
  #geom_point(data = outdat, aes(x=mean,  y=outdat$new_batch-outdat$mean,color=outdat$probe)) +
  labs(title = "Mean-Difference Plot of Antibodies", x="Mean of validation batches", y="New Batch - Mean of validation batches") +
  theme(legend.text=element_text(size=4))


print("saving QC plot")
ggsave(filename ="QC_antibody_linear_plot.pdf", device = "pdf", plot=ggReg,width = 8, height = 8)
ggsave(filename ="QC_antibody_meandiff_plot.pdf", device = "pdf", grid.arrange(baplot_ab, baplot_ctl),width = 8, height = 5)


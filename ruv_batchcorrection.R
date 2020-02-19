#!/usr/bin/env Rscript

##
## QC metrics, RUV and comparison with MBC cohort 
## Usage
## ./ruv_mbc.R -i 3_IGG_NORMALIZED.tsv --validation_file validation_samples_normalized.txt --mbc_md_file validation_mbc_metadata.txt --ab_ref_file ANTIBODY_REFERENCE.csv
##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ruv))


## ARGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
version="1.0"


parser <- ArgumentParser()
parser$add_argument("-i", help="", dest="input_file")
parser$add_argument("--validation_file", type="character", 
                    default="/Volumes/Histopathology\ Shared\ Resource/CLINICAL/Nanostring/REFERENCE_FILES/validation_samples_geosamp_normalized.txt",
                    dest="validation_file", help="validation file for controls comparison")
parser$add_argument("--mbc_md_file", type="character", default= "/Volumes/Histopathology\ Shared\ Resource/CLINICAL/Nanostring/REFERENCE_FILES/validation_mbc_metadata.txt",
                    dest="mbc_md_file", help="metastatic breast cancer metadata file")
parser$add_argument("--ab_ref_file", type="character", default= "/Volumes/Histopathology\ Shared\ Resource/CLINICAL/Nanostring/REFERENCE_FILES/ANTIBODY_REFERENCE.csv",
                    dest="ab_ref_file", help="ANTIBODY_REFERENCE.csv")
parser$add_argument("--include_ctrls", action="store_true", default=FALSE,
                    dest="include_ctrls", help="include all antibodies")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()
includeBCCL = args$includeBCCL

input_file = args$input_file
validation_file = args$validation_file
mbc_md_file = args$mbc_md_file
ab_ref_file = args$ab_ref_file
include_ctrls = args$include_ctrls


## TEST ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!grepl("GEOMEAN_NORM", input_file)){
  print("Error: The input file is not #2: 2_GEOMEAN_NORMALIZED.tsv")
  stop()
} else if (!grepl("geosamp", basename(validation_file))){
  print("Error: The validation file is not #2:validation_samples_geosamp_normalized.txt")
  stop()
}



## FUNCTIONS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## OUTPUT PREP ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir.create("ruv_figures", showWarnings = F)

## INPUT DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#validation file
validation = read.csv(validation_file, sep= "\t", row.names = 1, check.names = F)
colnames(validation)=make.names(gsub("\\.1$", "", colnames(validation)))
#validation=log2(validation+1)
#AB FILE
ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)
#NEW BATCH
new_batch = data.matrix(read.table(file = input_file, sep="\t", row.names=1, stringsAsFactors=F, header=T))
#MBC 
mbc_md= read.table(file = mbc_md_file, sep="\t", row.names=1, stringsAsFactors=F, header=T)



#Things to include and exclude
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468+EGF")
ab.ctrl = make.names(rownames(new_batch)[grepl("NEG|IgG|^S6$|^Histone", rownames(new_batch))])


## RUV DATASET ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#empty data frame
mbc_percentile=data.frame(row.names = rownames(validation))

for (samp in setdiff(colnames(new_batch), make.names(controls))) {
  print(samp)
  # COMBINING: rows combined with mbcctl, columns=ctrls and 1 samp only
  idx_controls = which(colnames(new_batch) %in% make.names(controls))
  newsampctl = new_batch[,c(idx_controls,which(colnames(new_batch)==samp))]
  colnames(newsampctl) = paste0("newbatch", "__", colnames(newsampctl))
  combined_norm = cbind(validation, newsampctl[match(rownames(validation), rownames(newsampctl)),])

  #make metadata to get replicate matrix
  combined_md = data.frame(batch = sapply(strsplit(as.character(colnames(combined_norm)), "__"), `[`, 1), 
                           samp =sapply(strsplit(as.character(colnames(combined_norm)), "__"), `[`, 2),
                           sampcolumn = c(colnames(combined_norm)),stringsAsFactors = F)
  combined_md$mbc= "notMBC"
  combined_md[combined_md$batch %in% make.names(mbc_md$BatchID) & combined_md$samp %in% make.names(controls),"mbc"] = "rep_control"
  combined_md[combined_md$batch %in% make.names(mbc_md$BatchID) & combined_md$samp %in% make.names(mbc_md$Sample),"mbc"] = "MBCsample"
  combined_md[combined_md$batch =="newbatch" & combined_md$samp %in% make.names(controls),"mbc"] = "rep_control"
  combined_md[combined_md$batch =="newbatch" & !combined_md$samp %in% make.names(controls),"mbc"] = "sample"

  # REPLICATE MATRIX: rep matrix based only on controls in MBC batches
  combined_md$reps=ifelse(combined_md$mbc == "rep_control", yes=combined_md$samp, no=combined_md$sampcolumn)
  repmat = replicate.matrix(combined_md$reps)
  ctlmat = repmat[,make.names(controls)]

  # PLOTTING ~ plotting pre-ruv figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  lctls = log2(as.matrix(combined_norm)[,which(rowSums(ctlmat)==1)]+1)
  limits_rle = max(as.matrix(lctls))-median(as.matrix(lctls))

  #rle 
  rle_orig = ruv_rle(Y = t(lctls), rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% colnames(lctls),]), 
                     ylim=c(-limits_rle,limits_rle)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label=samp), angle=90, hjust=0, size=2)+
    theme(legend.position = "right", legend.text = element_text(size=6)) + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression \nNormalized Data (IGG corrected) \n",  samp))
  #pca
  pca_orig = ruv_svdplot(Y.data = t(lctls), info = combined_md$samp[which(rowSums(ctlmat)==1)])
  
  clust = plot(hclust(dist(t(lctls), method ="euclidean")))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #~ RUV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  RUVcorrected = RUVIII(Y=t(log2(combined_norm +1)), ctl=make.names(rownames(combined_norm)) %in% ab.ctrl, M=repmat, k=1)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~ PLOTTING RUV processed figures~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ctl_ruv = t(RUVcorrected)[,(rownames(RUVcorrected) %in% colnames(lctls))]
  limits_rle = max(as.matrix(ctl_ruv))-median(as.matrix(ctl_ruv))
  
  rle_ruv = ruv_rle(Y =t(ctl_ruv), 
                    rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% colnames(ctl_ruv),]), 
                    ylim=c(-limits_rle,limits_rle)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label=samp), angle=90, hjust=0, size=2) +
    theme(legend.position = "right") + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression RUV processed\n", samp))

  pca_ruv = ruv_svdplot(Y.data=t(ctl_ruv), info=combined_md[combined_md$sampcolumn %in% colnames(ctl_ruv),])
  clust_ruv = plot(hclust(dist(t(ctl_ruv), method ="euclidean")))
  
  #~ MBC DATA PREP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (include_ctrls){
    print("keeping all antibodies: IgG and antibodies that did not perform well or did not have dynamic range")
    truv = t(RUVcorrected)
    other_abs=setdiff(colnames(RUVcorrected),ab_ref$X.AbID)
    ab_order = c(ab_ref$X.AbID[order(ab_ref$Target)], other_abs)

  } else {
    omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2", "NEG")
    omitregex = paste0(paste0("^", omit), collapse = "|")
    truv = t(RUVcorrected[,!grepl(omitregex, colnames(RUVcorrected))])
    #AB_ORDER
    ab_order = ab_ref$X.AbID[order(ab_ref$Target)]
    ab_order = ab_order[!grepl(omitregex, ab_order)]
    
  }
  
  #~ split these apart makes plotting easier
  #mbc only and newsamp only
  mbc_ruv = truv[,as.character(combined_md$sampcolumn[combined_md$mbc=="MBCsample"])]
  new_ruv = truv[,as.character(combined_md$sampcolumn[combined_md$mbc=="sample"]),drop=F]
  
  #~ MELTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mbcruv.m = melt(mbc_ruv,  id.vars=row.names)
  newruv.m = melt(new_ruv, id.vars=row.names)
  mbc_ecdf = tapply(mbcruv.m$value, mbcruv.m$Var1, ecdf)
  
  mat_newruv.m = as.matrix(newruv.m)
  
  for (idx in seq(1:nrow(newruv.m))){
    new_samp=mat_newruv.m[idx,]
    #print(new_samp)
    mbc_percentile[new_samp[1],samp] = mbc_ecdf[[new_samp[1]]](new_samp[3])
  }
#  mbc_percentile[,samp] =  as.vector(apply(newruv.m, 1, function(x) 
#    ((mbc_ecdf[[as.character(x["Var1"])]](x[["value"]]))))[rownames(mbc_percentile)])
  
  #max and min
  ruv = as.matrix(truv)
  mruv = melt(ruv)
  
  ruv_stats = data.frame(
    #do.call(rbind, (tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, summary))),  #quartiles
    "min" = tapply(mruv$value, mruv$Var1, min),
    "q1" = tapply(mruv$value, mruv$Var1, function(x) quantile(x, 0.25)),
    "mean" = tapply(mruv$value, mruv$Var1, median),
    "q3" = tapply(mruv$value, mruv$Var1, function(x) quantile(x, 0.75)),
    "max" = tapply(mruv$value, mruv$Var1, max))
  
  mbcruv.m$Var1 = factor(mbcruv.m$Var1, levels=ab_order)
  newruv.m$Var1 = factor(newruv.m$Var1, levels=ab_order)
  ruv_stats$ab = factor(rownames(ruv_stats), levels=ab_order)
  
  mbc_labels= paste0(as.character(levels(newruv.m$Var1))," (",round(newruv.m$value, 1), ",",round(mbc_percentile[as.character(levels(newruv.m$Var1)), samp]*100,0),")")
  
  
  #PLOTTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #~ boxplot of MBC
  bp_mbcruv = ggplot(data.frame(mbcruv.m)) +
    geom_point(data=newruv.m, mapping=aes(x=factor(Var1), y=value), colour=c('red'), shape=8, size=2) +
    geom_linerange(
      data=ruv_stats, aes(x=ab, ymin = min, ymax = max),
      color = "#808080", 
      size = 7, 
      alpha = 0.7) +
    #geom_point(data = ruv_stats, aes(x=ab, y=max),shape=93, fill="grey") +
    geom_boxplot(aes(factor(Var1), as.numeric(value)), outlier.colour = NA) +
    geom_point(data=newruv.m, mapping=aes(x=factor(Var1), y=value), colour=c('red'), shape=8, size=2) +
    scale_x_discrete(labels=paste0(as.character(levels(newruv.m$Var1))," (",
                                   round(mbc_percentile[as.character(levels(newruv.m$Var1)), samp]*100,0),")")) +
    labs(x="Antibody (Percentiles)", title=paste0(samp,  "\n within Distribution of Metastatic Breast Cancers"), y="RUVnormalized") +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5, vjust=0),
          legend.text=element_text(size=8),
          legend.position="right",
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black")) +
    coord_flip()
  
  print("saving RUV plots")
  #~ save all plots
  ggsave(file=paste0(samp, "_MBC.pdf"), bp_mbcruv, device="pdf", width = 8, height = 6)
  ggsave(file=paste0("ruv_figures/",samp, "_controls_norm_RLE.pdf"), 
         device="pdf", rle_orig, width = 9, height = 4.5)
  ggsave(file=paste0("ruv_figures/",samp, "_controls_norm_PCA.pdf"), 
         device="pdf", pca_orig, width = 8, height = 7)
  ggsave(file=paste0("ruv_figures/",samp, "_controls_RUV_RLE.pdf"), 
         device="pdf", rle_ruv, width = 9, height = 4.5)
  ggsave(file=paste0("ruv_figures/",samp, "_controls_RUV_PCA.pdf"), 
         device="pdf", pca_ruv, width = 8, height = 7)
  
  
  write.table(x = t(RUVcorrected), file=paste0("ruv_figures/",samp, "_RUVcorrected.tsv"), 
              sep="\t", quote = F, row.names = T, col.names = NA)
}

tar(tarfile=paste0("ruv_figures.tar.gz"), files=paste0("ruv_figures"), compression="gzip", tar="tar")
write.table(x = round(mbc_percentile, 2), file=paste0("MBC_percentiles.tsv"), sep="\t", quote = F, row.names = T, col.names = NA)







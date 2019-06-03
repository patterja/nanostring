#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

usage = "
./ruv_mbc.R input_NORMALIZED.xlsx validation_file validation_md

input_file= NORMALIZED.xlsx. Will look for igg_geosamp_corrected sheet. 
validation_dataset = txt file of validation data normalized
validation_metadata =xlsx with metastatic breast cancer list. 
"

argsLen <- length(args);
if (argsLen == 3) {
  input_file = args[1]
  validation_file = args[2]
  validation_md = args[3]
 
  print(paste0("Processing ", input_file))
} else if (argsLen == 1) {
  input_file = args[1]
  validation_file = "/Users/patterja/Box Sync/NANOSTRING/REFERENCE_FILES/ALLvalidation_samples_corrected.txt"
  validation_md = "/Users/patterja/Box Sync/NANOSTRING/VALIDATION/ALL_validation_samples.xlsx"
} else {
  stop(cat(usage))
}



library(xlsx)
library(ruv)
library(ggplot2)
library(reshape2)

pcaplt <- function (mat, title = "PCA Plot", repdf) {
  var = mat[apply(mat, 1, var, na.rm = TRUE) != 0, ]
  cc.var = var[complete.cases(var), ]
  pca_prcomp = prcomp(t(var), center = T, scale = F)
  PC1_and_PC2 = data.frame(PC1 = pca_prcomp$x[, 1], PC2 = pca_prcomp$x[,2], type = rownames(pca_prcomp$x))
  perc = (pca_prcomp$sdev^2)/sum(pca_prcomp$sdev^2) * 100
  labs <- sapply(seq_along(perc), function(i) {
    paste("PC ", i, " (", round(perc[i], 2), "%)", sep = "")})
  
  PCsmd = cbind(PC1_and_PC2, repdf[match(rownames(PC1_and_PC2), repdf$sampcolumn), c("samp", "batch")])
  p = ggplot(PCsmd,aes_string("PC1", "PC2", col="batch")) + 
    geom_point(size = 1.5) + 
    geom_text(aes(label = type), vjust = -1, size=2) + 
    labs(title = title,x = labs[1], y = labs[2]) + 
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "gray"), plot.title = element_text(hjust = 0.5), 
          legend.text = element_text(size = 4), legend.position = "right")
  return(p)
}
#Things to include and exclue
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468+EGF")
ctrlregex = gsub("\\+", "\\\\+", paste0(paste0(controls, collapse = "|"),"|", "MDA468"))
omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2")
omitregex = paste0(paste0("^", omit), collapse = "|")


#VALIDATION DATA
#validation_md = "/Users/patterja/Box Sync/NANOSTRING/VALIDATION/ALL_validation_samples.xlsx"
#validation_file = "/Users/patterja/Box Sync/NANOSTRING/REFERENCE_FILES/ALLvalidation_samples_corrected.txt"

mbc_md = read.xlsx(file= validation_md, sheetName = "MBC")
igg_geosamp = read.csv(validation_file, sep= "\t", row.names = 1)
igg_geosamp = igg_geosamp[!grepl(omitregex, rownames(igg_geosamp)),]
#PROCESS VALIDATION DATA
ctl_norm = igg_geosamp[,grepl(ctrlregex, colnames(igg_geosamp))]
mbc_norm = igg_geosamp[,make.names(as.character(mbc_md$Combined_Name))]
lctl_norm.m = melt(as.matrix(log2(ctl_norm+1)))
lctl_norm.m$ctl_probe = paste0(gsub("\\.1$", "", sapply(strsplit(as.character(lctl_norm.m$Var2), split = "_"), tail,1)), 
                               "_", make.names(lctl_norm.m$Var1))
lctl_norm.m$cellline = factor(gsub("\\.1$", "", sapply(strsplit(as.character(lctl_norm.m$Var2), split = "_"), tail,1)))


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


#INPUT NEW FILE
#NEW BATCH
#input_file =  "/Users/patterja/Box Sync/NANOSTRING/data/20190530_208420541020_SMMART7/20190530_208420541020_ST-05302019-KD0001_OUTPUT/20190530_208420541020_ST-05302019-KD0001_NORMALIZED.xlsx"

new_batch = read.xlsx(file= input_file, sheetName = "igg_geosamp_corrected")
row.names(new_batch) = as.character(new_batch[,1])
new_batch = new_batch[,-1]


#QC OUPUT

ctl_newbatch = melt(as.matrix(log2(new_batch[,make.names(controls)]+1)))

ctl_newbatch$ctl_probe = paste0(gsub("\\.1$", "", sapply(strsplit(as.character(ctl_newbatch$Var2), split = "_"), tail,1)), 
                               "_", make.names(ctl_newbatch$Var1))
ctl_newbatch$cellline = factor(gsub("\\.1$", "", sapply(strsplit(as.character(ctl_newbatch$Var2), split = "_"), tail,1)))

ctl_validation =cbind(validation_stats, "new_batch"=ctl_newbatch[match(rownames(validation_stats), 
                                                                       ctl_newbatch$ctl_probe), "value"])
ctl_validation$status = ifelse(ctl_validation$new_batch < ctl_validation$mean.minus.2sd, 
                               ctl_validation$mean.minus.2sd - ctl_validation$new_batch, 
                               ifelse(ctl_validation$new_batch > ctl_validation$mean.plus.2sd, 
                                      ctl_validation$new_batch - ctl_validation$mean.plus.2sd, "PASS"))


batch_name = tail(unlist(strsplit(dirname(normalizePath(input_file)), "/")),2)[1]
write.table(ctl_validation, file = paste0(dirname(normalizePath(input_file)), "/", batch_name, "_QC_CONTROLS.csv"), 
            sep=",", quote = F, row.names = T, col.names = NA)


#~ DATASET FOR RUV

# - mbc filtered from metadata sheet only, excludes some samples in igg_geosamp
mbcctl = cbind(ctl_norm, mbc_norm[match(rownames(ctl_norm), rownames(mbc_norm)),])
dir.create(file.path(paste0(dirname(normalizePath(input_file)), "/ruv_figures"), showWarnings = F))



#for sample in newbatch
mbc_percentile=data.frame(row.names = rownames(igg_geosamp))

for (samp in setdiff(colnames(new_batch), make.names(controls))) {
  print(samp)
  #COMBINING: rows combined with mbcctl, columns=ctrls and 1 samp only
  newsampctl = new_batch[,c(grep(ctrlregex, colnames(new_batch)),which(colnames(new_batch)==samp))]
  combined_norm = cbind(mbcctl, newsampctl[match(rownames(mbcctl), rownames(newsampctl)),])
  
  
  #~ replicated matrix 
  repdf = data.frame(samp = c(gsub("\\.1$", "", sapply(strsplit(as.character(colnames(mbcctl)), split = "_"), tail,1)),
                              colnames(newsampctl)),
                     batch = c(unlist(lapply(as.character(colnames(mbcctl)), function(x) substring(x, 0,29))),
                               rep(batch_name, length(colnames(newsampctl)))), 
                     sampcolumn = c(colnames(mbcctl), colnames(newsampctl)))
  repmat = replicate.matrix(c(make.names(repdf[,1])))
  
  #~ plotting pre-ruv figures
  lctls = log2(combined_norm[,grep(ctrlregex, colnames(combined_norm))]+1)
  
  limits_rle = max(as.matrix(lctls))-median(as.matrix(lctls))
  rle_orig = ruv_rle(Y = t(lctls), rowinfo = as.matrix(repdf[repdf$sampcolumn %in% colnames(lctls),]), 
                ylim=c(-limits_rle,limits_rle)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label=samp), angle=90, hjust=0, size=2)+
    theme(legend.position = "right", legend.text = element_text(size=6)) + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression \nNormalized Data (IGG corrected) \n", batch_name, "_", samp))
  
  pca_orig = pcaplt(lctls, title = paste0("log2 (Normalized Signal + 1)\n", batch_name, "_", samp), repdf)
  
  #~ RUV
  RUVcorrected = RUVIII(Y=t(log2(combined_norm +1)), ctl=c(1:nrow(combined_norm)), M=repmat, k=1)
  
  #~ plotting RUV processed figures
  ctl_ruv = t(RUVcorrected)[,grep(ctrlregex, colnames(t(RUVcorrected)))]
  rle_ruv = ruv_rle(Y = RUVcorrected[grep(ctrlregex, rownames(RUVcorrected)),], 
                    rowinfo = as.matrix(repdf[repdf$sampcolumn %in% colnames(ctl_ruv),]), 
                    ylim=c(-limits_rle,limits_rle)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label=samp), angle=90, hjust=0, size=2) +
    theme(legend.position = "right") + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression RUV processed\n", batch_name, "_", samp))
  
  pca_ruv = pcaplt(ctl_ruv, title = paste0("RUV Processed\n", batch_name, "_", samp), repdf)
  
  #~ split these apart makes plotting easier
  #ctl_ruv already assigned before graphing
  mbc_ruv = t(RUVcorrected)[,make.names(as.character(mbc_md$Combined_Name))]
  new_ruv = data.frame("Var1" = rownames(t(RUVcorrected)),
                       "Var2"= samp,
                       "value" = t(RUVcorrected)[rownames(t(RUVcorrected)), samp])
  #~ Get percentile using ecdf
  mbcruv.m = melt(mbc_ruv,  id.vars=row.names)
  mbc_ecdf = tapply(mbcruv.m$value, mbcruv.m$Var1, ecdf)
  #mbc_percentile[,samp] = new_ruv[match(rownames(mbc_stats), new_ruv$Var1),"value"]
  mbc_percentile[,samp] = as.vector(apply(new_ruv, 1, function(x) 
    ((mbc_ecdf[[x["Var1"]]](x[["value"]]))))[rownames(mbc_percentile)])
  
  #~ melt for plotting
  mbcruv.m = melt(mbc_ruv,  id.vars=row.names)
  mbcruv.m$ctl_probe = paste0(gsub("\\.1$", "", sapply(strsplit(as.character(mbcruv.m$Var2), split = "_"), tail,1)), 
                                 "_", make.names(mbcruv.m$Var1))
  mbcruv.m$cellline = factor(gsub("\\.1$", "", sapply(strsplit(as.character(mbcruv.m$Var2), split = "_"), tail,1)))
  

  #~ boxplot of MBC
  bp_mbcruv =ggplot(data.frame(mbcruv.m), aes(Var1, as.numeric(value))) +
    geom_boxplot() +
    geom_point(data=new_ruv, mapping=aes(y=value),  colour=c("red"), shape=8, size=2) +
    coord_flip() +
    labs(x="Antibody", title=paste0(samp,  "\n within Distribution of Metastatic Breast Cancers"), y="RUVnormalized") +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5, vjust=0),
          legend.text=element_text(size=8),
          legend.position="right",
          axis.text.x = element_text(hjust = 1, size=7, colour="black"),
          axis.text.y = element_text(size=6, colour="black"))
  
  #~ save all plots
  ggsave(file=paste0(dirname(normalizePath(input_file)), "/",samp, "_", batch_name, "_MBC.png"), bp_mbcruv, width = 8, height = 4)
  ggsave(file=paste0(dirname(normalizePath(input_file)), "/ruv_figures/",samp, "_", batch_name, "_controls_norm_RLE.png"), 
         rle_orig, width = 9, height = 4.5)
  ggsave(file=paste0(dirname(normalizePath(input_file)), "/ruv_figures/",samp, "_", batch_name, "_controls_norm_PCA.png"), 
         pca_orig, width = 8, height = 7)
  ggsave(file=paste0(dirname(normalizePath(input_file)), "/ruv_figures/",samp, "_", batch_name, "_controls_RUV_RLE.png"), 
         rle_ruv, width = 9, height = 4.5)
  ggsave(file=paste0(dirname(normalizePath(input_file)), "/ruv_figures/",samp, "_", batch_name, "_controls_RUV_PCA.png"), 
         pca_ruv, width = 8, height = 7)
  
  
  write.table(x = t(RUVcorrected), file=paste0(dirname(normalizePath(input_file)), "/ruv_figures/",samp, "_", batch_name, "_RUVcorrected.csv"), 
              sep=",", quote = F, row.names = T, col.names = NA)
}


write.table(x = round(mbc_percentile, 2), file=paste0(dirname(normalizePath(input_file)), "/", batch_name, "_MBCpercentile.csv"), sep=",", quote = F, row.names = T, col.names = NA)


  
  
  
  

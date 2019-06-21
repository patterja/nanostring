#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

usage = "
./ruv_batch_correction.R /Users/patterja/Box\ Sync/NANOSTRING/data/ ruv_sheet.txt validation_file validation_md

data_dir = where all your data lives
ruv_sheet = txt file samples and normalized matrices
validation_file = txt file of validation data normalized
validation_md =xlsx with metastatic breast cancer list. 
"


argsLen <- length(args);
if (argsLen ==1){
  cat("need ruv sheet")
  stop(cat(usage))
  
} else if (argsLen ==2) {
  data_dir = args[1]
  ruv_sheet = args[2]
  validation_file = "/Users/patterja/Box Sync/NANOSTRING/REFERENCE_FILES/ALL_validation_samples_normalized.txt"
  validation_md ="/Users/patterja/Box Sync/NANOSTRING/REFERENCE_FILES/ALL_validation_samples.xlsx"
  print(paste0("Processing ", ruv_sheet))
} else {
  stop(cat(usage))
}

library(xlsx)
library(ruv)
library(ggplot2)
library(reshape2)
library(pheatmap)

#Things to include and exclue
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468+EGF")
ctrlregex = gsub("\\+", "\\\\+", paste0(paste0(controls, collapse = "|"),"|", "MDA468"))
omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2")
omitregex = paste0(paste0("^", omit), collapse = "|")


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


samps2ruv = read.csv(file= ruv_sheet, sep="\t", stringsAsFactors = F)
files = as.character(unique(samps2ruv$Normalized.file.name))
sampname = gsub("_samplesheet.txt", "", basename(ruv_sheet))

allbatch_norms = data.frame()

for (f in 1:length(files)){
  file_name = list.files(path = data_dir,pattern = as.character(files[f]), recursive = T)
  batch_name = unlist(strsplit(file_name, "/"))[1]
  #get normalized data
  batch = read.xlsx2(file= paste0(data_dir, file_name), sheetName = "igg_geosamp_corrected",stringsAsFactors=F)
  rownames(batch) = batch[,1]
  batch[,1] = NULL
  #get only ruv samples
  sel_batch = batch[,c(grep(ctrlregex, colnames(batch)), which(colnames(batch) %in% make.names(samps2ruv$Sample.Name)))]
  md = data.frame("batch"=batch_name, "sample" = colnames(sel_batch), "fullname" = paste0(batch_name,"_", colnames(sel_batch)))
  colnames(sel_batch) = paste0(batch_name,"_", colnames(sel_batch))
  
  if (f==1){
    allbatch_norms = data.matrix(sel_batch)
    md_allbatch_norms = md
  }
  else{
    allbatch_norms = data.matrix(cbind(allbatch_norms, sel_batch[match(rownames(allbatch_norms), rownames(sel_batch)),]))
    md_allbatch_norms = rbind(md_allbatch_norms, md)
  }
}


#VALIDATION DATA

#validation_file = "/Users/patterja/Box Sync/NANOSTRING/REFERENCE_FILES/ALL_validation_samples_normalized.txt"
#validation_md = "/Users/patterja/Box Sync/NANOSTRING/REFERENCE_FILES/ALL_validation_samples.xlsx"
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

#MBC
mbc_md = read.xlsx(file= validation_md, sheetName = "MBC")
igg_geosamp = read.csv(validation_file, sep= "\t", row.names = 1)
igg_geosamp = igg_geosamp[!grepl(omitregex, rownames(igg_geosamp)),]

# - mbc filtered from metadata sheet only, excludes some samples in igg_geosamp
mbcctl = cbind(ctl_norm, mbc_norm[match(rownames(ctl_norm), rownames(mbc_norm)),])

#~ DATASET FOR RUV

#~ combining: rows combined with mbcctl, ctrls and samples
combined_mbcnorm = data.matrix(cbind(mbcctl, allbatch_norms[match(rownames(mbcctl), rownames(allbatch_norms)),]))

#~ replicated matrix
repdf = data.frame(samp = c(gsub("\\.1$", "", sapply(strsplit(as.character(colnames(combined_mbcnorm)), split = "_"), tail,1))),
                   batch = c(unlist(lapply(as.character(colnames(combined_mbcnorm)), function(x) substring(x, 0,30)))),
                   sampcolumn = c(colnames(combined_mbcnorm)))
repmat = replicate.matrix(c(make.names(repdf[,1])))

#~ box plot
lctls = log2(combined_mbcnorm[,grep(ctrlregex, colnames(combined_mbcnorm))]+1)
limits_rle = max(as.matrix(lctls))-median(as.matrix(lctls))
rle_orig = ruv_rle(Y = t(lctls), rowinfo = as.matrix(repdf[repdf$sampcolumn %in% colnames(lctls),]), 
                   ylim=c(-limits_rle,limits_rle)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  #geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label="PRE-RUV"), angle=90, hjust=0, size=2)+
  theme(legend.position = "right", legend.text = element_text(size=6)) + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression \nNormalized Data (IGG corrected) \n", sampname, "_", "PRE-RUV"))

pca_orig = pcaplt(lctls, title = paste0("log2 (Normalized Signal + 1)\n", sampname, "_", "PRE-RUV"), repdf)

#~ RUV
RUVcorrected = RUVIII(Y=t(log2(combined_mbcnorm +1)), ctl=c(1:nrow(combined_mbcnorm)), M=repmat, k=1)

#~ plotting RUV processed figures
ctl_ruv = t(RUVcorrected)[,grep(ctrlregex, colnames(t(RUVcorrected)))]
rle_ruv = ruv_rle(Y = RUVcorrected[grep(ctrlregex, rownames(RUVcorrected)),], 
                  rowinfo = as.matrix(repdf[repdf$sampcolumn %in% colnames(ctl_ruv),]), 
                  ylim=c(-limits_rle,limits_rle)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  #geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label="RUVcorrected"), angle=90, hjust=0, size=2) +
  theme(legend.position = "right") + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression RUV processed\n", sampname, "_", "RUVcorrected"))

pca_ruv = pcaplt(ctl_ruv, title = paste0("RUV Processed\n", sampname, "_", "RUVcorrected"), repdf)


#~ split these apart makes plotting easier
#ctl_ruv already assigned before graphing
mbc_ruv = t(RUVcorrected)[,make.names(as.character(mbc_md$Combined_Name))]
samps_ruv = t(RUVcorrected)[,setdiff(rownames(RUVcorrected), make.names(as.character(colnames(igg_geosamp))))]
samps_ruv = samps_ruv[,!grepl(ctrlregex, colnames(samps_ruv))]

#~ Get percentile using ecdf
mbcruv.m = melt(mbc_ruv,  id.vars=row.names)
mbc_ecdf = tapply(mbcruv.m$value, mbcruv.m$Var1, ecdf)

mbc_percentile=data.frame(row.names = rownames(igg_geosamp))

for (samp in colnames(samps_ruv)){
  samp_df = melt(samps_ruv[,samp, drop=FALSE])
  rownames(samp_df) = samp_df$Var1
  mbc_percentile[,samp] = as.vector(apply(samp_df, 1, function(x) 
    ((mbc_ecdf[[x["Var1"]]](x[["value"]]))))[rownames(mbc_percentile)])
}



#~ melt for boxplot plotting
mbcruv.m = melt(mbc_ruv,  id.vars=row.names)
sampsruv.m = melt(samps_ruv, id.vars=row.names)
mbcruv.m$ctl_probe = paste0(gsub("\\.1$", "", sapply(strsplit(as.character(mbcruv.m$Var2), split = "_"), tail,1)), 
                            "_", make.names(mbcruv.m$Var1))
mbcruv.m$cellline = factor(gsub("\\.1$", "", sapply(strsplit(as.character(mbcruv.m$Var2), split = "_"), tail,1)))

#~ boxplot of MBC
bp_mbcruv = ggplot(data.frame(mbcruv.m), aes(Var1, as.numeric(value))) +
  geom_boxplot() +
  geom_point(data=sampsruv.m, mapping=aes(y=value, colour=Var2),shape=8, size=2) +
  coord_flip() +
  labs(x="Antibody", title=paste0(samp,  "\n within Distribution of Metastatic Breast Cancers"), y="RUVnormalized") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major=element_line(colour="gray"),
        plot.title = element_text(hjust = 0.5, vjust=0),
        legend.text=element_text(size=8),
        legend.position="bottom",
        axis.text.x = element_text(hjust = 1, size=7, colour="black"),
        axis.text.y = element_text(size=6, colour="black")) +
  guides(col = guide_legend(ncol = 1))


o_samps = md_allbatch_norms[match(make.names(samps2ruv$Sample.Name), md_allbatch_norms$sample), "fullname"]
#~ heatmap

png(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_heatmap.png"), width = 350, height = 800)
pheatmap(mat = samps_ruv[,as.character(o_samps)],
                color             = colorRampPalette( c("green", "black", "red"), space="rgb")(10),
                cluster_cols = T,
                show_colnames     = TRUE,
                show_rownames     = TRUE,
                fontsize          = 9,
                #display_numbers = T,number_format = "%.2f", number_color="white",
                main              = sampname)
dev.off()

# ~ functions for pairs plot
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
scatter_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,...)
}

plt_samps = samps_ruv
#colnames(plt_samps) = gsub("^.*1020_", "", colnames(plt_samps))
png(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_controls_RUV_scatter.png"), width = 850, height = 850)
plt_scatter = pairs(plt_samps, lower.panel = scatter_line, upper.panel = panel.cor, pch = 19)
dev.off()


#~ save all plots
ggsave(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_MBC.png"), bp_mbcruv, width = 8, height = 6)

ggsave(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_controls_norm_RLE.png"), 
       rle_orig, width = 5, height = 4.5)
ggsave(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_controls_norm_PCA.png"), 
       pca_orig, width = 8, height = 7)
ggsave(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_controls_RUV_RLE.png"), 
       rle_ruv, width = 9, height = 4.5)
ggsave(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_controls_RUV_PCA.png"), 
       pca_ruv, width = 8, height = 7)
ggsave(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_controls_RUV_scatter.png"), 
       plt_scatter, width = 8, height = 7)

write.table(x = t(RUVcorrected), file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_RUVcorrected.csv"), 
            sep=",", quote = F, row.names = T, col.names = NA)

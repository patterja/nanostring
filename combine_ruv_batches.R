#!/usr/bin/env Rscript

##
## Usage
## ./ruv_batch_correction.R /Users/patterja/Box\ Sync/NANOSTRING/data/ ruv_sheet.txt validation_file mbc_md_file
## data_dir = where all your data lives
## ruv_sheet = txt file samples and normalized matrices
## validation_file = txt file of validation data normalized
## mbc_md_file =txt with metastatic breast cancer list. 
##

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(ruv))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))

version="2.0"


parser <- ArgumentParser()
parser$add_argument("--data_dir", type="character", default="/Users/patterja/Box\ Sync/NANOSTRING/output_subtracted_newfig", 
                    dest="data_dir", help="directory name")
parser$add_argument("--ruv_sheet", type="character", dest="ruv_sheet", help="ruv samplesheet")
parser$add_argument("--validation_file", type="character", default="/Users/patterja/Box\ Sync/NANOSTRING/REFERENCE_FILES/validation_samples_sub_normalized.txt",
                    dest="validation_file", help="validation file for controls comparison")
parser$add_argument("--mbc_md_file", type="character", default= "/Users/patterja/Box Sync/NANOSTRING/REFERENCE_FILES/validation_mbc_metadata.txt",
                    dest="mbc_md_file", help="metastatic breast cancer metadata file")
parser$add_argument("--ab_ref_file", type="character", default= "/Users/patterja/Box\ Sync/NANOSTRING/REFERENCE_FILES/ANTIBODY_REFERENCE.csv",
                    dest="ab_ref_file", help="ANTIBODY_REFERENCE.csv")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()
data_dir = args$data_dir
ruv_sheet = args$ruv_sheet
validation_file = args$validation_file
mbc_md_file = args$mbc_md_file
ab_ref_file = args$ab_ref_file

#Things to include and exclue
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468.control","MDA468+EGF")
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
  
  PCsmd = cbind(PC1_and_PC2, repdf[match(PC1_and_PC2$type, repdf$sampcolumn), c("samp", "batch")])
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
sampname = gsub("_samplesheet.txt", "", basename(ruv_sheet))

#~ MATRIX OF ALL BATCHES TO COMBINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allbatch_norms = data.frame()
for (f in 1:length(unique(samps2ruv$Batch))){
  file_name = as.character(paste0(data_dir,"/",samps2ruv$Batch[f] ,"/3_IGG_SUBTRACTED.tsv"))
  #get normalized data
  batch = read.table(file = file_name, sep="\t",stringsAsFactors=F, row.names=1, header=T)
  batch_name = samps2ruv$Batch[f]
  
  #get only samples you want to ruv and controls
  sel_batch = batch[,c(grep(ctrlregex, colnames(batch)),
                       which(colnames(batch) %in% make.names(samps2ruv$Sample.Name)))]
  
  md = data.frame("batch"=batch_name, 
                  "sample" = colnames(sel_batch), 
                  "fullname" = paste0(batch_name,"__", colnames(sel_batch)))
  colnames(sel_batch) = paste0(batch_name,"__", colnames(sel_batch))
  
  if (f==1){
    allbatch_norms = data.matrix(sel_batch)
    md_allbatch_norms = md
  }
  else{
    allbatch_norms = data.matrix(cbind(allbatch_norms, sel_batch[match(rownames(allbatch_norms), rownames(sel_batch)),]))
    md_allbatch_norms = rbind(md_allbatch_norms, md)
  }
}
lallbatch_norms = log2(allbatch_norms +1)


## VALIDATION DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mbc_md = read.table(file= mbc_md_file, sep="\t", header=T)
igg_geosamp = read.csv(validation_file, sep= "\t", row.names = 1)
igg_geosamp = igg_geosamp[!grepl(omitregex, rownames(igg_geosamp)),]

ligg_geosamp=log2(igg_geosamp+1)
colnames(ligg_geosamp)=make.names(gsub("\\.1$", "", colnames(ligg_geosamp)))


#~ PREP DATASET FOR RUV~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ combining: rows combined with mbcctl, ctrls and samples
combined_allvalid = data.matrix(cbind(ligg_geosamp, lallbatch_norms[match(rownames(ligg_geosamp), rownames(lallbatch_norms)),]))
combined_md = data.frame(batch = sapply(strsplit(as.character(colnames(combined_allvalid)), "__"), `[`, 1), 
                         samp = gsub("\\.1$", "", sapply(strsplit(as.character(colnames(combined_allvalid)), "__"), `[`, 2)),
                         sampcolumn = c(colnames(combined_allvalid)),stringsAsFactors = F)
combined_md$mbc= ifelse(combined_md$batch %in% make.names(as.character(mbc_md$BatchID)), #MBC batch only
                        yes=ifelse(combined_md$samp %in% make.names(controls), yes="MBCcontrol", #controls MBC 
                                   no=ifelse(combined_md$samp %in% make.names(mbc_md$Sample), yes="MBCsample",no="notMBC")),#in mbc md sample
                        no=ifelse(combined_md$batch %in% unique(samps2ruv$Batch), yes="new", no="notMBC"))

#~ replicated matrix
#rep matrix based only on controls in MBC batches and newsamps
combined_md$reps=ifelse(grepl("MBCcontrol|new", combined_md$mbc), yes=combined_md$samp, no=combined_md$sampcolumn)
#~ replicated matrix 
repmat = replicate.matrix(combined_md$reps)

#~ PRE RUV RLEPLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~ box plot
lctls = combined_allvalid[,grep(ctrlregex, colnames(combined_allvalid))]
limits_rle = max(as.matrix(lctls))-median(as.matrix(lctls))
rle_orig = ruv_rle(Y = t(lctls), rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% colnames(lctls),]), 
                   ylim=c(-limits_rle,limits_rle)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  #geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label="PRE-RUV"), angle=90, hjust=0, size=2)+
  theme(legend.position = "right", legend.text = element_text(size=6)) + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression \nNormalized Data (IGG corrected) \n", sampname, "_", "PRE-RUV"))

pca_orig = pcaplt(lctls, title = paste0("log2 (Normalized Signal + 1)\n", sampname, "_", "PRE-RUV"), combined_md)



#~ RUV~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RUVcorrected = RUVIII(Y=t(log2(combined_allvalid +1)), ctl=c(1:nrow(combined_allvalid)), M=repmat, k=1)

#~ POST RUV RLE PLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ctl_ruv = t(RUVcorrected)[,grep(ctrlregex, colnames(t(RUVcorrected)))]
rle_ruv = ruv_rle(Y = RUVcorrected[grep(ctrlregex, rownames(RUVcorrected)),], 
                  rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% colnames(ctl_ruv),]), 
                  ylim=c(-limits_rle,limits_rle)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  #geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label="RUVcorrected"), angle=90, hjust=0, size=2) +
  theme(legend.position = "right") + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression RUV processed\n", sampname, "_", "RUVcorrected"))

pca_ruv = pcaplt(ctl_ruv, title = paste0("RUV Processed\n", sampname, "_", "RUVcorrected"), combined_md)

#~ RUV MELTING & ECDF~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ split these apart makes plotting easier

mbc_ruv = t(RUVcorrected)[,combined_md$sampcolumn[combined_md$mbc=="MBCsample"]]
samp_ruv = t(RUVcorrected)[,combined_md$sampcolumn[combined_md$mbc=="new"]]

#get samples only, exclude controls. Tricky because some samples have control names. dumb
idx_controls = which(!combined_md[combined_md$mbc=="new","reps"] %in% make.names(controls))
samp_ruv = samp_ruv[,idx_controls]
sampruv.m = melt(samp_ruv,  id.vars=row.names)


#~ Get percentile using ecdf
mbcruv.m = melt(mbc_ruv,  id.vars=row.names)
mbc_ecdf = tapply(mbcruv.m$value, mbcruv.m$Var1, ecdf)
mbc_percentile=data.frame(row.names = rownames(igg_geosamp))
mbcruv.m$ctl_probe = paste0(gsub("\\.1$", "", sapply(strsplit(as.character(mbcruv.m$Var2), split = "_"), tail,1)), 
                            "_", make.names(mbcruv.m$Var1))
mbcruv.m$cellline = factor(gsub("\\.1$", "", sapply(strsplit(as.character(mbcruv.m$Var2), split = "_"), tail,1)))



for (samp in colnames(samp_ruv)){
  samp_df = melt(samp_ruv[,samp, drop=FALSE])
  rownames(samp_df) = samp_df$Var1
  mbc_percentile[,samp] = as.vector(apply(samp_df, 1, function(x) 
    ((mbc_ecdf[[x["Var1"]]](x[["value"]]))))[rownames(mbc_percentile)])
}
#~ MAX AND MIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ruv = as.matrix(t(RUVcorrected))
mruv = melt(ruv)
ruv_stats = data.frame(
  #do.call(rbind, (tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, summary))),  #quartiles
  "min" = tapply(mruv$value, mruv$Var1, min),
  "q1" = tapply(mruv$value, mruv$Var1, function(x) quantile(x, 0.25)),
  "mean" = tapply(mruv$value, mruv$Var1, median),
  "q3" = tapply(mruv$value, mruv$Var1, function(x) quantile(x, 0.75)),
  "max" = tapply(mruv$value, mruv$Var1, max))


#~ ORDERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)
ab_order = ab_ref$X.AbID[order(ab_ref$Target)]
ab_order = ab_order[!grepl(omitregex, ab_order)]

mbcruv.m$Var1 = factor(mbcruv.m$Var1, levels=ab_order)
sampruv.m$Var1 = factor(sampruv.m$Var1, levels=ab_order)
ruv_stats$ab = factor(rownames(ruv_stats), levels=ab_order)


palette = c("#FF0000FF","#004CFFFF","#00A600FF","#984ea3","#ff7f00","#a65628")


#~ BOXPLOT MBC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bp_mbcruv = ggplot() +
  geom_linerange(
    data=ruv_stats, aes(x=ab, ymin = min, ymax = max),
    color = "#808080", 
    size = 7, 
    alpha = 0.7) +
  geom_boxplot(data = data.frame(mbcruv.m), mapping = aes(x = Var1, y=as.numeric(value)), outlier.colour = NA) +
  geom_point(data=sampruv.m, mapping=aes(x=Var1, y=value, colour=Var2),shape=8, size=2) +
  scale_color_manual(values=palette[1:length(unique(sampruv.m$Var2))]) +
  #scale_x_discrete(labels=paste0(as.character(levels(sampruv.m$Var1))," (",round(sampruv.m$value, 1), ",",round(mbc_percentile[as.character(levels(sampruv.m$Var1)), samp]*100,0),")")) +
  coord_flip() +
  labs(x="Antibody", title=paste0(sampname,  "\n within Distribution of Metastatic Breast Cancers"),y="RUVnormalized") +
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
pheatmap(mat = samp_ruv[,as.character(o_samps)],
                color             = colorRampPalette( c("green", "black", "red"), space="rgb")(10),
                cluster_cols = T,
                show_colnames     = TRUE,
                show_rownames     = TRUE,
                fontsize          = 9,
                #display_numbers = T,number_format = "%.2f", number_color="white",
                main              = sampname)
dev.off()
print("done with heatmaps making scatterplots")
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

plt_samps = samp_ruv
#colnames(plt_samps) = gsub("^.*1020_", "", colnames(plt_samps))
png(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_controls_RUV_scatter.png"), width = 850, height = 850)
plt_scatter = pairs(plt_samps, lower.panel = scatter_line, upper.panel = panel.cor, pch = 19)
dev.off()


#~ save all plots
ggsave(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_MBC.png"), bp_mbcruv, width = 8, height = 6)

ggsave(file=paste0(dirname(normalizePath(ruv_sheet)), "/", sampname, "_controls_norm_RLE.png"), 
       rle_orig, width = 9, height = 4.5)
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
write.table(x = round(mbc_percentile, 2), file=paste0(dirname(normalizePath(ruv_sheet)), "/", "combined_MBC_percentiles.tsv"), sep="\t", quote = F, row.names = T, col.names = NA)


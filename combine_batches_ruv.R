#!/usr/bin/env Rscript

##
## Usage
## ./gm_batch_correction.R /Users/patterja/Box\ Sync/NANOSTRING/data/ comb_sheet.txt validation_file mbc_md_file
## data_dir = where all your data lives
## comb_sheet = txt file samples and normalized matrices
## validation_file = txt file of validation data normalized
## mbc_md_file =txt with metastatic breast cancer list. 
##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(ruv))
suppressPackageStartupMessages(library(nanostring))


version="3.0"
mat_version = "20200320"


parser <- ArgumentParser()
parser$add_argument("--data_dir", type="character", default="/Volumes/OHSU/CLINICAL/Nanostring/output", 
                    dest="data_dir", help="directory name")
parser$add_argument("--comb_sheet", type="character", dest="comb_sheet", help="ruv samplesheet")
parser$add_argument("--validation_file", type="character",default=paste0("/Volumes/OHSU/CLINICAL/Nanostring/REFERENCE_FILES/validation_samples_rawdata_", mat_version, ".txt"),
                    dest="validation_file", help="validation file for controls comparison")
parser$add_argument("--md_file", type="character", default= "/Users/patterja/Box Sync/NANOSTRING/nanostring_metadata.xlsx",
                    dest="mbc_md_file", help="metastatic breast cancer metadata file")
parser$add_argument("--ab_ref_file", type="character", default= "/Volumes/OHSU/CLINICAL/Nanostring/REFERENCE_FILES/ANTIBODY_REFERENCE.csv",
                    dest="ab_ref_file", help="ANTIBODY_REFERENCE.csv")
parser$add_argument("--include_ctrls", action="store_true", default=FALSE,
                    dest="include_ctrls", help="include all antibodies")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()
data_dir = args$data_dir
comb_sheet = args$comb_sheet
validation_file = args$validation_file
md_file = args$mbc_md_file
ab_ref_file = args$ab_ref_file
include_ctrls = args$include_ctrls

#Things to include and exclue
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468.control","MDA468+EGF")
omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2")
omitregex = paste0(paste0("^", omit), collapse = "|")
ab.ctrl = "IgG|POS|NEG|^S6|^Histone"



#FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#todo: add these too nanostring package
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

#~ SAMPLESHEET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
outdir = dirname(normalizePath(comb_sheet))
print(paste0("output files will be in current dir: ", getwd()))

samps2batchcorr = read.table(file= comb_sheet, sep="\t", stringsAsFactors = F, header=T)
sampname = gsub("_samplesheet.txt", "", basename(comb_sheet))

# METADATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

md = read.xlsx(file=md_file, sheetName = "nansostring_metadata", check.names=T, stringsAsFactors=F)
md$sampcolumn = make.names(paste0(md$Batch, "__", md$Sample.Name))

# antibody metadata
ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)

#~ PARSE ALL SAMPLES FROM ALL BATCHES TO COMBINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
raw_dats = c()
raw_samps = c()
for (f in 1:length(unique(samps2batchcorr$Batch))){
  
  #get raw data
  batch_name = samps2batchcorr$Batch[f]
  file_name = as.character(paste0(data_dir,"/",batch_name ,"/rawdata.txt"))
  batch = read.table(file = file_name, sep="\t", row.names=2, stringsAsFactors=F, header=T, check.names = T)
  batch[,c("CodeClass", "Accession")] <- NULL
  
  #TEST to see if sample names in batch directories
  if(any(colnames(batch) %in% make.names(samps2batchcorr$Sample.Name))){
    print(paste0("Sample:", samps2batchcorr$Sample.Name[f], " in ", file_name))}
  else{
    print(paste0("Sample:", samps2batchcorr$Sample.Name[f], " NOT in ", file_name, "\n Check sample name."))
    stop()
  }
  #TEST to make sure the controls are in each batch
  for (ctrl in make.names(controls)){
    if (!ctrl %in% colnames(batch)){
      print(paste("WARNING: There is a control missing: ", ctrl))
    } else{
      print(paste0(ctrl, " controls accounted for."))
    }
  }
  
  #get only samples you want to batch and the controls
  idx_controls = which(make.names(colnames(batch)) %in% make.names(controls))
  idx_sample = which(make.names(colnames(batch)) %in% make.names(samps2batchcorr$Sample.Name))
  sel_batch = batch[,c(idx_controls, idx_sample)]
  sel_md = data.frame("batch"=batch_name,
                      "sample" = colnames(sel_batch),
                      "fullname" = make.names(paste0(batch_name,"__", colnames(sel_batch))), stringsAsFactors = F)
  
  colnames(sel_batch) = make.names(paste0(batch_name,"__", colnames(sel_batch)))
  
  # combine with all batches if not first one from samplesheet
  if (f==1){
    raw_dats = data.matrix(sel_batch)
    raw_samps = sel_md
  } else{
    raw_dats = cbind(raw_dats, sel_batch[match(rownames(raw_dats), rownames(sel_batch)),])
    raw_samps = rbind(raw_samps, sel_md)
  }
}
colnames(raw_dats) = paste0("combining", colnames(raw_dats))
raw_samps$fullname = paste0("combining", raw_samps$fullname)

# LOAD VALIDATION DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

validation = read.csv(validation_file, sep= "\t", row.names = 1, check.names = T)

# COMBINING DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## combine for replicate structure: rows combined with validation, ctrls and samples
combined_raw = cbind(validation, raw_dats[match(rownames(validation), rownames(raw_dats)),])
lcombined = log(combined_raw+1)

# COMBINED DATA METADATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## create a new metadata for replicate structure
combined_md = data.frame(batch = sapply(strsplit(as.character(colnames(lcombined)), "__"), `[`, 1), 
                         samp =sapply(strsplit(as.character(colnames(lcombined)), "__"), `[`, 2),
                         sampcolumn = c(colnames(lcombined)),stringsAsFactors = F)

combined_md$reps = ifelse(combined_md$samp %in% make.names(controls), yes=combined_md$samp, no=combined_md$sampcolumn)
combined_md$deid = ifelse(combined_md$samp %in% make.names(controls), yes=combined_md$samp, no="sample")
repmat = replicate.matrix(combined_md$reps)

# LIST OF SAMPLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## because i got sick of calling it all the time
combined_md_ordered = combined_md[order(combined_md$samp),]
celllines = combined_md_ordered$sampcolumn[combined_md_ordered$samp %in% make.names(controls)]
valid_celllines = setdiff(celllines, raw_samps$fullname)

valid_samps = intersect(md$sampcolumn[md$cohort=="validation"], combined_md$sampcolumn[!combined_md$samp %in% c(make.names(controls), raw_samps$sample)])
samps = raw_samps$fullname[!raw_samps$sample %in% make.names(controls)]
samp_celllines = setdiff(raw_samps$fullname, samps)


#~ RUV~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idx_ctl = grep(ab.ctrl, rownames(lcombined))

RUVcorrected = RUVIII(Y=t(lcombined), ctl=idx_ctl, M=repmat, k=2, include.intercept = FALSE, inputcheck = F)
RUVcorrected = RUVcorrected[,!grepl("POS|NEG", rownames(lcombined))]

bcdat = t(RUVcorrected)


#~ QC PLOTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~ RLE box plot
lctls=lcombined[!grepl("POS|NEG", rownames(lcombined)),celllines]

limits_rle_raw = max(as.matrix(lctls))-median(as.matrix(lctls))
rle_rawbox = ruv_rle(Y = t(lctls), rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% celllines,]), 
                  ylim=c(-limits_rle_raw,limits_rle_raw)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  geom_text(aes(x = rle.x.factor, y=-limits_rle_raw-0.5, label=colnames(lctls)), angle=90, hjust=0, size=2)+
  theme(legend.position = "right", legend.text = element_text(size=6)) + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression \nNormalized Data (IGG corrected) \n", sampname, "_", "PRE-RUV"))

ctl_bc = bcdat[,celllines]
limits_rle_bc = max(as.matrix(ctl_bc))-median(as.matrix(ctl_bc))

rle_ruvbox = ruv_rle(Y = t(ctl_bc), 
                   rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% celllines,]), 
                   ylim=c(-limits_rle_bc,limits_rle_bc)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  geom_text(aes(x = rle.x.factor, y=-limits_rle_bc-0.5, label="RUVcorrected"), angle=90, hjust=0, size=2) +
  theme(legend.position = "right") + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression RUV processed\n", sampname, "_", "RUVcorrected"))



# NORMAL HEATMAP
## rle raw
rownames(combined_md) = combined_md$sampcolumn
pheat_raw = pheatmap(mat = t(lctls),
                     color             = colorRampPalette(c("green", "black", "red"))(50),
                     border_color      = NA,
                     show_colnames     = TRUE,
                     show_rownames     = TRUE,
                     annotation_row    = combined_md[colnames(lctls),c("reps"), drop=F], 
                     drop_levels       = TRUE,
                     fontsize          = 8,
                     fontsize_row      = 5,
                     fontsize_col      = 5,
                     main              = "log raw", 
                     cluster_rows      = T,
                     cluster_cols      = T)

## ruv raw
pheat_ruv = pheatmap(mat = t(ctl_bc),
                     color             = colorRampPalette(c("green", "black", "red"))(50),
                     border_color      = NA,
                     show_colnames     = TRUE,
                     show_rownames     = TRUE,
                     annotation_row    = combined_md[colnames(lctls),c("reps"), drop=F], 
                     drop_levels       = TRUE,
                     fontsize          = 8,
                     fontsize_row      = 5,
                     fontsize_col      = 5,
                     main              = "RUV corrected controls", 
                     cluster_rows      = T,
                     cluster_cols      = T)

# RLE HEATMAP
## rle raw
rleraw = rel_log_exp(t(lctls))
pheat_rleraw = pheatmap(
  mat = rleraw,
  color             = colorRampPalette(c("green", "black", "red"))(50),
 border_color      = NA,
 show_colnames     = TRUE,
 show_rownames     = TRUE,
 annotation_row    = combined_md[colnames(lctls),c("reps"), drop=F], 
 drop_levels       = TRUE,
 fontsize          = 8,
 fontsize_row      = 5,
 fontsize_col      = 5,
 main              = "rle raw", 
 cluster_rows      = T,
 cluster_cols      = T)

## ruv rle
rleruv = rel_log_exp(t(ctl_bc))

pheat_rleruv = pheatmap(mat = (rleruv),
                     color             = colorRampPalette(c("green", "black", "red"))(50),
                     border_color      = NA,
                     show_colnames     = TRUE,
                     show_rownames     = TRUE,
                     annotation_row    = combined_md[colnames(lctls),c("reps"), drop=F], 
                     drop_levels       = TRUE,
                     fontsize          = 8,
                     fontsize_row      = 5,
                     fontsize_col      = 5,
                     main              = "RLE of RUV corrected controls", 
                     cluster_rows      = T,
                     cluster_cols      = T)
#~ PCA plots
pcacol= combined_md$deid[match(colnames(lctls), combined_md$sampcolumn)]
pca_raw = pcaplot(mat = lctls, title = paste0("log2 (Normalized Signal + 1)\n", sampname, "_", "PRE-RUV"), col=pcacol)

pcacol= combined_md$deid[match(colnames(ctl_bc), combined_md$sampcolumn)]
pca_ruv = pcaplot(ctl_bc, title = paste0("RUV Processed\n", sampname, "_", "RUVcorrected"), pcacol)

#~ TRA 
#t(RUVcorrected) is post norm, comb is pre norm, using no antibody filtered data
eruv = t(exp(ctl_bc))
rawctls = t(combined_raw)[rownames(eruv),colnames(eruv)]
tra = matrix(NA, nrow = 0, ncol = 4)

for (ctrl in make.names(controls)){
  #get matrix of celllines for controls not including batches of interest
  valid.ctrl_names = combined_md$sampcolumn[combined_md$samp == ctrl]

  selv.ctrl_raw = rawctls[valid.ctrl_names,, drop=F]
  selv.ctrl_ruv = eruv[valid.ctrl_names,, drop=F]
  tra_ctrlraw =c()
  tra_ctrlruv =c()
  for (r in 1:nrow(selv.ctrl_ruv)){
    tra_ctrlraw = rbind(tra_ctrlraw, log(sweep(as.matrix(selv.ctrl_raw[-r,,drop=F]), 2, as.numeric(selv.ctrl_raw[r,,drop=F]), `/`)))
    tra_ctrlruv = rbind(tra_ctrlruv, log(sweep(as.matrix(selv.ctrl_ruv[-r,,drop=F]), 2, as.numeric(selv.ctrl_ruv[r,,drop=F]), `/`)))
  }
  tra = rbind(tra, data.frame(melt(as.matrix(tra_ctrlraw)), "sample"="raw"))
  tra = rbind(tra, data.frame(melt(as.matrix(tra_ctrlruv)), "sample"="ruv"))
}

tra_btwn_batch = matrix(NA, nrow = 0, ncol = 4)
selraw = t(combined_raw[!grepl("POS|NEG", rownames(lcombined)),samps])
selruv = t(exp(bcdat[,samps]))
for (r in 1:(length(samps)-1)){
  print(r)
  selraw = log(sweep(as.matrix(selraw[-r,,drop=F]), 2, as.numeric(selraw[r,,drop=F]), `/`))
  selruv = log(sweep(as.matrix(selruv[-r,,drop=F]), 2, as.numeric(selruv[r,,drop=F]), `/`))
  tra_btwn_batch = rbind(tra_btwn_batch, data.frame(melt(as.matrix(selraw)), "sample"="raw"))
  tra_btwn_batch = rbind(tra_btwn_batch, data.frame(melt(as.matrix(selruv)), "sample"="ruv"))
}

ptra=ggplot(tra, aes(x=sample, y=value, fill=sample)) + 
  geom_boxplot() +
  facet_wrap(~Var2, scale="free") +
  geom_hline(yintercept =0, color="red") +
  labs(title=paste0("Distribution of TRA (technical replicate agreement) RAW and RUV Corrected\nTRAs calculated for new batch control versus each of the validation replicates\n"),
       y="log(validation count/new batch count)") +
  theme(legend.position = "bottom")

ptra_btwn=ggplot(tra_btwn_batch, aes(x=sample, y=value, fill=sample)) + 
  geom_boxplot() +
  geom_hline(yintercept =0, color="red") +
  labs(title=paste0("Distribution of TRA between RAW and RUV Corrected between samplesheet batches\n"),
       y="log(validation count/new batch count)") +
  theme(legend.position = "bottom")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create("qc_figures", showWarnings = F)


print("saving batch correction plots")
pdf(file=paste0("qc_figures/",sampname,"_TRA.png"), width = 8, height = 8)
print(plot(ptra))
print(plot(ptra_btwn))
dev.off()

pdf(paste0("qc_figures/",sampname,"_RLE_boxplots.pdf"), width = 8, height = 4)
print(rle_rawbox)
print(rle_ruvbox)
dev.off()

pdf(paste0("qc_figures/",sampname,"_signal_heatmap.pdf"), width = 7, height = 8)
grid::grid.newpage()
grid::grid.draw(pheat_raw)
grid::grid.newpage()
grid::grid.draw(pheat_ruv)
dev.off()

pdf(paste0("qc_figures/",sampname,"_RLE_heatmap.pdf"), width = 7, height = 8)
grid::grid.newpage()
grid::grid.draw(pheat_rleraw)
grid::grid.newpage()
grid::grid.draw(pheat_rleruv)
dev.off() 

pdf(paste0("qc_figures/",sampname,"_PCA.pdf"), width = 8, height = 8)
print(pca_raw)
print(pca_ruv)
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~ VALIDATION DATA FILTERING and AB ORDER FOR PLOTTING~~~~~~~~~~~~~~~~~~~~~~~~~~
if (include_ctrls){
  print("keeping all antibodies: IgG and antibodies that did not perform well or did not have dynamic range")
  fbcdat = bcdat
  other_abs=setdiff(colnames(bcdat),ab_ref$X.AbID)
  ab_order = c(ab_ref$X.AbID[order(ab_ref$Target)], other_abs)
  
} else {
  fbcdat = bcdat[!grepl(omitregex, rownames(bcdat)),]
  #AB_ORDER
  ab_order = ab_ref$X.AbID[order(ab_ref$Target)]
  ab_order = ab_order[!grepl(omitregex, ab_order)]
  
}
#~ SPLIT AND MELT FOR PLOTTING #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

validsamp_dat = fbcdat[,valid_samps]
samps_dat = fbcdat[,samps]
validsamp.datm = melt(validsamp_dat,  id.vars=row.names)
samps.datm = data.frame(melt(samps_dat, id.vars=row.names), "detectable"=TRUE)


# ANTIBODY THRESHOLDING setting with detectable flag ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## if signal for antibody below sample igg then turned to minimum of validation  cohort
## if signal above do nothing
ab_file="/Volumes/OHSU/CLINICAL/Nanostring/REFERENCE_FILES/ANTIBODY_REFERENCE.csv"  
ab_ref = read.csv(ab_file, sep=",", stringsAsFactors=F)


for (samp in samps){
  newsamp = combined_raw[,samp, drop=F]
  
  rbigg = newsamp[which(rownames(newsamp)=="RbAb-IgG"),]
  mmigg = newsamp[which(rownames(newsamp)=="MmAb-IgG1"),]
  for (i in seq(1:length(ab_ref$X.AbID))){
    ab = (ab_ref$X.AbID)[i]
    if (ab_ref$Host[which(ab_ref$X.AbID==ab)]=="rabbit"){
      val = newsamp[ab,]-rbigg
      print("rabbit")
    } else if (ab_ref$Host[which(ab_ref$X.AbID==ab)]=="mouse"){
      val = newsamp[ab,]-mmigg
      print("mouse")
    } else {
      print("double check antibody name")
      stop()
    }
    if (val < 0){
      #get index of samples with ab and samp
      idx = which(samps.datm$Var1==ab & samps.datm$Var2==samp)
      samps.datm[idx,"detectable"] = FALSE
    }
  }
}

#~ ECDF~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~ Get percentile using ecdf
valid_ecdf = tapply(validsamp.datm$value, validsamp.datm$Var1, ecdf)

samp_percentile=data.frame(row.names = rownames(samps_dat))
validsamp.datm$batch = sapply(strsplit(as.character(validsamp.datm$Var2), split = "__"), head,1)
validsamp.datm$samp = sapply(strsplit(as.character(validsamp.datm$Var2), split = "__"), tail,1)


for (samp in colnames(samps_dat)){
  samp_df = melt(samps_dat[,samp, drop=FALSE])
  rownames(samp_df) = samp_df$Var1
  samp_percentile[,samp] = as.vector(apply(samp_df, 1, function(x) 
    ((valid_ecdf[[x["Var1"]]](x[["value"]]))))[rownames(samp_percentile)])
}



#~ ORDERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
validsamp.datm$Var1 = factor(validsamp.datm$Var1, levels=ab_order)
samps.datm$Var1 = factor(samps.datm$Var1, levels=ab_order)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURES FOR SAMPLE COMPARISON
#~ BOXPLOT MBC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
palette = c("#FF0000FF","#004CFFFF","#00A600FF","#984ea3","#ff7f00","#a65628")

minval = min(as.numeric(validsamp.datm[, "value"]), na.rm=T)-1
samps.datm$newvalue = ifelse(samps.datm$detectable, samps.datm$value, minval)


bp = ggplot(validsamp.datm, aes(x=Var1, y=value)) +
  geom_boxplot() 
bp = bp + geom_text(data=samps.datm[samps.datm$detectable==FALSE,], mapping=aes(x=Var1, y=newvalue, colour=Var2), label="ND", size=3,position=position_jitter(width=c( 0.1)))
bp = bp + geom_point(data=samps.datm[samps.datm$detectable==TRUE,], mapping=aes(x=Var1, y=newvalue, colour=Var2), shape=8, size=2)

bp = bp + labs(x="Antibody", title=sampname) +
  scale_color_manual(values=c("red", "blue")) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major=element_line(colour="gray"),
        plot.title = element_text(hjust = 0.5, vjust=0),
        legend.text=element_text(size=8),
        legend.position="bottom",
        axis.text.x = element_text(size=9, colour=c("black"), angle = -90, hjust=0),
        axis.text.y = element_text(size=9, colour="black")) +
  coord_flip()

#~ HEATMAP #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pheat_comparison = pheatmap(mat = samps_dat[,samps],
         color             = colorRampPalette( c("green", "black", "red"), space="rgb")(10),
         cluster_cols = T,
         show_colnames     = TRUE,
         show_rownames     = TRUE,
         fontsize          = 9,
         main              = sampname)


png(file=paste0(sampname, "_comparison_heatmap.png"), width = 350, height = 800)
print(pheat_comparison)
dev.off()


#~ SCATTER PAIRS 
#colnames(plt_samps) = gsub("^.*1020_", "", colnames(plt_samps))
pdf(file=paste0(sampname, "_RUV_scatter.png"), width = 8, height = 6)
plt_scatter = pairs(samps_dat, lower.panel = scatter_line, upper.panel = panel.cor, pch = 19, main=paste0("Correlation of ", sampname))
dev.off()

ggsave(file=paste0(sampname, "_boxplots.png"), bp, width = 8, height = 6)

print("done with comparison plots saving tables")

#~ TABLES #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(x = t(bcdat), file=paste0(sampname, "_RUVcorrected.csv"), 
            sep=",", quote = F, row.names = T, col.names = NA)
write.table(x = t(combined_raw), file=paste0(sampname, "_raw.csv"), 
            sep=",", quote = F, row.names = T, col.names = NA)
write.table(x = round(samp_percentile, 2), file=paste0(sampname, "_combined_samp_percentiles.tsv"), sep="\t", quote = F, row.names = T, col.names = NA)


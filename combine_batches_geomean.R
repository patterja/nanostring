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

version="3.0"
mat_version = "20200320"


parser <- ArgumentParser()
parser$add_argument("--data_dir", type="character", default="/Volumes/Histopathology Shared Resource/CLINICAL/Nanostring/output", 
                    dest="data_dir", help="directory name")
parser$add_argument("--comb_sheet", type="character", dest="comb_sheet", help="ruv samplesheet")
parser$add_argument("--validation_file", type="character",default=paste0("/Volumes/Histopathology Shared Resource/CLINICAL/Nanostring/REFERENCE_FILES/validation_samples_rawdata_", mat_version, ".txt"),
                    dest="validation_file", help="validation file for controls comparison")
parser$add_argument("--md_file", type="character", default= "/Users/patterja/Box Sync/NANOSTRING/nanostring_metadata.xlsx",
                    dest="mbc_md_file", help="metastatic breast cancer metadata file")
parser$add_argument("--ab_ref_file", type="character", default= "/Volumes/Histopathology Shared Resource/CLINICAL/Nanostring/REFERENCE_FILES/ANTIBODY_REFERENCE.csv",
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
ab.ctrl = "IgG|NEG|^S6|^Histone"

#FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gm_mean = function(x){
  #' @param x: vector
  #' @export
  g=exp(mean(log(x+1, base = exp(1))))
  g=g-1
  return(g)
}

pcaplt <- function (mat, title = "PCA Plot", col=rownames(mat)) {
  #' pca
  #'
  #' @param mat (matrix/dataframe):  mat
  #' @param title (character) : 
  #' @param subtype (dataframe metadata) : rownames correspnond to
  #' @param labe (character) : rownames correspnond to
  #' @return ggplot 
  col = c(col)
  var = mat[apply(mat, 1, var, na.rm = TRUE) != 0, ]
  cc.var = var[complete.cases(var), ]
  pca_prcomp = prcomp(t(var), center = T, scale = F)
  PC1_and_PC2 = data.frame(PC1 = pca_prcomp$x[, 1], PC2 = pca_prcomp$x[,2], type = rownames(pca_prcomp$x))
  perc = (pca_prcomp$sdev^2)/sum(pca_prcomp$sdev^2) * 100
  labs <- sapply(seq_along(perc), function(i) {
    paste("PC ", i, " (", round(perc[i], 2), "%)", sep = "")})
  
  PCsmd = cbind(PC1_and_PC2, col=col)
  levs = levels(factor(col))
  cols =c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#8DD3C7","#FFFFB3",
          "#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F",
          "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#8DD3C7","#FFFFB3",
          "#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")
  p = ggplot(PCsmd,aes_string("PC1", "PC2", col="col")) + 
    geom_point(size = 1.5) + 
    geom_text(aes(label = PCsmd$type), vjust = -1, size=2) + 
    labs(title = title,x = labs[1], y = labs[2]) + 
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "gray90"), 
          panel.border = element_rect(colour = "gray90", fill=NA),
          plot.title = element_text(hjust = 0.5), 
          legend.text = element_text(size = 4), legend.position = "right") +
    scale_colour_manual(values =cols[1:length(levs)]) +
    xlim(-20,20) +ylim(-9,9)
  return(p)
}

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
## OUTPUT PREP ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir.create(paste0(dirname(normalizePath(comb_sheet)),"/batchcorrection_figures"), showWarnings = F)

#~ SAMPLESHEET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samps2batchcorr = read.csv(file= comb_sheet, sep="\t", stringsAsFactors = F)
sampname = gsub("_samplesheet.txt", "", basename(comb_sheet))

# METADATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

md = read.xlsx(file=md_file, sheetName = "nansostring_metadata", check.names=T, stringsAsFactors=F)
md$sampcolumn = make.names(paste0(md$Batch, "__", md$Sample.Name))

# antibody metadata
ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)

#~ PARSE ALL SAMPLES FROM ALL BATCHES TO COMBINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allbatch_norms = c()
md_allbatch_norms = c()
for (f in 1:length(unique(samps2batchcorr$Batch))){
  
  #get raw data
  batch_name = samps2batchcorr$Batch[f]
  file_name = as.character(paste0(data_dir,"/",batch_name ,"/rawdata.txt"))
  batch = read.table(file = file_name, sep="\t", row.names=2, stringsAsFactors=F, header=T, check.names = T)
  batch[,c("CodeClass", "Accession")] <- NULL
  
  #TEST to see if sample names in batch
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
  
  #get only samples you want to batch and controls
  idx_controls = which(make.names(colnames(batch)) %in% make.names(controls))
  idx_sample = which(make.names(colnames(batch)) %in% make.names(samps2batchcorr$Sample.Name))
  sel_batch = batch[,c(idx_controls, idx_sample)]
  sel_md = data.frame("batch"=batch_name,
                      "sample" = colnames(sel_batch),
                      "fullname" = make.names(paste0(batch_name,"__", colnames(sel_batch))), stringsAsFactors = F)

  colnames(sel_batch) = make.names(paste0(batch_name,"__", colnames(sel_batch)))
  
  # combine with all batches if not first one from samplesheet
  if (f==1){
    allbatch_norms = data.matrix(sel_batch)
    md_allbatch_norms = sel_md
  } else{
    allbatch_norms = cbind(allbatch_norms, sel_batch[match(rownames(allbatch_norms), rownames(sel_batch)),])
    md_allbatch_norms = rbind(md_allbatch_norms, sel_md)
  }
}

# VALIDATION DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

validation = read.csv(validation_file, sep= "\t", row.names = 1, check.names = T)
validation = validation[,md$sampcolumn[md$cohort=="validation"]]
## remove samples from validation cohort if they exists. 
validation = validation[,!colnames(validation) %in% make.names(colnames(allbatch_norms))]

# COMBINING DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## combine for metadata to repmat: rows combined with validation, ctrls and samples
combined_allvalid = cbind(validation, allbatch_norms[match(rownames(validation), rownames(allbatch_norms)),])
combined_allvalid = combined_allvalid[!grepl("NEG|POS", rownames(combined_allvalid)),]
lcombined = log(combined_allvalid+1)

# COMBINING METADATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## create a new metadata
combined_md = data.frame(batch = sapply(strsplit(as.character(colnames(combined_allvalid)), "__"), `[`, 1), 
                         samp =sapply(strsplit(as.character(colnames(combined_allvalid)), "__"), `[`, 2),
                         sampcolumn = c(colnames(combined_allvalid)),stringsAsFactors = F)

# replicated matrix
combined_md$reps = ifelse(combined_md$samp %in% make.names(controls), yes=combined_md$samp, no=combined_md$sampcolumn)
combined_md$deid = ifelse(combined_md$samp %in% make.names(controls), yes=combined_md$samp, no="sample")
ordered_combmd = combined_md[order(combined_md$samp),]


# LIST OF SAMPLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## because i got sick of calling it all the time
celllines = combined_md$sampcolumn[combined_md$samp %in% make.names(controls)]
valid_celllines = setdiff(celllines, md_allbatch_norms$fullname)
valid_samps = combined_md$sampcolumn[!combined_md$samp %in% c(make.names(controls), md_allbatch_norms$sample)]
samps = md_allbatch_norms$fullname[!md_allbatch_norms$sample %in% make.names(controls)]
samp_celllines = setdiff(md_allbatch_norms$fullname, samps)


#~ BATCH CORRECTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCALE BY GEOMEAN
idx_controls = which(combined_md$samp %in% controls)
comb.noercc = combined_allvalid[!grepl("NEG|POS", rownames(combined_allvalid)),]
comb.noercc.ctrls = combined_allvalid[!grepl("NEG|POS", rownames(combined_allvalid)),idx_controls]

gm_cf = gm_mean(apply(comb.noercc.ctrls, 2, gm_mean))/apply(comb.noercc, 2, gm_mean)
comb.norm = t(t(comb.noercc)* gm_cf)

#LOG
bcdat = log2(comb.norm+1)

#~ PRE RUV RLEPLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~ box plot
lctls=lcombined[!grepl("NEG", rownames(lcombined)),celllines]

limits_rle = max(as.matrix(lctls))-median(as.matrix(lctls))
rle_pre = ruv_rle(Y = t(lctls), rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% celllines,]), 
                  ylim=c(-limits_rle,limits_rle)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label=colnames(lctls)), angle=90, hjust=0, size=2)+
  theme(legend.position = "right", legend.text = element_text(size=6)) + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression \nNormalized Data (IGG corrected) \n", sampname, "_", "PRE-RUV"))
#~ pca plot
pcacol= combined_md$deid[match(colnames(lctls), combined_md$sampcolumn)]
pca_pre = pcaplt(mat = lctls, title = paste0("log2 (Normalized Signal + 1)\n", sampname, "_", "PRE-RUV"), col=pcacol)

#~ POST BC RLE PLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ctl_bc = bcdat[,celllines]
rle_post = ruv_rle(Y = t(ctl_bc), 
                  rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% celllines,]), 
                  ylim=c(-limits_rle,limits_rle)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label="RUVcorrected"), angle=90, hjust=0, size=2) +
  theme(legend.position = "right") + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression RUV processed\n", sampname, "_", "RUVcorrected"))
#~ pca plot
pcacol= combined_md$deid[match(colnames(ctl_bc), combined_md$sampcolumn)]
pca_post = pcaplt(ctl_bc, title = paste0("RUV Processed\n", sampname, "_", "RUVcorrected"), pcacol)

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
samps.datm = melt(samps_dat, id.vars=row.names)


#~ TRA & PLOTTING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
edat = (exp(fbcdat))-1
prenorm = comb.noercc[ab_order,]

for (ctrl in make.names(controls)){
  #get matrix of celllines for controls not including batches of interest
  v.ctrl_names = setdiff(combined_md$sampcolumn[combined_md$samp == ctrl], samp_celllines)
  s.ctrl_names = intersect(combined_md$sampcolumn[combined_md$samp == ctrl], samp_celllines)
  
  selv.ctrl_pre = prenorm[,v.ctrl_names, drop=F]
  sels.ctrl_pre = prenorm[,s.ctrl_names, drop=F]
  selv.ctrl_post = data.frame(edat[,v.ctrl_names, drop=F])
  sels.ctrl_post = data.frame(edat[,s.ctrl_names, drop=F])
  
  tra_ctrlpre = do.call(rbind,apply(selv.ctrl_pre, 2, function(x) log(x/sels.ctrl_pre)))
  tra_ctrlpost = do.call(rbind,apply(selv.ctrl_post, 2, function(x) log(x/sels.ctrl_post)))
  
  tra = rbind(data.frame(melt(as.matrix(tra_ctrlpre)), "sample"="0pre"),data.frame(melt(as.matrix(tra_ctrlpost)), "sample"="1post"))
  tra$ab = sapply(strsplit(as.character(tra$Var1), ".", fixed=T), tail, 1)
  tra$samp_norm_label = paste0(tra$Var2, "_", tra$sample)
  
  p=ggplot(tra, aes(x=samp_norm_label, y=value, fill=samp_norm_label)) + 
    geom_violin() +
    facet_wrap(~ab, scale="free") +
    geom_hline(yintercept =0, color="red") +
    labs(title=paste0("Distribution of TRA (technical replicate agreement) Pre and Post Normalization\nCompared To Each Controls In Validation Batch\n", sampname, "_",ctrl),
         y="log(validation count/sample count)") +
    theme(legend.position = "bottom",
          legend.direction ="vertical",
          axis.text.x = element_blank()) +
    guides(fill=guide_legend(ncol=2))
  png(file=paste0(dirname(normalizePath(comb_sheet)), "/batchcorrection_figures/",ctrl, "_TRA.png"), width = 800, height = 800)
  plot(p)
  dev.off()
  write.table(tra, file=paste0(dirname(normalizePath(comb_sheet)),"/batchcorrection_figures/",ctrl, "_TRA_table.txt"),sep="\t", quote = F, row.names = T, col.names = NA)
  #TRA of pre/post
  tra_samp = matrix(NA, nrow = length(ab_order), ncol=length(s.ctrl_names), dimnames = list(ab_order, s.ctrl_names))
  for (n in s.ctrl_names){
    tra_samp[,n] = log(sels.ctrl_pre[ab_order,n]/sels.ctrl_post[ab_order,n])
  }
  write.table(tra_samp, file=paste0(dirname(normalizePath(comb_sheet)),"/batchcorrection_figures/",ctrl, "_TRA_btwn_batch.txt"),sep="\t", quote = F, row.names = T, col.names = NA)
  
}


# ANTIBODY THRESHOLDING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# if signal for antibody below sample igg then turned to minimum of validation  cohort
# if signal above do nothing

for (samp in samps){
  rbigg = allbatch_norms[which(rownames(allbatch_norms)=="RbAb-IgG"),samp]
  mmigg = allbatch_norms[which(rownames(allbatch_norms)=="MmAb-IgG1"),samp]
  for (i in seq(1:length(ab_order))){
    ab = ab_order[i]
    if (ab_ref$Host[which(ab_ref$X.AbID==ab)]=="rabbit"){
      val = allbatch_norms[ab,samp]-mmigg
    } else if (ab_ref$Host[which(ab_ref$X.AbID==ab)]=="mouse"){
      val = allbatch_norms[ab,samp]-rbigg
    } else {
      print("double check antibody name")
      stop()
    }
    if (val < 0){
      print(paste0("replacing ", samp, " for ", ab))
      print(paste0("from:", val, ", to:", min(fbcdat[ab,])))
      val = min(fbcdat[ab,])
      fbcdat[ab,samp] = val
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

# VALID STAT SUMMARY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fbcdatm = melt(as.matrix(fbcdat))

bcstats = data.frame(
  #do.call(rbind, (tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, summary))),  #quartiles
  "min" = tapply(fbcdatm$value, fbcdatm$Var1, min),
  "q1" = tapply(fbcdatm$value, fbcdatm$Var1, function(x) quantile(x, 0.25)),
  "median" = tapply(fbcdatm$value, fbcdatm$Var1, median),
  "q3" = tapply(fbcdatm$value, fbcdatm$Var1, function(x) quantile(x, 0.75)),
  "max" = tapply(fbcdatm$value, fbcdatm$Var1, max))

#~ ORDERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
validsamp.datm$Var1 = factor(validsamp.datm$Var1, levels=ab_order)
samps.datm$Var1 = factor(samps.datm$Var1, levels=ab_order)
bcstats$ab = factor(rownames(bcstats), levels=ab_order)


#~ BOXPLOT MBC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
palette = c("#FF0000FF","#004CFFFF","#00A600FF","#984ea3","#ff7f00","#a65628")

bp = ggplot(bcstats, aes(x=ab, y=q1)) +
  #geom_boxplot(aes(factor(Var1), as.numeric(value)), outlier.colour = NA) +
  geom_crossbar(aes(ymin = min, ymax = q1), width = 0.9, color="#606060") +
  geom_crossbar(aes(ymin = q1, ymax = q3), width = 0.9, color="#606060", fatten=0.5, fill="#808080") +
  geom_crossbar(aes(ymin = q3, ymax = max), width = 0.9, color="#606060", fatten=0.5) +
  
  #geom_linerange(data=bcstats, aes(x=ab, ymin = min, ymax = max),color = "#808080",size = 7, alpha = 0.7) +
  geom_point(data=samps.datm, mapping=aes(x=Var1, y=value, colour=Var2),shape=8, size=2) +
  scale_color_manual(values=palette[1:length(unique(samps.datm$Var2))]) +
  #scale_x_discrete(labels=paste0(as.character(levels(samps.datm$Var1))," (",round(samps.datm$value, 1), ",",round(mbc_percentile[as.character(levels(samps.datm$Var1)), samp]*100,0),")")) +
  coord_flip() +
  labs(x="Antibody", title=paste0(sampname,  "\n within Distribution of Metastatic Breast Cancers"),y="gmnormalized") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major=element_line(colour="gray"),
        plot.title = element_text(hjust = 0.5, vjust=0),
        legend.text=element_text(size=8),
        legend.position="bottom",
        axis.text.x = element_text(hjust = 1, size=7, colour="black"),
        axis.text.y = element_text(size=6, colour="black")) +
  guides(col = guide_legend(ncol = 1))


#~ HEATMAP #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

png(file=paste0(dirname(normalizePath(comb_sheet)), "/", sampname, "_heatmap.png"), width = 350, height = 800)
pheatmap(mat = samps_dat[,samps],
                color             = colorRampPalette( c("green", "black", "red"), space="rgb")(10),
                cluster_cols = T,
                show_colnames     = TRUE,
                show_rownames     = TRUE,
                fontsize          = 9,
                #display_numbers = T,number_format = "%.2f", number_color="white",
                main              = sampname)
dev.off()
print("done with heatmaps making scatterplots")

#~ SCATTER PAIRS #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#colnames(plt_samps) = gsub("^.*1020_", "", colnames(plt_samps))
png(file=paste0(dirname(normalizePath(comb_sheet)), "/", sampname, "_controls_GM_scatter.png"), width = 850, height = 850)
plt_scatter = pairs(samps_dat, lower.panel = scatter_line, upper.panel = panel.cor, pch = 19)
dev.off()


#~ save all plots
ggsave(file=paste0(dirname(normalizePath(comb_sheet)), "/", sampname, "_boxrange.png"), bp, width = 8, height = 6)

ggsave(file=paste0(dirname(normalizePath(comb_sheet)), "/", sampname, "_controls_raw_RLE.png"), 
       rle_pre, width = 9, height = 4.5)
ggsave(file=paste0(dirname(normalizePath(comb_sheet)), "/", sampname, "_controls_raw_PCA.png"), 
       pca_pre, width = 8, height = 7)
ggsave(file=paste0(dirname(normalizePath(comb_sheet)), "/", sampname, "_controls_GM_RLE.png"), 
       rle_post, width = 9, height = 4.5)
ggsave(file=paste0(dirname(normalizePath(comb_sheet)), "/", sampname, "_controls_GM_PCA.png"), 
       pca_post, width = 8, height = 7)
ggsave(file=paste0(dirname(normalizePath(comb_sheet)), "/", sampname, "_controls_GM_scatter.png"), 
       plt_scatter, width = 8, height = 7)

write.table(x = t(bcdat), file=paste0(dirname(normalizePath(comb_sheet)), "/", sampname, "_batchcorrected.csv"), 
            sep=",", quote = F, row.names = T, col.names = NA)
write.table(x = round(samp_percentile, 2), file=paste0(dirname(normalizePath(comb_sheet)), "/", "combined_MBC_percentiles.tsv"), sep="\t", quote = F, row.names = T, col.names = NA)


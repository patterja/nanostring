#!/usr/bin/env Rscript

## This function takes the rawdata and runs RUV using the rawdata validation matrix. 
## QC metrics, RUV and comparison with MBC cohort 
## Usage
## ./ruv_batchcorrection.R -i rawdata.txt --validation_file validation_samples_rawdata.txt --md_file nanostring_metadata.xlsx --ab_ref_file ANTIBODY_REFERENCE.csv
##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pheatmap))

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ruv))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(nanostring))


## ARGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
version="3.0"
mat_version = "20200320"

parser <- ArgumentParser()

parser$add_argument("-i", help="rawdata.txt", dest="input_file")
parser$add_argument("--validation_file", type="character", 
                    default=paste0("/Volumes/OHSU/CLINICAL/Nanostring/REFERENCE_FILES/validation_samples_rawdata_", mat_version, ".txt"),
                    dest="validation_file", help="validation file for controls comparison")
parser$add_argument("--md_file", type="character", default= "/Users/patterja/Box Sync/NANOSTRING/nanostring_metadata.xlsx",
                    dest="md_file", help="metadata file")
parser$add_argument("--ab_ref_file", type="character", default= "/Volumes/OHSU/CLINICAL/Nanostring/REFERENCE_FILES/ANTIBODY_REFERENCE.csv",
                    dest="ab_ref_file", help="ANTIBODY_REFERENCE.csv")
parser$add_argument("--include_ctrls", action="store_true", default=FALSE,
                    dest="include_ctrls", help="include all antibodies")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()
includeBCCL = args$includeBCCL

input_file = args$input_file
validation_file = args$validation_file
md_file = args$md_file
ab_ref_file = args$ab_ref_file
include_ctrls = args$include_ctrls

## TEST ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!grepl("raw", input_file)){
  print("Error: The input file is not rawdata.txt")
  stop()
}

## OUTPUT PREP ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir.create("ruv_figures", showWarnings = F)

## INPUT DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# validation file
validation = read.csv(validation_file, sep= "\t", row.names = 1, check.names = T)

# antibody metadata
ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)

# NEW BATCH
new_batch = read.table(file = input_file, sep="\t", row.names=2, stringsAsFactors=F, header=T, check.names = T)
new_batch[,c("CodeClass", "Accession")] <- NULL

# metadata
md = read.xlsx(file=md_file, sheetName = "nanostring_metadata", check.names=T, stringsAsFactors=F)
md$sampcolumn = make.names(paste0(md$Batch, "__", md$Sample.Name))

#Things to include and exclude
controls=make.names(c("MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468+EGF"))
ab.ctrl = make.names(rownames(new_batch)[grepl("IgG|POS|NEG|^S6|^Histone", rownames(new_batch))])
omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2", "NEG", "POS")
omitregex = paste0(paste0("^", omit), collapse = "|")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOOP THRU DATASET #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#empty data frame
samp_percentile=data.frame(row.names = rownames(validation)[!grepl(omitregex, rownames(validation))])
for (samp in setdiff(colnames(new_batch), controls)) {
  print(samp)
  
  # COMBINING: rows combined with validation  cohort, columns=ctrls and 1 samp only
  idx_controls = which(colnames(new_batch) %in% controls)
  newsampctl = new_batch[,c(idx_controls,which(colnames(new_batch)==samp))]
  colnames(newsampctl) = paste0("newbatch", "__", colnames(newsampctl))
  comb = cbind(validation, newsampctl[match(rownames(validation), rownames(newsampctl)),])
  # METADATA: adjust metadata to match
  comb_md = data.frame(batch = make.names(sapply(strsplit(as.character(colnames(comb)), "__"), `[`, 1)), 
                       samp =make.names(sapply(strsplit(as.character(colnames(comb)), "__"), `[`, 2)),
                       sampcolumn = c(colnames(comb)),stringsAsFactors = F, check.names = T)
  md = md[md$sampcolumn %in% comb_md$sampcolumn,]
  valid_controls=c(md$sampcolumn[md$cohort=="validation" & md$Study=="control"],paste0("newbatch__", controls))
  # REPLICATE STRUCTURE: repmat based only on controls in validation cohort for RUV input
  comb_md$reps = ifelse(comb_md$sampcolumn %in% valid_controls, yes=comb_md$samp, no=comb_md$sampcolumn)
  sort_comb_md = comb_md[order(comb_md$samp),]
  
  #LOG
  lcomb=log(comb+1)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~ BATCH CORRECTION: RUV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  repmat = replicate.matrix(comb_md$reps)
  negctl = which(make.names(rownames(lcomb)) %in% ab.ctrl)
  # r by c matrix, r=observations, c=features
  RUVcorrected = RUVIII(Y=t(lcomb), ctl=negctl, M=repmat, k=2, include.intercept = FALSE)
  RUVcorrected = RUVcorrected[,!grepl("NEG|POS", colnames(RUVcorrected))]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~ PLOTTING QC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LOG RAW   ## filter out negatives and ab_order
  lctls = lcomb[!grepl("NEG|POS", rownames(lcomb)),] 
  lctls = lctls[!make.names(rownames(lctls)) %in% ab.ctrl, sort_comb_md$sampcolumn[sort_comb_md$samp %in% controls]]
  limits_rle_raw = max(as.matrix(lctls))-median(as.matrix(lctls))
  # RUV 
  ctl_ruv = t(RUVcorrected)[, sort_comb_md$sampcolumn[sort_comb_md$samp %in% controls]]
  limits_rle_ruv = max(as.matrix(ctl_ruv))-median(as.matrix(ctl_ruv))
  
  # RLE BOXPLOTS
  ## rle raw
  rle_raw = ruv_rle(Y = t(lctls), 
                     rowinfo = as.matrix(sort_comb_md[sort_comb_md$samp %in% controls,]), 
                     ylim=c(-limits_rle_raw,limits_rle_raw)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle_raw-0.5, label=samp), angle=90, hjust=0, size=2)+
    theme(legend.position = "right", legend.text = element_text(size=6)) + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression Raw Data",  samp))
  
  ## rle batch corrected
  rle_ruv = ruv_rle(Y = t(ctl_ruv[!make.names(rownames(lctls)) %in% ab.ctrl,]), 
                    rowinfo = as.matrix(sort_comb_md[sort_comb_md$samp %in% controls,]), 
                    ylim=c(-limits_rle_ruv,limits_rle_ruv)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle_ruv-0.5, label=samp), angle=90, hjust=0, size=2) +
    theme(legend.position = "right") + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression Batch Corrected\n", samp))
  
  
  # NORMAL HEATMAP
  rownames(comb_md) = comb_md$sampcolumn
  pheat_raw = pheatmap(mat = t(lctls),
    color             = colorRampPalette(c("green", "black", "red"))(50),
    border_color      = NA,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    annotation_row    = comb_md[colnames(lctls),c("reps"), drop=F], 
    drop_levels       = TRUE,
    fontsize          = 8,
    fontsize_row      = 5,
    fontsize_col      = 5,
    main              = "log raw", 
    cluster_rows      = T,
    cluster_cols      = T)
  dev.off()
  
  ## NORMAL ruv
  pheat_ruv = pheatmap(mat = t(ctl_ruv),
                       color             = colorRampPalette(c("green", "black", "red"))(50),
                       border_color      = NA,
                       show_colnames     = TRUE,
                       show_rownames     = TRUE,
                       annotation_row    = comb_md[colnames(lctls),c("reps"), drop=F], 
                       drop_levels       = TRUE,
                       fontsize          = 8,
                       fontsize_row      = 5,
                       fontsize_col      = 5,
                       main              = "log raw", 
                       cluster_rows      = T,
                       cluster_cols      = T)
  # RLE HEATMAP
  rownames(comb_md) = comb_md$sampcolumn
  pheat_raw = pheatmap(mat = t(lctls),
                       color             = colorRampPalette(c("green", "black", "red"))(50),
                       border_color      = NA,
                       show_colnames     = TRUE,
                       show_rownames     = TRUE,
                       annotation_row    = comb_md[colnames(lctls),c("reps"), drop=F], 
                       drop_levels       = TRUE,
                       fontsize          = 8,
                       fontsize_row      = 5,
                       fontsize_col      = 5,
                       main              = "log raw", 
                       cluster_rows      = T,
                       cluster_cols      = T)
  dev.off()
  
  ## rle ruv
  pheat_ruv = pheatmap(mat = t(ctl_ruv),
                       color             = colorRampPalette(c("green", "black", "red"))(50),
                       border_color      = NA,
                       show_colnames     = TRUE,
                       show_rownames     = TRUE,
                       annotation_row    = comb_md[colnames(lctls),c("reps"), drop=F], 
                       drop_levels       = TRUE,
                       fontsize          = 8,
                       fontsize_row      = 5,
                       fontsize_col      = 5,
                       main              = "log raw", 
                       cluster_rows      = T,
                       cluster_cols      = T)
  #~ pca
  pca_raw = pcaplot(mat = lctls, 
                    title="Normalized, Pre-RUV", 
                    col = sort_comb_md$samp[sort_comb_md$samp %in% controls])
  pca_ruv = pcaplot(mat = ctl_ruv, 
                   title="RUV processed", 
                   col = sort_comb_md$samp[sort_comb_md$samp %in% controls])
 
  #~ TRA 
  #t(RUVcorrected) is post norm, comb is pre norm, using no antibody filtered data
  eruv = t(exp(ctl_ruv))
  rawctls = t(comb[colnames(eruv),rownames(eruv)])
  tra = matrix(NA, nrow = 0, ncol = 4)
  
  for (ctrl in make.names(controls)){
    #get matrix of celllines for controls not including batches of interest
    valid.ctrl_names = comb_md$sampcolumn[comb_md$samp == ctrl & !comb_md$batch=="newbatch"]
    samp.ctrl_names = comb_md$sampcolumn[comb_md$samp == ctrl & comb_md$batch=="newbatch"]
    
    selv.ctrl_raw = rawctls[valid.ctrl_names,, drop=F]
    sels.ctrl_raw = rawctls[samp.ctrl_names,, drop=F]
    selv.ctrl_ruv = (eruv[valid.ctrl_names,, drop=F])
    sels.ctrl_ruv = (eruv[samp.ctrl_names,, drop=F])
    
    tra_ctrl_raw = log(sweep(as.matrix(selv.ctrl_raw), 2, as.numeric(sels.ctrl_raw), `/`))
    tra_ctrl_ruv = log(sweep(as.matrix(selv.ctrl_ruv), 2, as.numeric(sels.ctrl_ruv), `/`))
    
    
    tra = rbind(tra, data.frame(melt(as.matrix(tra_ctrl_raw)), "sample"="raw"))
    tra = rbind(tra, data.frame(melt(as.matrix(tra_ctrl_ruv)), "sample"="ruv"))
    }
    
  ptra=ggplot(tra, aes(x=sample, y=value, fill=sample)) + 
    geom_boxplot() +
    facet_wrap(~Var2, scale="free") +
    geom_hline(yintercept =0, color="red") +
    labs(title=paste0("Distribution of TRA (technical replicate agreement) RAW and RUV Corrected\nTRAs calculated for new batch control versus each of the validation replicates\n", samp, "_",ctrl),
         y="log(validation count/new batch count)") +
    theme(legend.position = "bottom")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print("saving batch correction plots")
  png(file=paste0("ruv_figures/",samp,"_", ctrl, "_TRA.png"), width = 800, height = 800)
  plot(ptra)
  dev.off()
  
  pdf(paste0("ruv_figures/",samp,"_RLE_boxplots.pdf"), width = 8, height = 4)
  print(rle_raw)
  print(rle_ruv)
  dev.off()
  
  pdf(paste0("ruv_figures/",samp,"_RLE_heatmap.pdf"), width = 7, height = 8)
  print(pheat_raw)
  print(pheat_ruv)
  dev.off()
  
  pdf(paste0("ruv_figures/",samp,"_PCA.pdf"), width = 8, height = 8)
  print(pca_raw)
  print(pca_ruv)
  dev.off()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #~ VALIDATION DATA PREP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (include_ctrls){
      print("keeping all antibodies: IgG and antibodies that did not perform well or did not have dynamic range")
      nmat = t(RUVcorrected)
      other_abs=setdiff(colnames(RUVcorrected),ab_ref$X.AbID)
      ab_order = c(ab_ref$X.AbID[order(ab_ref$Target)], other_abs)
    } else {
      nmat = t(RUVcorrected)
      nmat = nmat[!grepl(omitregex, rownames(nmat)),]
      #AB_ORDER
      ab_order = ab_ref$X.AbID[order(ab_ref$Target)]
      ab_order = ab_order[!grepl(omitregex, ab_order)]
    }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  #antibody threshold ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # if signal for antibody below sample igg then turned to minimum of validation  cohort
  # if signal above do nothing
  newsamp = comb[,paste0("newbatch", "__", (samp)),drop=F]
  rbigg = newsamp[which(rownames(newsamp)=="RbAb-IgG"),]
  mmigg = newsamp[which(rownames(newsamp)=="MmAb-IgG1"),]
  for (i in seq(1:length(ab_order))){
    ab = ab_order[i]
    if (ab_ref$Host[which(ab_ref$X.AbID==ab)]=="rabbit"){
      val = newsamp[ab,]-rbigg
    } else if (ab_ref$Host[which(ab_ref$X.AbID==ab)]=="mouse"){
      val = newsamp[ab,]-mmigg
    } else {
      print("double check antibody name")
      stop()
    }
    if (val < 0){
      val = min(nmat[ab,])
      nmat[ab,paste0("newbatch", "__", (samp))] = val
      print(val)
    }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    #~ split these apart makes plotting easier
    #validation only and newsamp only
    val_nmat = nmat[,as.character(comb_md$sampcolumn[md$cohort=="validation" & !md$Study=="control"] ),drop=F]
    new_nmat = nmat[,as.character(comb_md$sampcolumn[comb_md$batch=="newbatch" & !comb_md$samp %in% controls]),drop=F]
    #~ MELTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    val.m = melt(val_nmat,  id.vars=row.names)
    new.m = melt(new_nmat, id.vars=row.names)
    valid_ecdf = tapply(val.m$value, val.m$Var1, ecdf)
    mat_new.m = as.matrix(new.m)
    
    for (idx in seq(1:nrow(mat_new.m))){
      new_samp=mat_new.m[idx,]
      samp_percentile[new_samp[1],samp] = valid_ecdf[[new_samp[1]]](new_samp[3])
    }
    
    m.mat = melt(as.matrix(nmat))
    
    norm_stats = data.frame(
      "min" = tapply(m.mat$value, m.mat$Var1, min),
      "q1" = tapply(m.mat$value, m.mat$Var1, function(x) quantile(x, 0.25)),
      "median" = tapply(m.mat$value, m.mat$Var1, median),
      "q3" = tapply(m.mat$value, m.mat$Var1, function(x) quantile(x, 0.75)),
      "max" = tapply(m.mat$value, m.mat$Var1, max))
    
    val.m$Var1 = factor(val.m$Var1, levels=ab_order)
    new.m$Var1 = factor(new.m$Var1, levels=ab_order)
    norm_stats$ab = factor(rownames(norm_stats), levels=ab_order)
    
    ab_labels= paste0(as.character(levels(new.m$Var1))," (",round(new.m$value, 1), ",",round(samp_percentile[as.character(levels(new.m$Var1)), samp]*100,0),")")
    
    #PLOTTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #~ boxplot of MBC
    bp = ggplot(norm_stats, aes(x=ab, y=q1)) +
      geom_crossbar(aes(ymin = min, ymax = q1), width = 0.9, color="#606060",fill="#ececec") +
      geom_crossbar(aes(ymin = q1, ymax = q3), width = 0.9, color="#606060", fatten=0.5, fill="#c0c0c0") +
      geom_crossbar(aes(ymin = q3, ymax = max), width = 0.9, color="#606060", fatten=0.5,fill="#808080") +
      geom_point(data=new.m, mapping=aes(x=factor(Var1), y=value), colour=c('red'), shape=8, size=2) +
      #geom_linerange(
      #  data=norm_stats, aes(x=ab, ymin = min, ymax = max),
      #  color = "#808080", 
      #  size = 7, 
      #  alpha = 0.7) +
      #geom_boxplot(aes(factor(Var1), as.numeric(value)), outlier.colour = NA) +
      #scale_x_discrete(labels=paste0(as.character(levels(new.m$Var1))," (",
      #                               round(samp_percentile[as.character(levels(new.m$Var1)), samp]*100,0),")")) +
      labs(x="Antibody", title=paste0(samp,  "\n within Distribution of Metastatic Breast Cancers"), y="RUV corrected") +
      theme(panel.background = element_rect(fill = "white"),
            panel.grid.major=element_line(colour="gray"),
            plot.title = element_text(hjust = 0.5, vjust=0),
            legend.text=element_text(size=8),
            legend.position="right",
            axis.text.x = element_text(size=7, colour=c("black")),
            axis.text.y = element_text(size=7, colour="black")) +
      coord_flip()
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~ save all plots
    print("saving boxplot plots and tables")
    
    ggsave(file=paste0(samp, "_boxplot.pdf"), bp, device="pdf", width = 8, height = 6)
    write.table(lcomb[!grepl("NEG|POS", colnames(lcomb)),], file=paste0("ruv_figures/",samp,"_lograw.txt"),sep="\t", quote = F, row.names = T, col.names = NA)
    write.table(t(RUVcorrected), file=paste0("ruv_figures/",samp,"_ruvcorrected.txt"),sep="\t", quote = F, row.names = T, col.names = NA)
  
  }
  
  tar(tarfile=paste0("ruv_figures.tar.gz"), files=paste0("ruv_figures"), compression="gzip", tar="tar")
  write.table(x = round(samp_percentile, 2), file=paste0("samp_percentiles.tsv"), sep="\t", quote = F, row.names = T, col.names = NA)
  
  
  
  
  
  
  
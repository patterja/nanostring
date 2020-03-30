#!/usr/bin/env Rscript

## This function takes the rawdata and runs RUV using the rawdata validation matrix. 
## QC metrics, RUV and comparison with MBC cohort 
## Usage
## ./ruv_batchcorrection.R -i rawdata.txt --validation_file validation_samples_rawdata.txt --md_file nanostring_metadata.xlsx --ab_ref_file ANTIBODY_REFERENCE.csv
##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ruv))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(gridExtra))

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


## FUNCTIONS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
md = read.xlsx(file=md_file, sheetName = "nansostring_metadata", check.names=T, stringsAsFactors=F)
md$sampcolumn = make.names(paste0(md$Batch, "__", md$Sample.Name))

#Things to include and exclude
controls=make.names(c("MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468+EGF"))
ab.ctrl = make.names(rownames(new_batch)[grepl("IgG|NEG|^S6|^Histone", rownames(new_batch))])
omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2", "NEG", "POS")
omitregex = paste0(paste0("^", omit), collapse = "|")

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
  comb = comb[!grepl("POS", rownames(comb)),]
  # METADATA: adjust metadata to match
  comb_md = data.frame(batch = make.names(sapply(strsplit(as.character(colnames(comb)), "__"), `[`, 1)), 
                       samp =make.names(sapply(strsplit(as.character(colnames(comb)), "__"), `[`, 2)),
                       sampcolumn = c(colnames(comb)),stringsAsFactors = F, check.names = T)
  md = md[md$sampcolumn %in% comb_md$sampcolumn,]
  valid_controls=c(md$sampcolumn[md$cohort=="validation" & md$Study=="control"],paste0("newbatch__", controls))
  # REPLICATE MATRIX: rep matrix based only on controls in validation cohort
  comb_md$reps = ifelse(comb_md$sampcolumn %in% valid_controls, yes=comb_md$samp, no=comb_md$sampcolumn)
  sort_comb_md = comb_md[order(comb_md$samp),]
  
  #LOG
  lcomb=log2(comb+1)
  #~ BATCH CORRECTION: RUV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  repmat = replicate.matrix(comb_md$reps)
  negctl = which(make.names(rownames(lcomb)) %in% ab.ctrl)
  # r by c matrix, r=observations, c=features
  RUVcorrected = RUVIII(Y=t(lcomb), ctl=negctl, M=repmat, k=1, eta=1)
  RUVcorrected = RUVcorrected[,!grepl("NEG", colnames(RUVcorrected))]
  
  #~ RLE PLOT ~ pre-ruv figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## filter out negatives and ab_order
  lctls = lcomb[!grepl("NEG|POS", rownames(lcomb)),] 
  lctls = lctls[!make.names(rownames(lctls)) %in% ab.ctrl, sort_comb_md$sampcolumn[sort_comb_md$samp %in% controls]]
  limits_rle = max(as.matrix(lctls))-median(as.matrix(lctls))
  
  #rle 
  rle_orig = ruv_rle(Y = t(lctls), 
                     rowinfo = as.matrix(sort_comb_md[sort_comb_md$samp %in% controls,]), 
                     ylim=c(-limits_rle,limits_rle)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label=samp), angle=90, hjust=0, size=2)+
    theme(legend.position = "right", legend.text = element_text(size=6)) + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression Raw Data",  samp))
  #pca
  pca_orig = pcaplt(mat = lctls, 
                    title="Normalized, Pre-RUV", 
                    col = sort_comb_md$samp[sort_comb_md$samp %in% controls])
  
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  #~ PLOTTING RUV processed figures
  ctl_ruv = t(RUVcorrected)[, sort_comb_md$sampcolumn[sort_comb_md$samp %in% controls]]
  limits_rle = max(as.matrix(ctl_ruv))-median(as.matrix(ctl_ruv))
  
  rle_ruv = ruv_rle(Y = t(ctl_ruv[!make.names(rownames(lctls)) %in% ab.ctrl,]), 
                    rowinfo = as.matrix(sort_comb_md[sort_comb_md$samp %in% controls,]), 
                    ylim=c(-limits_rle,limits_rle)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label=samp), angle=90, hjust=0, size=2) +
    theme(legend.position = "right") + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression RUV processed\n", samp))
  
  pca_ruv = pcaplt(mat = ctl_ruv, 
                   title="RUV processed", 
                   col = sort_comb_md$samp[sort_comb_md$samp %in% controls])
  #~ RLE TABLE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  prerle = t(apply(lctls, 1, function(x) x-median(x)))
  postrle = t(apply(ctl_ruv, 1, function(x) x-median(x)))
  
  #~ VALIDATION DATA PREP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  #~ TRA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #t(RUVcorrected) is post norm, comb is pre norm, using no antibody filtered data
  edat = (exp(t(RUVcorrected)))
  prenorm = comb[rownames(edat),]
  for (ctrl in make.names(controls)){
    tra = matrix(NA, nrow = 0, ncol = 3)
    png(file=paste0("ruv_figures/",samp,"_", ctrl, "_TRA.png"), width = 800, height = 800)
    #get matrix of celllines for controls not including batches of interest
    v.ctrl_names = comb_md$sampcolumn[comb_md$samp == ctrl & !comb_md$batch=="newbatch"]
    s.ctrl_names = comb_md$sampcolumn[comb_md$samp == ctrl & comb_md$batch=="newbatch"]
    
    selv.ctrl_pre = prenorm[,v.ctrl_names, drop=F]
    sels.ctrl_pre = prenorm[,s.ctrl_names, drop=F]
    selv.ctrl_post = data.frame(edat[,v.ctrl_names, drop=F])
    sels.ctrl_post = data.frame(edat[,s.ctrl_names, drop=F])
    
    tra_ctrlpre = do.call(cbind, apply(selv.ctrl_pre, 2, function(x) log(x/sels.ctrl_pre)))
    tra = data.frame(melt(as.matrix(tra_ctrlpre)), "sample"="pre")
    tra_ctrlpost = do.call(cbind, apply(selv.ctrl_post, 2, function(x) log(x/sels.ctrl_post)))
    tra = rbind(tra, data.frame(melt(as.matrix(tra_ctrlpost)), "sample"="post"))
    
    
    p=ggplot(tra, aes(x=sample, y=value, fill=sample)) + 
      geom_violin() +
      facet_wrap(~Var1, scale="free") +
      geom_hline(yintercept =0, color="red") +
      labs(title=paste0("Distribution of TRA (technical replicate agreement) Pre and Post Normalization\nCompared To Each Controls In Validation Batch\n", samp, "_",ctrl),
           y="log(validation count/sample count)") +
      theme(legend.position = "bottom")
    
    plot(p)
    dev.off()
    tratab_pre = setNames(tra_ctrlpre, colnames(selv.ctrl_pre))
    tratab_post = setNames(tra_ctrlpost, colnames(selv.ctrl_post))
    write.table(tratab_pre, file=paste0("ruv_figures/",samp,"_", ctrl, "_TRA_table.txt"),sep="\t", quote = F, row.names = T, col.names = NA)
    write.table(tratab_post, file=paste0("ruv_figures/",samp,"_", ctrl,  "_TRA_table.txt"),sep="\t", quote = F, row.names = T, col.names = NA)
  }
  
  #antibody threshold ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # if signal for antibody below sample igg then turned to minimum of validation  cohort
  # if signal above do nothing
  newsamp = comb[,paste0("newbatch", "__", (samp)),drop=F]
  rbigg = newsamp[which(rownames(newsamp)=="RbAb-IgG"),]
  mmigg = newsamp[which(rownames(newsamp)=="MmAb-IgG1"),]
  for (i in seq(1:length(ab_order))){
    ab = ab_order[i]
    if (ab_ref$Host[which(ab_ref$X.AbID==ab)]=="rabbit"){
      val = newsamp[ab,]-mmigg
    } else if (ab_ref$Host[which(ab_ref$X.AbID==ab)]=="mouse"){
      val = newsamp[ab,]-rbigg
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
    #  mbc_percentile[,samp] =  as.vector(apply(newruv.m, 1, function(x) 
    #    ((mbc_ecdf[[as.character(x["Var1"])]](x[["value"]]))))[rownames(mbc_percentile)])
    #max and min
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
      geom_crossbar(aes(ymin = min, ymax = q1), width = 0.9, color="#606060") +
      geom_crossbar(aes(ymin = q1, ymax = q3), width = 0.9, color="#606060", fatten=0.5, fill="#808080") +
      geom_crossbar(aes(ymin = q3, ymax = max), width = 0.9, color="#606060", fatten=0.5) +
      geom_point(data=new.m, mapping=aes(x=factor(Var1), y=value), colour=c('red'), shape=8, size=2) +
      #geom_linerange(
      #  data=norm_stats, aes(x=ab, ymin = min, ymax = max),
      #  color = "#808080", 
      #  size = 7, 
      #  alpha = 0.7) +
      #geom_boxplot(aes(factor(Var1), as.numeric(value)), outlier.colour = NA) +
      #scale_x_discrete(labels=paste0(as.character(levels(new.m$Var1))," (",
      #                               round(samp_percentile[as.character(levels(new.m$Var1)), samp]*100,0),")")) +
      labs(x="Antibody (Percentiles)", title=paste0(samp,  "\n within Distribution of Metastatic Breast Cancers"), y="Geomean Normalized") +
      theme(panel.background = element_rect(fill = "white"),
            panel.grid.major=element_line(colour="gray"),
            plot.title = element_text(hjust = 0.5, vjust=0),
            legend.text=element_text(size=8),
            legend.position="right",
            axis.text.x = element_text(size=7, colour=c("black")),
            axis.text.y = element_text(size=7, colour="black")) +
      coord_flip()
    print("saving RUV plots")
    
    print("saving normalization plots")
    pdf(paste0("ruv_figures/",samp,"_hierarchical_clustering_plot.pdf"), width = 7, height = 7)
    par(cex=0.7)
    plot(hclust(dist(t(lctls[!make.names(rownames(lctls)) %in% ab.ctrl,]), method ="euclidean")),
                 main="Raw Data")
    plot(hclust(dist(t(ctl_ruv), method ="euclidean")), main="RUV Clustering")
    
    dev.off()
    
    
    #~ save all plots
    ggsave(file=paste0(samp, "_MBC.pdf"), bp, device="pdf", width = 8, height = 6)
    write.table(prerle, file=paste0("ruv_figures/",samp,"_", ctrl, "_raw_RLEtable.txt"),sep="\t", quote = F, row.names = T, col.names = NA)
    write.table(postrle, file=paste0("ruv_figures",samp,"_", ctrl,  "_bc_RLEtable.txt"),sep="\t", quote = F, row.names = T, col.names = NA)
    ggsave(file=paste0("ruv_figures/",samp, "_controls_norm_RLE.pdf"),device="pdf", rle_orig, width = 9, height = 4.5)
    ggsave(file=paste0("ruv_figures/",samp, "_controls_norm_PCA.pdf"),device="pdf", pca_orig, width = 8, height = 7)
    ggsave(file=paste0("ruv_figures/",samp, "_controls_RUV_RLE.pdf"),device="pdf", rle_ruv, width = 9, height = 4.5)
    ggsave(file=paste0("ruv_figures/",samp, "_controls_RUV_PCA.pdf"),device="pdf", pca_ruv, width = 8, height = 7)
    
    write.table(x = nmat, file=paste0("ruv_figures/",samp, "_RUVcorrected.tsv"), sep="\t", quote = F, row.names = T, col.names = NA)
  }
  
  tar(tarfile=paste0("ruv_figures.tar.gz"), files=paste0("ruv_figures"), compression="gzip", tar="tar")
  write.table(x = round(samp_percentile, 2), file=paste0("samp_percentiles.tsv"), sep="\t", quote = F, row.names = T, col.names = NA)
  
  
  
  
  
  
  
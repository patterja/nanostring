#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)

usage = "
./norm_geomean.R input_file <normal/all> 

      input_file= txt file after processing RCC
      normal= controls are MCF7,HCC1954,BT474,HeyA8,MDA468 control, MDA468+EGF
      all= all column in the matrix are used as controls
"

argsLen <- length(args);
if (argsLen == 2) {
  input_file = args[1]
  ctrl_flag = args[2]
  if (!ctrl_flag %in% c("normal", "all")) {
    cat("second parameter must be normal or all")
    stop(cat(usage))
  }
  print(paste0("Processing ", input_file))
} else {
  stop(cat(usage))
}
  
#Functions and Libraries

library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(openxlsx)

#Functions
gm_mean = function(x){
  #' @param x: vector
  #' @export
  g=exp(mean(log(x)))
  return(g)
}

gm_sd = function(x){
  #' standard deviation of row
  #' @param x: vector
  #' @param gm: geomean
  #' @export
  gm=exp(mean(log(x)))
  g = exp(sqrt(sum(log(x/gm)^2)/(length(x))))
  #g = exp(sd(log(x)))
  return(g)
}


get_output <- function(input_file){
  #' Geomean normalized
  #' @param input_file: nanonstring input file name
  #' @output 
  #' @return 
  #' @export
  dat = read.csv(input_file, sep="\t", check.names = F, stringsAsFactors = F)
  metadat = dat[,c(1:3)]
  #colnames(dat) = c(colnames(dat)[1:ncol(metadat)], as.character(dat[1,(ncol(metadat)+1):ncol(dat)]))
  #dat = dat[-1,]
  samples = apply((dat[!grepl("NEG|POS", dat[,"Name"]),c((ncol(metadat)+1):ncol(dat))]), 2, as.numeric)
  rownames(samples) = dat$Name[!grepl("NEG|POS", dat[,"Name"])]
  
  pos= apply(as.matrix(dat[grepl("POS", dat[,"Name"]),c((ncol(metadat)+1):ncol(dat))]), 2, as.numeric)
  rownames(pos) = dat[grepl("POS", dat[,"Name"]),"Name"]
  
  neg= apply(as.matrix(dat[grepl("NEG", dat[,"Name"]),c((ncol(metadat)+1):ncol(dat))]), 2, as.numeric)
  rownames(neg) = dat[grepl("NEG", dat[,"Name"]),"Name"]
  
  return(list("pos"=pos, "neg"=neg, "samples"=samples))
  
}


geomean_norm <- function(samples, controls, pos, mm, rb){
  #' Geomean normalized
  #' @param samples: matrix. First column is gene names, rows are genes, columns are samples.
  #' @param controls: list of column names that are the control
  #' @param pos: matrix. of positive controls
  #' @param mm: list of mouse antibodies
  #' @param rb: list of rabbit antibodies
  #' @return pcaplot
  #' @export
  
  #normalize lize by positive control
  pos_cf=mean(apply(pos[,(colnames(pos) %in% controls)], 2, sum))/apply(pos[,controls], 2, sum)
  mat_poscf = t(t(samples)*pos_cf)
  
  #normalize by column/Sample
  geomean_cfcol =mean(apply(mat_poscf[,(colnames(mat_poscf) %in% controls)], 2, gm_mean))/apply(mat_poscf[,controls], 2, gm_mean)
  mat_gmcf = t(t(mat_poscf)*geomean_cfcol)
  
  #get mouse and rabbit antibody background
  threshold_MmAb = (mat_gmcf[grepl("MmAb-IgG1", rownames(mat_gmcf)),])
  threshold_RbAb = (mat_gmcf[grepl("RbAb-IgG", rownames(mat_gmcf)),])
  
  #proportional increase from IGG=1
  mmnorm = t(t(mat_gmcf[which(rownames(mat_gmcf) %in% mm),])/(threshold_MmAb))
  rbnorm = t(t(mat_gmcf[which(rownames(mat_gmcf) %in% rb),])/(threshold_RbAb))
  igg_corrected = rbind(mmnorm, rbnorm)
  
  #get mean and sd and se
  geomean_cfrow = (apply(igg_corrected[,controls], 1, gm_mean))
  geomsd_cfrow = (apply(igg_corrected[,controls], 1, gm_sd))
  se_logAb = (apply(igg_corrected[,controls], 1, function(x) 1.96*(sd(log(x))/sqrt(length(x)))))
  
  #zscore table
  #zmat = t(t(log(igg_corrected)-log(geomean_cfrow))/log(geomsd_cfrow))
  zmat = t(t(log(igg_corrected)-log(geomean_cfrow)))
  
  #upper confidence interval
  up_iggc = exp(log(igg_corrected) + se_logAb)
  
  #lower confidence interval
  down_iggc = exp(log(igg_corrected) - se_logAb)
  zup = t(t(log(up_iggc)-log(geomean_cfrow)))
  #zup = t(t(log(up_iggc)-log(geomean_cfrow))/log(geomsd_cfrow))
  zdown = t(t(log(down_iggc)-log(geomean_cfrow)))
  #zdown = t(t(log(down_iggc)-log(geomean_cfrow))/log(geomsd_cfrow))
  
  return(list("ercc_normalized" = mat_poscf[,c(controls, setdiff(colnames(samples), controls))], 
              "geomean_normbysample"=mat_gmcf[,c(controls, setdiff(colnames(samples), controls))], 
              "igg_geosamp_corrected"= igg_corrected[,c(controls, setdiff(colnames(samples), controls))], 
              "igg_geosamp_corrected_up" = up_iggc[,c(controls, setdiff(colnames(samples), controls))], 
              "igg_geosamp_corrected_down" = down_iggc[,c(controls, setdiff(colnames(samples), controls))], 
              "residual"=zmat[,c(controls, setdiff(colnames(samples), controls))], 
              "residual_up" = zup[,c(controls, setdiff(colnames(samples), controls))], 
              "residual_down" = zdown[,c(controls, setdiff(colnames(samples), controls))], 
              "stats"=data.frame("geomean_ab" = geomean_cfrow, "se_ab" = se_logAb)))
}


#process input file

raw_data = get_output(input_file)

#set controls
if (ctrl_flag == "normal"){
  controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468+EGF")
} else {
  controls = colnames(raw_data$samples)
}

mouseAb =c("Ki-67", "Pan-Keratin", "S6 Ribosomal", "p53", "MmAb-IgG1")
mm = c(rownames(raw_data$samples)[grepl(paste0(paste0("^", mouseAb), collapse="|"),  rownames(raw_data$samples))])
rb = c(setdiff(rownames(raw_data$samples), mm))

norm_dat = geomean_norm(samples = raw_data$samples, controls = controls, pos = raw_data$pos, mm, rb)

outfile_name = gsub(".txt", "", input_file)

write.xlsx(norm_dat, file =paste0(outfile_name, "_NORMALIZED.xlsx"))


mzmat = melt(norm_dat$residual)
mzmat$key = paste0(mzmat$Var1, mzmat$Var2)
mzup = melt(norm_dat$residual_up)
mzup$key = paste0(mzup$Var1, mzup$Var2)
mzdown = melt(norm_dat$residual_down)
mzdown$key = paste0(mzdown$Var1, mzdown$Var2)


mz_ci = data.frame(mzmat[,c("Var1", "Var2")], 
                   lower = mzdown[match(mzmat$key, mzdown$key),"value"], 
                   upper = mzup[match(mzmat$key, mzup$key),"value"], 
                   value = mzmat$value)

mz_ci$Var2 = factor(mz_ci$Var2, levels=unique(sort(as.character(mz_ci$Var2))))
mz_ci$status = ifelse(mz_ci$lower < 0 & mz_ci$upper < 0 & mz_ci$value < 0, "low", 
                      ifelse(mz_ci$lower > 0 & mz_ci$upper > 0 & mz_ci$value > 0, "high", "ambiguous"))
mz_ci$id = paste0(round(mz_ci$lower, 1),",",round(mz_ci$value, 1),",",round(mz_ci$upper, 1))

#for coloring because this color scheme is lame but requested by chris
z_rescale = c(min(mz_ci$value), quantile(mz_ci$value, probs = 0.25), 0, quantile(mz_ci$value, probs = 0.25), max(mz_ci$value))
z_color = (z_rescale - range(z_rescale)[1]) / diff(range(z_rescale)) * diff(c(0,1)) + c(0,1)[1]
pmzmat = ggplot(mz_ci, aes(Var2, Var1)) +
  geom_tile(aes(fill = value), colour="white") +
  geom_text(aes(label = round(value, 2)), size=3, colour="white") +
  #scale_fill_gradient2(low = "green", mid = "black", high = "red", midpoint = 0, guide = "legend")
  scale_fill_gradientn(colours=(colorRampPalette( c("green", "black", "red"))(5)),
                       values= z_color, limits=c(min(mz_ci$value),max(mz_ci$value))) +
  labs(y="Antibody", y="Sample", 
       title=paste0(gsub(".txt", "", tail(unlist(strsplit(split="/", x=input_file)), n=1)), ":\nResidual from Mean of Ab")) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major=element_line(colour="gray"),
        plot.title = element_text(hjust = 0.5, vjust=0),
        legend.text=element_text(size=8),
        legend.position="right",
        axis.text.x = element_text(angle=90,hjust = 1, size=10, colour="black"),
        axis.text.y = element_text(size=10, colour="black"))


pm_ci= ggplot(mz_ci, aes(Var2, Var1)) +
  geom_tile(aes(fill = value), colour="white") +
  geom_text(aes(label = id), size=2, colour="white") +
  scale_fill_gradientn(colours=(colorRampPalette( c("green", "black", "red"))(5)),
                       values= z_color, limits=c(min(mz_ci$value),max(mz_ci$value))) +
  labs(y="Antibody", y="Sample",
       title=paste0(gsub(".txt", "", tail(unlist(strsplit(split="/", x=input_file)), n=1)), ":\n Residual of Confidence Intervals")) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major=element_line(colour="gray"),
        plot.title = element_text(hjust = 0.5, vjust=0),
        legend.text=element_text(size=8),
        legend.position="bottom",
        axis.text.x = element_text(angle=90,hjust = 1, size=8, colour="black"),
        axis.text.y = element_text(size=8, colour="black"))

#migg
migg = melt(norm_dat$igg_geosamp_corrected)
migg$key = paste0(migg$Var1, migg$Var2)
miggup = melt(norm_dat$igg_geosamp_corrected_up)
miggup$key = paste0(miggup$Var1, miggup$Var2)
miggdown = melt(norm_dat$igg_geosamp_corrected_down)
miggdown$key = paste0(miggdown$Var1, miggdown$Var2)

migg_ci = data.frame(migg[,c("Var1", "Var2")], 
                   lower = miggdown[match(migg$key, miggdown$key),"value"], 
                   upper = miggup[match(migg$key, miggup$key),"value"], 
                   value = migg$value)
migg_ci$id = paste0(round(migg_ci$lower, 1),",",round(migg_ci$value, 1),",",round(migg_ci$upper, 1))
migg_ci$Var2 = factor(migg_ci$Var2, levels=unique(sort(as.character(migg_ci$Var2))))

pm_migg= ggplot(migg_ci, aes(Var2, Var1)) +
  geom_tile(aes(fill = value), colour="black") +
  geom_text(aes(label = id), size=1.7, colour="black") +
  scale_fill_gradientn(colours=(colorRampPalette( c("green", "black", "red"))(10))) +
  labs(y="Antibody", y="Sample",
       title=paste0(gsub(".txt", "", tail(unlist(strsplit(split="/", x=input_file)), n=1)), ":\n Normalized Confidence Interval")) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major=element_line(colour="gray"),
        plot.title = element_text(hjust = 0.5, vjust=0),
        legend.text=element_text(size=8),
        legend.position="bottom",
        axis.text.x = element_text(angle=90,hjust = 1, size=8, colour="black"),
        axis.text.y = element_text(size=8, colour="black"))

ggsave(filename = paste0(outfile_name, "_residual_HEATMAP.png"), plot=pmzmat,width = 11, height = 8.5)
ggsave(filename = paste0(outfile_name, "_residual_HEATMAP_CI.png"), plot=pm_ci,width = 11, height = 7)
ggsave(filename = paste0(outfile_name, "_normalized_HEATMAP_CI.png"), plot=pm_migg,width = 11, height = 7)
#dev.off()

#dev.off()


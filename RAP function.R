#read data

#setwd("~/Documents/Studie/Major/Artikel Scriptie/DATA/v19 mQTLs") #mac
setwd("~/Genetics/Artikel Scriptie/DATA/v19 mQTLs") #windows

library(RACER)
library(readr)
library(ggplot2)
library(dplyr)

mQTLs_shear_stress_v20 <- read_csv("mQTLs_shear_stress_nominal_v19_filtered_Illumina_v2.csv")

mQTLs_shear_stress_nominal_v20 <- mQTLs_shear_stress_v20[order(mQTLs_shear_stress_v20$Nominal_P), ]

TEST <- mQTLs_shear_stress_nominal_v20[!duplicated(mQTLs_shear_stress_nominal_v20$VARIANT), ]

#Regional association plot 
#enter the columns of your chromosome, chromosome position, and p-value column

mQTLs_plot_new <- RACER::formatRACER(assoc_data = TEST, 
                                     chr_col = 3, 
                                     pos_col = 4, 
                                     p_col = 32)

#Add LD data (if you have a lead SNP (SNP with highest p-value), you can look for SNPs in LD)

mQTLs_plot_new_2 <- RACER::ldRACER(assoc_data = mQTLs_plot_new, 
                                   rs_col = 2, pops = "EUR", 
                                   lead_snp = "rs7235957") 

#plotting data

#personalized function to change plot layour and colours (just run this once to add function to your environment)

rasplot <- function (assoc_data, chr, build = "hg19", set = "protein_coding", 
          plotby, gene_plot = NULL, snp_plot = NULL, start_plot = NULL, 
          end_plot = NULL, label_lead = FALSE) 
{
  if (missing(assoc_data)) {
    stop("Please provide a data set to plot.")
  }
  else if (missing(chr)) {
    stop("Please specify which chromosome you wish to plot.")
  }
  else if (missing(plotby)) {
    stop("Please specify the method by which you wish to plot.")
  }
  else if (plotby == "gene") {
    if (is.null(gene_plot)) {
      stop("Please specify a gene to plot by.")
    }
  }
  else if (plotby == "snp") {
    if (is.null(snp_plot)) {
      stop("Please specify a snp to plot by.")
    }
  }
  else if (plotby == "coord") {
    if (is.null(start_plot) | is.null(end_plot)) {
      stop("Please specify start coordinate for plot.")
    }
  }
  else {
    message("All inputs are go.")
  }
  reqs = c("CHR", "POS", "LOG10P")
  cols = colnames(assoc_data)
  if (sum(reqs %in% cols) == 3) {
  }
  else {
    stop("Association Data Set is missing a required column, please format your data set using formatRACER.R.")
  }
  reqs_2 = c("LD", "LD_BIN")
  if (sum(reqs_2 %in% cols) == 2) {
  }
  else {
    message("Association Data Set is missing LD data, the resulting plot won't have LD information, but you can add it using the ldRACER.R function.")
  }
  `%>%` <- magrittr::`%>%`
  if (build == "hg38") {
    utils::data(hg38)
    chr_in = chr
    colnames(hg38) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", 
                       "LENGTH", "GENE_NAME", "TYPE")
    gene_sub = hg38[hg38$CHR == chr_in, ]
  }
  else if (build == "hg19") {
    utils::data(hg19)
    chr_in = chr
    colnames(hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", 
                       "LENGTH", "GENE_NAME", "TYPE")
    gene_sub = hg19[hg19$CHR == chr_in, ]
  }
  if (set == "protein_coding") {
    gene_sub = gene_sub[gene_sub$TYPE == "protein_coding", 
                        ]
  }
  else {
    gene_sub = gene_sub
  }
  if (sum(is.null(plotby)) == 1) {
    stop("Please specify a method by which to plot.")
  }
  if (sum(is.null(plotby)) == 0) {
    message("Plotting by...")
    if ((plotby == "coord") == TRUE) {
      message("coord")
      start = start_plot
      end = end_plot
    }
    else if ((plotby == "gene") == TRUE) {
      message(paste("gene:", gene_plot))
      if (sum(is.null(gene_plot)) == 0) {
        p = subset(gene_sub, gene_sub$GENE_NAME == gene_plot)
        start = min(p$TRX_START) - 5e+05
        end = max(p$TRX_END) + 5e+05
      }
      else {
        message("No gene specified.")
      }
    }
    else if ((plotby == "snp") == TRUE) {
      message(paste("snp", snp_plot))
      q = assoc_data[assoc_data$RS_ID == snp_plot, ]
      w = q$POS
      w = as.numeric(as.character(w))
      start = w - 5e+05
      end = w + 5e+05
    }
  }
  gene_sub = subset(gene_sub, gene_sub$TRX_START > (start - 
                                                      50000))
  gene_sub = subset(gene_sub, gene_sub$TRX_END < (end + 50000))
  gene_sub = gene_sub[, c(3, 4, 6)]
  gene_sub = reshape2::melt(gene_sub, id.vars = "GENE_NAME")
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
  plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")
  message("Reading in association data")
  in.dt <- as.data.frame(assoc_data)
  in.dt$POS = as.numeric(as.character(in.dt$POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter_(in.dt, ~CHR == chr_in)
  in.dt = dplyr::filter_(in.dt, ~POS > start) %>% dplyr::filter_(~POS < 
                                                                   end)
  if (label_lead == TRUE) {
    lsnp_row = which(in.dt$LABEL == "LEAD")
    label_data = in.dt[lsnp_row, ]
    if (dim(label_data)[1] == 0) {
      lsnp_row = in.dt[in.dt$LOG10P == max(in.dt$LOG10P), 
                       ]
      label_data = lsnp_row[1, ]
    }
  }
  message("Generating Plot")
  
  uithof_color = c("#FBB820","#F59D10","#E55738","#DB003F","#E35493","#D5267B",
                   "#CC0071","#A8448A","#9A3480","#8D5B9A","#705296","#686AA9",
                   "#6173AD","#4C81BF","#2F8BC9","#1290D9","#1396D8","#15A6C1",
                   "#5EB17F","#86B833","#C5D220","#9FC228","#78B113","#49A01D",
                   "#595A5C","#A2A3A4", "#D7D8D7", "#ECECEC", "#FFFFFF", "#000000")
  
  if ("LD" %in% colnames(in.dt) && "LD_BIN" %in% colnames(in.dt)) {
    c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", 
                                                      y = "y_value")) + 
      geom_line(aes_string(group = "GENE_NAME"), 
                size = 2.5, 
                color = "black",
                alpha = 0.8) +
      geom_line(ggplot2::aes_string(group = "GENE_NAME"), 
                         size = 2, 
                         color = uithof_color[23],
                         alpha = 0.6) +
      ggplot2::theme_minimal() + ggplot2::geom_text(data = plot_lab, 
                                               ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"), 
                                               hjust = -0.1, vjust = 0.3, size = 2) + 
      ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) + 
      ggplot2::coord_cartesian(xlim = c(start, end), ylim = c(0, (max(gene_sub$y_value) + 1))) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 26, vjust = 1, hjust = 1),
                     axis.title.x = element_text(size = 10),
                     axis.text.y = ggplot2::element_blank(), 
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.grid.major.y = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major.x = element_line(colour = "#DEDEDE", size = 0.5, linetype = "longdash", lineend = "round"),
                     axis.line.x = element_line(colour = uithof_color[30]),
                     axis.text.x = element_text(colour = uithof_color[30], size = 6),
                     plot.margin=unit(c(0.2,0,0.2,0.2), "cm"))
    
    
    
    b = ggplot2::ggplot(in.dt, ggplot2::aes_string(x = "POS", 
                                                   y = "LOG10P", color = "LD_BIN")) + ggplot2::geom_point(size = 2) + 
      ggplot2::scale_colour_manual(expression(paste(R^2, sep = "")), values = c(`1.0-0.8` = uithof_color[11], 
                                              `0.8-0.6` = uithof_color[2], `0.6-0.4` = uithof_color[8], 
                                              `0.4-0.2` = uithof_color[24], `0.2-0.0` = uithof_color[18],
                                              `NA` = "grey"), drop = FALSE) + ggplot2::theme_minimal() +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
                     legend.text = element_text(size = 6),
                     legend.title.align = 0.15,
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major.y = element_line(colour = "#DEDEDE", size = 0.5, linetype = "longdash", lineend = "round"),
                     axis.line = element_line(colour = uithof_color[30]),
                     axis.text.x = element_blank(), #axis.text.x = element_text(colour = uithof_color[30], size = 6)
                     axis.text.y = element_text(colour = uithof_color[30], size = 6),
                     axis.title.y = element_text(size = 10),
                     axis.title.x = element_blank(),
                     plot.margin=unit(c(0.2,0,-0.2,0.2), "cm")) + #add more stuff!!!
      ggplot2::xlab("Chromosome Position") + ggplot2::ylab("-log10(p-value)") + 
      ggplot2::coord_cartesian(xlim = c(start, end), ylim = c(min(in.dt$LOG10P), 
                                                              max(in.dt$LOG10P)))
  }
  else {
    c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", 
                                                      y = "y_value")) + ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), 
                                                                                           size = 2) + ggplot2::theme_bw() + ggplot2::geom_text(data = plot_lab, 
                                                                                                                                                ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"), 
                                                                                                                                                hjust = -0.1, vjust = 0.3, size = 2.5) + ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", 
                                                                                                                                                                                                                                             size = 28), axis.text.y = ggplot2::element_blank(), 
                                                                                                                                                                                                        axis.ticks.y = ggplot2::element_blank()) + ggplot2::xlab(paste0("Chromosome ", 
                                                                                                                                                                                                                                                                        chr_in, " Position")) + ggplot2::coord_cartesian(xlim = c(start, 
                                                                                                                                                                                                                                                                                                                                  end), ylim = c(0, (max(gene_sub$y_value) + 1)))
    b = ggplot2::ggplot(in.dt, ggplot2::aes_string(x = "POS", 
                                                   y = "LOG10P")) + ggplot2::geom_point() + ggplot2::theme_bw() + 
      ggplot2::xlab("Chromosome Position") + ggplot2::ylab("-log10(p-value)") + 
      ggplot2::coord_cartesian(xlim = c(start, end), ylim = c(min(in.dt$LOG10P), 
                                                              max(in.dt$LOG10P)))
  }
  if (label_lead == TRUE) {
    b = b + ggplot2::geom_point(data = label_data, aes_string(x = "POS", 
                                                     y = "LOG10P"), color = "black", size = 2.5)
    b = b + ggplot2::geom_text(data = label_data, aes_string(label = "RS_ID"), 
                      color = "black", size = 3, hjust = 1.05, nudge_y = 1.5, nudge_x = 2)
  }
  ggpubr::ggarrange(b, c, heights = c(3, 1), nrow = 2, ncol = 1, 
                    common.legend = TRUE, legend = "right")
}




ENOSF1_region_plot <- rasplot(assoc_data = mQTLs_plot_new_2, 
                       chr = 18, 
                       build = "hg19", #genome build
                       plotby = "coord", 
                       start_plot = 463009, #region you are interested in
                       end_plot = 965639,
                       label_lead = T) #label your lead SNP


plot(ENOSF1_region_plot)



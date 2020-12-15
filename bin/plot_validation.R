#!/usr/bin Rscript
#Call the script like this: Rscript plot_validation.R "firm5_1" perc_id.txt 

options(warn=1)
rm(list=ls())
require("ggplot2")


density_plot <- function(infile, sdp) {
	 data <- read.table(infile, h=T)
	 df_1 <- subset(data, data$Closest_SDP == sdp)
	 df_2 <- subset(data, data$Closest_SDP != sdp)
	 length_sdp <- length(df_1[,1])
	 length_other <- length(df_2[,1])
	 tot_length = length_sdp + length_other
	 plot_title <- paste(sdp, ". Tot nb. orfs: ",tot_length, sep="")
	 all_perc <- c(df_1$Perc_id, df_2$Perc_id)
	 all_names <- c(rep(sdp,length_sdp), rep("other", length_other))
	 melt.df <- data.frame(all_perc, all_names)
	 melt.df$all_names <- factor(melt.df$all_names, levels=c(sdp, "other"))
	 p <- ggplot() + geom_density(data=melt.df, aes(x=all_perc, group=all_names, fill=all_names), alpha=0.5)+ xlim(60,100) +xlab("% identity") + ggtitle(plot_title) + scale_fill_manual(values=c("blue","grey")) +theme(legend.position="none") + geom_vline(xintercept=95, linetype="dashed")
	 return(p)	 
}

##data acquisition

Args <- commandArgs(trailingOnly=TRUE)
sdp <- Args[1]
infile <- Args[2]
p <- density_plot(infile,sdp)
outfile <- paste0(sdp,"_recruitment_plot.pdf")
cairo_pdf(outfile, h=9.5,w=7.3)
p
dev.off()


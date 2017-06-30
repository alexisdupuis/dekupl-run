#!/usr/bin/Rscript

setwd("/home/alexis/Documents/dekupl-run/")

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop("Two arguments must be supplied (input file and output name).\n", call.=FALSE)
}

sample <- read.table(args[1], sep="\t", dec=".", strip.white=TRUE)

library(reshape)

sample <- rename(sample, c(V1="chr", V2="locus", V3="depth"))

library(ggplot2)

plot <- ggplot(data=sample, aes(x=locus, y=depth, colour=chr, shape=chr)) + geom_smooth(method="loess", span=0.3)
# plot <- plot + geom_point(size=1, alpha=0.3)
plot <- plot + geom_vline(aes(xintercept=3902350), color="green", linetype="dashed", size=1)
plot <- plot + geom_text(aes(x=3918000, label="\nContig position", y=25), colour="green", angle=0)
plot <- plot + scale_color_manual(values=c('Y-healthy'="yellow", 'Y-sick'="purple", 'Y-no-contig'="red", 'Y-contig'="blue"))
plot <- plot + ggtitle("Profondeur de couverture par position - ChrY")
plot <- plot + ylim(15, 25)

ggsave(plot=plot, file=args[2])
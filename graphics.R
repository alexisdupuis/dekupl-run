#!/usr/bin/Rscript

setwd("/home/alexis/Documents/dekupl-run/")

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop("Two arguments must be supplied (input file and output name).\n", call.=FALSE)
}

sample <- read.table(args[1], sep=" ", dec=".", strip.white=TRUE)

library(reshape)

sample <- rename(sample, c(V1="chr", V2="locus", V3="depth"))

library(ggplot2)

plot <- ggplot(data=sample, aes(x=locus, y=depth, colour=chr))
plot <- plot + geom_point(size=2)
plot <- plot + ggtitle("Profondeur de couverture par position - ChrY ; individu 12227")
plot <- plot + stat_smooth(colour="blue")
plot <- plot + ylim(0, 60)

ggsave(plot=plot, file=args[2])
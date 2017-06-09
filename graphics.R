#!/usr/bin/Rscript

setwd("/home/alexis/Documents/dekupl-run/")

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop("Two arguments must be supplied (input file and output name).\n", call.=FALSE)
}

sample <- read.table(args[1], sep=" ", dec=".", strip.white=TRUE)

library(reshape)

sample <- rename(sample, c(V1="chr", V2="locus", V3="depth"))

# xyplot(depth ~ locus, type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), scales=list(x=list(relation='same'), y=list(relation='same')), data=sample, main="depth by locus - ChrY")

library(ggplot2)

# plot <- ggplot(data=sample, aes(x=locus, y=depth, colour=chr, shape=chr))
# plot <- plot + geom_point(size=4)
# plot <- plot + ggtitle("depth by locus - ChrY")
# plot <- plot + geom_line(size=2)

# plot <- ggplot(data=sample, aes(x=locus, y=depth, fill=chr)) + geom_histogram(stat='identity')

plot <- ggplot(data=sample, aes(x=locus, y=depth, colour=chr))
plot <- plot + geom_point(size=2)
plot <- plot + ggtitle("depth by locus - ChrY ; 22382 (malade)")
plot <- plot + stat_smooth(colour="blue")
plot <- plot + ylim(0, 60)


ggsave(plot=plot, file=args[2])
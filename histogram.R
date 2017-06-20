#!/usr/bin/Rscript

setwd("/home/alexis/Documents/dekupl-run/")

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop("Two arguments must be supplied (input file and output name).\n", call.=FALSE)
}

sample <- read.table(args[1], sep=",", dec=".", strip.white=TRUE)

library(reshape)

sample <- rename(sample, c(V1="name", V2="nb_kmers"))

library(ggplot2)


plot <- ggplot(data=sample, aes(x=nb_kmers)) + geom_histogram(color="darkblue", fill="lightblue", breaks=seq(200000000, 1200000000, by=18000000))
plot <- plot + geom_vline(aes(xintercept=mean(nb_kmers)), color="blue", linetype="dashed", size=1)
plot <- plot + geom_vline(aes(xintercept=625000000), color="red", linetype="dashed", size=1)
plot <- plot + geom_text(aes(x=mean(nb_kmers)+45000000, label="\nMean", y=26), colour="blue", angle=0)
plot <- plot + ggtitle("Répartition des individus selon le nombre total de k-mers comptés")

ggsave(plot=plot, file=args[2])
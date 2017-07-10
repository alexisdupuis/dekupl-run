#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop("Two arguments must be supplied (input file and output name).\n", call.=FALSE)
}

setwd("/home/alexis/Documents/dekupl-run/")

sample <- read.table(args[1], sep="\t", dec=".", strip.white=TRUE)

library(reshape)

sample <- rename(sample, c(V1="chr", V2="locus", V3="depth", V4="sd"))

library(ggplot2)

# Affiche une ligne reliant tous les points obtenues avec les données
plot <- ggplot(data=sample, aes(x=locus, y=depth, group=chr, colour=chr, shape=chr)) + geom_line()

# Affiche un intervalle de confiance correspondant à la valeur associée à un point + ou - l'écart-type
plot <- plot + geom_ribbon(aes(ymin=(depth-sd), ymax=(depth+sd), fill=chr, linetype=NA), alpha=0.2)

# Trace une courbe de régression passant au plus proche de tous les points
# plot <- plot + geom_smooth(method="loess", span=0.3)

# Affiche une barre verticale représentant la position du contig
plot <- plot + geom_vline(aes(xintercept=3902350), color="green", linetype="dashed", size=1)

# plot <- plot + geom_text(aes(x=3918000, label="\nContig position", y=25), colour="green", angle=0)

# Change divers paramètres pour améliorer la lisibilité de la figure
plot <- plot + scale_color_manual(values=c('Y-healthy'="black", 'Y-sick'="brown", 'Y-no-contig'="red", 'Y-contig'="blue"))
plot <- plot + scale_fill_manual(values=c("black", "brown"))

# Ajoute un titre et modifie les limites de l'axe des ordonnées
plot <- plot + ggtitle("Profondeur de couverture par position - ChrY")
plot <- plot + ylim(0, 50)

# Enregistre la figure
ggsave(plot=plot, file=args[2])
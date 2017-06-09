#!/usr/bin/Rscript

setwd("/home/alexis/Documents/dekupl-run/")

ratios <- data.frame(sains=c(6.38, 13.91), malades=c(51.01, 23.96), row.names=c("c > 0","c = 0"))

chiSquared <- chisq.test(ratios)

chiSquared$observed

chiSquared

# chiSquared$expected

# chiSquared$residual
pkgs_list <- c("ggplot2", "ggpubr")

idx <- sapply(pkgs_list, function(x) x %in% installed.packages())
install.packages(pkgs_list[idx])

library(pkgs_list)
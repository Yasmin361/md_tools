library(ggplot2)
library(gridExtra)
library(grid)


#####bulk-reading####
setwd("~/Desktop/run002/analyses_complete/alpha_helicity")        #supply path here
rama <- read.csv("fingerloop_ramachandran.list", header = FALSE, sep = "")
avg_helicity <- read.csv("freq_fingerloop_helicity.list", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

fl_helicity <- data.frame(residues = c(sort(avg_helicity$V1,decreasing = TRUE)),
                  freq=c(avg_helicity$V3))
fl_rama <- data.frame(phi = c(rama$V1), psi = c(rama$V2), res = c(rama$V3))

fl_helicity_plot <- ggplot(fl_helicity, aes(x=factor(rev(residues), labels = c("TYR-63", "GLY-64", "ARG-65", "GLU-66", 
                                                                          "ASP-67", "LEU-68", "ASP-69", "VAL-70", 
                                                                          "LEU-71", "GLY-72", "LEU-73", "THR-74", 
                                                                          "PHE-75", "ARG-76", "LYS-77", "ASP-78")), y=freq)) + 
  geom_line(aes(group = 1)) + geom_point() +
  ylab("Average helicity of finger loop residues for last 500 ns of all simulations") +
  xlab("ßarr1 residues") + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=14)) + 
  labs(labels=as.character(avg_helicity$V1)) + coord_flip()
fl_helicity_plot

fl_rama_plot <- ggplot(subset(fl_rama, (res == "ASP-70" 
                                       # & phi<=-35 & phi>=-95 & psi<=-10 & psi>=-70
                                        )), aes(x=phi, y=psi)) + geom_point(alpha = 0.05) + ylim(c(-180,180)) + xlim(c(-180,180))
fl_rama_plot
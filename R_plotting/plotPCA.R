library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#####bulk-reading####
setwd("~/Desktop/run002/analyses_complete/pca")        #supply path
temp <- list.files(pattern = "*cAlpha_proj1v2.xvg", recursive = TRUE) #list files in subfolders
pca_arr <- lapply(temp,function(i){
  read.csv(i,header = FALSE,skip = 17,sep = "",stringsAsFactors = FALSE)         #read with skipping gromacs output header
})
names(pca_arr) <- gsub("/","-",list.files(pattern = "*cAlpha_proj1v2.xvg", recursive = TRUE)) #strip parent folder names
plots <- list()
for (i in c(1:5)) {
  pca_arr[[i]] <- cbind(pca_arr[[i]], seq(0,length(pca_arr[[i]]$V1)-1)) 
  colnames(pca_arr[[i]]) <- c("eigv1","eigv2","runtime")
  plots[[i]] <- ggplot(pca_arr[[i]],aes(x=as.numeric(eigv1),y=as.numeric(eigv2),color=runtime/5)) + 
    geom_point(shape=".") + scale_color_viridis_c() + labs(title= paste("Rep", i,sep = " "), color = "time [ns]") + 
    ylim(c(-11,11)) + xlim(c(-11,11)) + coord_fixed() + xlab("eigenvector 1") + ylab("eigenvector 2") +
    theme(legend.text=element_text(size=14), legend.title=element_text(size = 14))
  legend <- get_legend(plots[[i]])
  plots[[i]] <- plots[[i]] + theme(legend.position = "none")
}

overview <- grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]], legend,
                         nrow=2, ncol=3 ,top=textGrob("Projection of first two eigenvectors of trajectories", 
                                                      gp=gpar(fontsize=16)))
overview


####spectrum###
setwd("~/Desktop/run002/analyses_complete/pca/step7_production_ext_all_complete_200ps_clustered")        #supply path
evspectrum <- read.csv("cAlpha/all_cAlpha_ev_spectrum.xvg",header = FALSE, skip = 17,sep = "",stringsAsFactors = FALSE)
evspectrum$V1 <- as.numeric(evspectrum$V1)
evspectrum$V2 <- as.numeric(evspectrum$V2)
sum(evspectrum$V2[1:3])/sum(evspectrum$V2)
eva <- ggplot(evspectrum,aes(x=V1,y=V2)) + geom_bar(stat="identity")
eva
####unified trj####

setwd("~/Desktop/run002/analyses_complete/pca/step7_production_ext_all_complete_200ps_clustered")        #supply path
temp <- list.files(pattern = "*cAlpha_proj1v2.xvg", recursive = TRUE) #list files in subfolders
pca_all <- lapply(temp,function(i){
  read.csv(i,header = FALSE,skip = 17,sep = "",stringsAsFactors = FALSE)         #read with skipping gromacs output header
})
names(pca_all) <- gsub("/","-",list.files(pattern = "*cAlpha_proj1v2.xvg", recursive = TRUE)) #strip parent folder names
expnames <- c(rep("Rep3",14956),rep("Rep2",14945),rep("Rep5",14930),rep("Rep4",14922),rep("Rep1",14918))
exptimes <- c(seq(1,14956),seq(1,14945),seq(1,14930),seq(1,14922),seq(1,14918))
pca_all[[1]] <- cbind(pca_all[[1]],expnames,exptimes)
head(pca_all[[1]])
all <- ggplot(pca_all[[1]],aes(x=as.numeric(V1),y=as.numeric(V2),color=expnames)) + 
  geom_point(shape=".",alpha=0.5) + 
  labs(color = "Trajectories") + 
  coord_fixed() + xlab("eigenvector 1") + ylab("eigenvector 2") +
  theme(legend.text=element_text(size=8), legend.title=element_text(size = 10), 
        axis.title = element_text(size=10), axis.text = element_text(size=8)) + 
  scale_color_manual(values=rainbow(5)) + guides(color = guide_legend(override.aes = list(shape = 19)))


ggsave("20201204_alltry_PCA_ev1_vs_ev2.png", plot=all, device=png(), 
       path="~/Desktop/figures/", width=16, height = 16, units = "cm", dpi=300)
dev.off()

alltime <- ggplot(pca_all[[1]],aes(x=as.numeric(V1),y=as.numeric(V2),color=exptimes/5)) + 
  geom_point(shape=".",alpha=0.5)  + labs(color = "time [ns]") + 
  coord_fixed() + xlab("eigenvector 1") + ylab("eigenvector 2") +
  theme(legend.text=element_text(size=8), legend.title=element_text(size = 10), 
        axis.title = element_text(size=10), axis.text = element_text(size=8)) +   
  scale_color_viridis_c(limits=c(1,3000))

ggsave("20201204_alltry_PCA_ev1_vs_ev2_bytime.png", plot=alltime, device=png(), 
       path="~/Desktop/figures/", width=16, height = 16, units = "cm", dpi=300)
dev.off()
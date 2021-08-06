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

setwd("~/Desktop/run002/analyses_complete/VMD_stride")        #supply path
temp <- list.files() #list files in subfolders
stride <- lapply(temp,function(i){
  read.csv(i,header = FALSE,skip = 9,sep = " ",stringsAsFactors = FALSE)         #read with skipping gromacs output header
})
names(stride)<- temp

rep1 <- data.frame(cbind(stride[[1]]$V1,stride[[1]]$V4,stride[[1]]$V5))
colnames(rep1) <- c("residue","runtime","secstr")
rep2 <- data.frame(cbind(stride[[2]]$V1,stride[[2]]$V4,stride[[2]]$V5))
colnames(rep2) <- c("residue","runtime","secstr")
rep3 <- data.frame(cbind(stride[[3]]$V1,stride[[3]]$V4,stride[[3]]$V5))
colnames(rep3) <- c("residue","runtime","secstr")
rep4 <- data.frame(cbind(stride[[4]]$V1,stride[[4]]$V4,stride[[4]]$V5))
colnames(rep4) <- c("residue","runtime","secstr")
rep5 <- data.frame(cbind(stride[[5]]$V1,stride[[5]]$V4,stride[[5]]$V5))
colnames(rep5) <- c("residue","runtime","secstr")


transition1 <- ggplot(rep1, aes(x=as.numeric(runtime), y=residue, fill= secstr)) + 
  geom_tile() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  ylab("Rep 1") + xlab("") +
  scale_y_discrete(labels= c("Y63",rep("",6),"V70",rep("",7),"D78")) +  
  scale_fill_brewer(palette="Accent", name = "Secondary Structure", 
                    labels= c("isolated bridge","coil","extended","3-10 helix","\u03b1-helix","turn")) + 
  theme(legend.text = element_text(size=8),legend.title = element_text(size=10),
        legend.key.size = unit(.7,"line")) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))
legend1 <- get_legend(transition1)
transition1 <- transition1 + theme(legend.position = "none")

transition2 <- ggplot(rep2, aes(x=as.numeric(runtime), y=residue, fill= secstr)) + 
  geom_tile() +
  scale_fill_brewer(palette="Accent", name = "Secondary Structure", labels= c("isolated bridge","coil","extended","3-10 helix","\u03b1-helix","turn")) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  ylab("Rep 2") + theme(legend.position = "none") +  xlab("") +
  scale_y_discrete(labels= c("Y63",rep("",6),"V70",rep("",7),"D78")) +  
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))

transition3 <- ggplot(rep3, aes(x=as.numeric(runtime), y=residue, fill= secstr)) + 
  geom_tile() +
  scale_fill_brewer(palette="Accent", name = "Secondary Structure", labels= c("isolated bridge","coil","extended","3-10 helix","\u03b1-helix","turn")) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  ylab("Rep 3") + theme(legend.position = "none") +  xlab("Time [ns]") +
  scale_y_discrete(labels= c("Y63",rep("",6),"V70",rep("",7),"D78")) +  
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))

transition4 <- ggplot(rep4, aes(x=as.numeric(runtime), y=residue, fill= secstr)) + 
  geom_tile() +
  scale_fill_brewer(palette="Accent", name = "Secondary Structure", labels= c("isolated bridge","coil","extended","3-10 helix","\u03b1-helix","turn")) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  ylab("Rep 4") + theme(legend.position = "none") + xlab("") +
  scale_y_discrete(labels= c("Y63",rep("",6),"V70",rep("",7),"D78")) +  
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))

transition5 <- ggplot(rep5, aes(x=as.numeric(runtime), y=residue, fill= secstr)) + 
  scale_fill_brewer(palette="Accent", name = "Secondary Structure", labels= c("isolated bridge","coil","extended","3-10 helix","\u03b1-helix","turn")) + 
  geom_tile() + ylab("Rep 5") + theme(legend.position = "none") + xlab("Time [ns]") +
  scale_y_discrete(labels= c("Y63",rep("",6),"V70",rep("",7),"D78")) +  
  scale_x_continuous(expand = c(0.01,0.01))
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))

# lay <- rbind(c(1,1,1,6),
#              c(2,2,2,6),
#              c(3,3,3,6),
#              c(4,4,4,6),
#              c(5,5,5,6))
# 
# overview_transitions <- grid.arrange(arrangeGrob(transition1,transition2,transition3,transition4,transition5,nrow=5),legend1,
#                                      ncol=2, widths=c(0.8,0.2), top="Secondary Structure Transitions of Finger Loop")
emptyt <- textGrob("")
  
overview_transitions <- grid.arrange(arrangeGrob(transition1,transition2,transition3,nrow=3),
                                     emptyt,
                                     arrangeGrob(transition4,transition5,legend1,nrow=3),
                                     emptyt,
                                     widths=c(0.48,0.02,0.48,0.02))
ggsave("20201109_secondary_structure_transitions_fingerloop.png", plot=overview_transitions, device=png(), 
       path="~/Desktop/figures/", width=16, height = 10.6, units = "cm", dpi=300)
dev.off()


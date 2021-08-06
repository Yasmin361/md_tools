library(ggplot2)
library(gridExtra)
library(grid)

#####bulk-reading####
setwd("~/Desktop/run002/analyses_complete/density")        #supply path here
temp <- list.files(pattern = "*.xvg", recursive = TRUE) #list files in subfolders
densitydis <- lapply(temp,function(i){
  read.csv(i,header = FALSE,skip = 24,sep = "")         #read with skipping gromacs output header
})
names(densitydis) <- basename(list.files(pattern = "*.xvg", recursive = TRUE)) #strip parent folder names
replicats <- 5
expnum <- length(densitydis)/replicats

index <- 0
for (j in 1:replicats) {
  for (i in 1:expnum) {
    densitydis[[index+i]][["run"]] <- rep(as.character(j),length(densitydis[[index+i]]$V1))  #add feature run number for tidy plotting
    densitydis[[index+i]][["comp"]] <- rep(as.character(i),length(densitydis[[index+i]]$V1))  #add feature run number for tidy plotting
  }
  index <- index + expnum
}

water <- rbind(densitydis$density_rep1_waterL.xvg,
               densitydis$density_rep2_waterL.xvg,
               densitydis$density_rep3_waterL.xvg,
               densitydis$density_rep4_waterL.xvg,
               densitydis$density_rep5_waterL.xvg)

protein <- rbind(densitydis$density_rep1_protL.xvg,
               densitydis$density_rep2_protL.xvg,
               densitydis$density_rep3_protL.xvg,
               densitydis$density_rep4_protL.xvg,
               densitydis$density_rep5_protL.xvg)

lipids <- rbind(densitydis$density_rep1_lipidsL.xvg,
                densitydis$density_rep2_lipidsL.xvg,
                densitydis$density_rep3_lipidsL.xvg,
                densitydis$density_rep4_lipidsL.xvg,
                densitydis$density_rep5_lipidsL.xvg)

cedge <- rbind(densitydis$density_rep1_cedgecAlphaL.xvg,
               densitydis$density_rep2_cedgecAlphaL.xvg,
               densitydis$density_rep3_cedgecAlphaL.xvg,
               densitydis$density_rep4_cedgecAlphaL.xvg,
               densitydis$density_rep5_cedgecAlphaL.xvg)

phosphates <- rbind(densitydis$density_rep1_headgroupsL.xvg,
                    densitydis$density_rep2_headgroupsL.xvg,
                    densitydis$density_rep3_headgroupsL.xvg,
                    densitydis$density_rep4_headgroupsL.xvg,
                    densitydis$density_rep5_headgroupsL.xvg)

####phosphate headgroup maximum####
phosphatesll <- subset(phosphates,phosphates$V1<0)
phosphatesll["product"] <- phosphatesll$V1*10*phosphatesll$V2
ll <- sum(phosphatesll$product)/sum(phosphatesll$V2)

phosphatesul <- subset(phosphates,phosphates$V1>0)
phosphatesul["product"] <- phosphatesul$V1*10*phosphatesul$V2
ul <- sum(phosphatesul$product)/sum(phosphatesul$V2)

#calculate max of distribution
repmax <- subset(cedge,cedge$run=="2")
repmax["product"] <- repmax$V1*10*repmax$V2
repmaximum <- sum(repmax$product)/sum(repmax$V2)
repmaximum

####plots####
sys_dis <- ggplot() +
  ggtitle("System components")+
  xlab("Distance from membrane center [\u00C5]") + ylab("Density [kg*m^-3]")  +
  geom_vline(xintercept=ll, linetype="dashed",size=0.5) +
  geom_vline(xintercept=ul, linetype="dashed",size=0.5)+
  geom_line(data=protein,aes(x=V1*10,y=V2,color=run,linetype=comp),size=0.5) +
  geom_line(data=water,aes(x=V1*10,y=V2,color=run,linetype=comp),size=0.5) + 
  scale_color_manual(values=rainbow(5)) +
  geom_line(data=lipids,aes(x=V1*10,y=V2,color=run,linetype=comp),size=0.5) + 
  geom_line(data=phosphates,aes(x=V1*10,y=V2,color=run,linetype=comp),size=0.5) + 
  labs(color="trajectories") + 
  scale_linetype_manual(name ="component",
                        values = c("dotdash","longdash","solid","dotted"), 
                        labels=c("headgroup P","POPC","protein","water")) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.text = element_text(size=8), title = element_text(size=10)) +
  coord_flip()
sys_dis <- arrangeGrob(sys_dis,top = textGrob("a", 
                                                  x=unit(0,"npc"), 
                                                  y = unit(1,"npc"),
                                                  just=c("left","top")))
cedge_dis <- ggplot() + 
  ggtitle("\u03B2arr-1 C-edge C\u03B1")+
  xlab("Distance from membrane center [\u00C5]") + ylab("Density [kg*m^-3]")  +
  geom_line(data=cedge,aes(x=V1*10,y=V2,color=run),size=0.5) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) + 
  theme(legend.text = element_text(size=8), title = element_text(size=10)) +
  geom_vline(xintercept=ll, linetype="dashed",size=0.5) +
  geom_vline(xintercept=ul, linetype="dashed",size=0.5)+
  coord_flip() + theme(legend.position = "none")
cedge_dis <- arrangeGrob(cedge_dis,top = textGrob("b", 
                                                  x=unit(0,"npc"), 
                                                  y = unit(1,"npc"),
                                                  just=c("left","top")))

all <- grid.arrange(sys_dis,cedge_dis,nrow=1,widths=c(0.65,0.35))
ggsave("20201119_densities.png", plot=all, device=png(), 
       path="~/Desktop/figures/", width=16, height = 8, units = "cm", dpi=300)
dev.off()
library(ggplot2)
library(gridExtra)
library(grid)
library(zoo)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#####bulk-reading####
setwd("~/Desktop/run002/analyses_complete/energy")        #supply path here
temp <- list.files(pattern = "*.xvg", recursive = TRUE) #list files in subfolders
ene_calc <- lapply(temp,function(i){
  read.csv(i,header = FALSE,skip = 24,sep = "")         #read with skipping gromacs output header
})
names(ene_calc) <- basename(list.files(pattern = "*.xvg", recursive = TRUE)) #strip parent folder names

replicats <- 5
expnum <- length(ene_calc)/replicats

index <- 0
for (j in 1:replicats) {
  for (i in 1:expnum) {
    ene_calc[[index+i]] <- cbind(ene_calc[[index+i]],run=rep(as.character(j),length(ene_calc[[index+i]]$V1)))  #add feature run number for tidy plotting
  }
  index <- index + expnum
}

index <- 0
for (j in 1:replicats) {
  for (i in 1:expnum) {
    ene_calc[[index+i]] <- cbind(ene_calc[[index+i]],exp=rep(as.character(i),length(ene_calc[[index+i]]$V1)))  #add feature exp number for tidy plotting
  }
  index <- index + expnum
}

ene_calc <- do.call("rbind", ene_calc)
ene_calc_reduced <- ene_calc[seq(1,nrow(ene_calc),1000),] # reduce from 10 ps to 10 ns steps
#y= rollmean(V2*10,k=100,na.pad = TRUE)-rollapply(V2*10, width=100, FUN=sd, fill=0, align="r")
allmeans <- c(mean(ene_calc[which(ene_calc$exp==1),]$V2*10),
              mean(ene_calc[which(ene_calc$exp==2),]$V2),
              mean(ene_calc[which(ene_calc$exp==3),]$V2),
              mean(ene_calc[which(ene_calc$exp==4),]$V2),
              mean(ene_calc[which(ene_calc$exp==5),]$V2),
              mean(ene_calc[which(ene_calc$exp==6),]$V2))
allsd <- c(sd(ene_calc[which(ene_calc$exp==1),]$V2*10),
           sd(ene_calc[which(ene_calc$exp==2),]$V2),
           sd(ene_calc[which(ene_calc$exp==3),]$V2),
           sd(ene_calc[which(ene_calc$exp==4),]$V2),
           sd(ene_calc[which(ene_calc$exp==5),]$V2),
           sd(ene_calc[which(ene_calc$exp==6),]$V2))

dev.off()
####PLOTS####
update_geom_defaults("smooth", list(size=.5))
boxx <- ggplot(ene_calc_reduced[which(ene_calc_reduced$exp==1),],aes(x = V1/1000, y = V2*10)) +
  xlab("Time [ns]") + ylab("X-Y Box size [\u00C5]") + 
  scale_color_manual(values=rainbow(5)) +
  geom_smooth(span=0.2,aes(color=run),se=FALSE) +
#  ggtitle(paste(as.character(round(allmeans[1],2)),"±",
#                as.character(round(allsd[1],2)))) +
  labs(color="trajectories") + 
  annotation_custom(textGrob(paste(as.character(round(allmeans[1],2)),"±",
                                   as.character(round(allsd[1],2))),
                             gp=gpar(fontsize=8),x=0.05,y=0.05,hjust=0)) +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10)) +
  theme(legend.position = "bottom", legend.title = element_text(size=8), 
        legend.text = element_text(size=10))
legend1 <- get_legend(boxx)
boxx <- boxx + theme(legend.position = "none")

dens <- ggplot(ene_calc_reduced[which(ene_calc_reduced$exp==2),],aes(x = V1/1000, y = V2)) +
  xlab("Time [ns]") + ylab("Density [kg/m³]") + 
  geom_smooth(span=0.2,aes(color=run),se=FALSE) +
  scale_color_manual(values=rainbow(5)) +
#  ggtitle(paste(as.character(round(allmeans[2],2)),"±",
#                as.character(round(allsd[2],2)))) +
  annotation_custom(textGrob(paste(as.character(round(allmeans[2],2)),"±",
                                   as.character(round(allsd[2],2))),
                             gp=gpar(fontsize=8),x=0.05,y=0.05,hjust=0)) + 
  labs(color="trajectories") + theme(legend.position = "none")  +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10))


ekin <- ggplot(ene_calc_reduced[which(ene_calc_reduced$exp==3),],aes(x = V1/1000, y = V2)) +
  xlab("Time [ns]") + ylab("Kinetic Energy [kJ/mol]") + 
  geom_smooth(span=0.2,aes(color=run),se=FALSE) +
  scale_color_manual(values=rainbow(5)) +
#  ggtitle(paste(as.character(round(allmeans[3],2)),"±",
#                as.character(round(allsd[3],2)))) +
  annotation_custom(textGrob(paste(as.character(round(allmeans[3],2)),"±",
                                   as.character(round(allsd[3],2))),
                             gp=gpar(fontsize=8),x=0.05,y=0.05,hjust=0)) + 
  labs(color="trajectories") + theme(legend.position = "none")  +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10))


epot <- ggplot(ene_calc_reduced[which(ene_calc_reduced$exp==4),],aes(x = V1/1000, y = V2)) +
  xlab("Time [ns]") + ylab("Potential Energy [kJ/mol]") + 
  geom_smooth(span=0.2,aes(color=run),se=FALSE) +
  scale_color_manual(values=rainbow(5)) +
#  ggtitle(paste(as.character(round(allmeans[4],2)),"±",
#                as.character(round(allsd[4],2)))) +
  annotation_custom(textGrob(paste(as.character(round(allmeans[4],2)),"±",
                                   as.character(round(allsd[4],2))),
                             gp=gpar(fontsize=8),x=0.05,y=0.05,hjust=0)) + 
  labs(color="trajectories") + theme(legend.position = "none")  +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10))


pres <- ggplot(ene_calc_reduced[which(ene_calc_reduced$exp==5),],aes(x = V1/1000, y = V2)) +
  xlab("Time [ns]") + ylab("Pressure [bar]") + 
  geom_smooth(span=0.2,aes(color=run),se=FALSE) +
  scale_color_manual(values=rainbow(5)) +
#  ggtitle(paste(as.character(round(allmeans[5],2)),"±",
#                as.character(round(allsd[5],2)))) +
  annotation_custom(textGrob(paste(as.character(round(allmeans[5],2)),"±",
                                   as.character(round(allsd[5],2))),
                             gp=gpar(fontsize=8),x=0.05,y=0.05,hjust=0)) + 
  labs(color="trajectories") + theme(legend.position = "none")  +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10))


tempe <- ggplot(ene_calc_reduced[which(ene_calc_reduced$exp==6),],aes(x = V1/1000, y = V2)) +
  xlab("Time [ns]") + ylab("Temperature [K]") + 
  geom_smooth(span=0.2,aes(color=run),se=FALSE) +
  scale_color_manual(values=rainbow(5)) +
#  ggtitle(paste(as.character(round(allmeans[6],2)),"±",
#                as.character(round(allsd[6],2)))) +
  annotation_custom(textGrob(paste(as.character(round(allmeans[6],2)),"±",
                                   as.character(round(allsd[6],2))),
                             gp=gpar(fontsize=8),x=0.05,y=0.05,hjust=0)) + 
  labs(color="trajectories") + theme(legend.position = "none")  +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10))


overview_energy <- grid.arrange(arrangeGrob(tempe,pres,dens,boxx,epot,ekin,
                                            nrow=2, ncol=3), legend1,
                                nrow=2,heights=c(0.9,0.1))

#overview_energy <- grid.arrange(tempe,pres,dens,boxx,epot,ekin,nrow=2, ncol=3)


ggsave("20201128_quality_assurance.png", plot=overview_energy, device=png(), 
       path="~/Desktop/figures/", width=16, height = 10.6, units = "cm", dpi=300)


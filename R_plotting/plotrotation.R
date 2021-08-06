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

setwd("~/Desktop/run002/analyses_complete/VMD_rotation")        #supply path
temp <- list.files(pattern = "*.agr", recursive = TRUE) #list files in subfolders
rotation <- lapply(temp,function(i){
  read.csv(i,header = FALSE,skip = 3,sep = "")         #read with skipping gromacs output header
})
names(rotation) <- basename(list.files(pattern = "*.agr", recursive = TRUE)) #strip parent folder names
replicats <- 5
expnum <- length(rotation)/replicats

index <- 0
for (j in 1:replicats) {
  for (i in 1:expnum) {
    rotation[[index+i]] <- cbind(rotation[[index+i]],run=rep(as.character(j),length(rotation[[index+i]]$V1)))  #add feature run number for tidy plotting
  }
  index <- index + expnum
}

index <- 0
for (j in 1:replicats) {
  for (i in 1:expnum) {
    rotation[[index+i]] <- cbind(rotation[[index+i]],exp=rep(as.character(i),length(rotation[[index+i]]$V1)))  #add feature exp number for tidy plotting
  }
  index <- index + expnum
}

rotation <- do.call("rbind", rotation)
rotation$V1 <- as.numeric(rotation$V1)
rotation$V2 <- -as.numeric(rotation$V2)


rotation1 <- data.frame(cbind(rotation[which(rotation$run ==1 & rotation$exp == 1 & is.na(rotation$V2) == FALSE),]$V2,
                              rotation[which(rotation$run ==1 & rotation$exp == 2 & is.na(rotation$V2) == FALSE),]$V2,
                              rotation[which(rotation$run ==1 & rotation$exp == 3 & is.na(rotation$V2) == FALSE),]$V2
                   ))

rotation2 <- data.frame(cbind(rotation[which(rotation$run ==2 & rotation$exp == 1 & is.na(rotation$V2) == FALSE),]$V2,
                              rotation[which(rotation$run ==2 & rotation$exp == 2 & is.na(rotation$V2) == FALSE),]$V2,
                              rotation[which(rotation$run ==2 & rotation$exp == 3 & is.na(rotation$V2) == FALSE),]$V2
))

rotation3 <- data.frame(cbind(rotation[which(rotation$run ==3 & rotation$exp == 1 & is.na(rotation$V2) == FALSE),]$V2,
                              rotation[which(rotation$run ==3 & rotation$exp == 2 & is.na(rotation$V2) == FALSE),]$V2,
                              rotation[which(rotation$run ==3 & rotation$exp == 3 & is.na(rotation$V2) == FALSE),]$V2
))

rotation4 <- data.frame(cbind(rotation[which(rotation$run ==4 & rotation$exp == 1 & is.na(rotation$V2) == FALSE),]$V2,
                              rotation[which(rotation$run ==4 & rotation$exp == 2 & is.na(rotation$V2) == FALSE),]$V2,
                              rotation[which(rotation$run ==4 & rotation$exp == 3 & is.na(rotation$V2) == FALSE),]$V2
))

rotation5 <- data.frame(cbind(rotation[which(rotation$run ==5 & rotation$exp == 1 & is.na(rotation$V2) == FALSE),]$V2,
                              rotation[which(rotation$run ==5 & rotation$exp == 2 & is.na(rotation$V2) == FALSE),]$V2,
                              rotation[which(rotation$run ==5 & rotation$exp == 3 & is.na(rotation$V2) == FALSE),]$V2
))

#calculate mean/sd of theta angle
mean(rowMeans(rotation1))
sd(rowMeans(rotation1))

####plots####
ylimits <- c(-15,32)
activeangle <- mean(c(13.540,24.882,34.186))
startangle <- mean(c(15.147,20.625,20.339))

meanrot1 <- ggplot() +
  xlab("Time [µs]") +   theme(legend.position = "none", axis.title.y = element_blank()) +
  geom_hline(yintercept=0,size=0.5,linetype="dotted") +
  geom_hline(yintercept=activeangle,size=0.5,color="red")+
  geom_hline(yintercept=startangle, linetype="dashed",size=0.5,color="red") +
  geom_line(aes(x=(as.numeric(rownames(rotation1))+24)/1000,
                y=rollmean(rowMeans(rotation1), 200, fill = NA)), size=.5) +
  labs(color="trajectories") + ylim(ylimits) + ggtitle("Rep1") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10), 
        plot.title = element_text(size=10,hjust = 0.5)) 

meanrot2 <- ggplot() +
  xlab("Time [µs]") +   theme(legend.position = "none", axis.title.y = element_blank()) +
  geom_hline(yintercept=0,size=0.5,linetype="dotted") +
  geom_hline(yintercept=activeangle,size=0.5,color="red")+
  geom_hline(yintercept=startangle, linetype="dashed",size=0.5,color="red") +
  geom_line(aes(x=(as.numeric(rownames(rotation2))+20)/1000,
                y=rollmean(rowMeans(rotation2), 200, fill = NA)), size=.5) +
  labs(color="trajectories") + ylim(ylimits) + ggtitle("Rep2") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10), 
        plot.title = element_text(size=10,hjust = 0.5))

meanrot3 <- ggplot() +
  xlab("Time [µs]") +   theme(legend.position = "none", axis.title.y = element_blank()) +
  geom_hline(yintercept=0,size=0.5,linetype="dotted") +
  geom_hline(yintercept=activeangle,size=0.5,color="red")+
  geom_hline(yintercept=startangle, linetype="dashed",size=0.5,color="red") +
  geom_line(aes(x=(as.numeric(rownames(rotation3))+17)/1000,
                y=rollmean(rowMeans(rotation3), 200, fill = NA)), size=.5) +
  labs(color="trajectories") + ylim(ylimits) +  ggtitle("Rep3") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10), 
        plot.title = element_text(size=10,hjust = 0.5))

meanrot4 <- ggplot() +
  xlab("Time [µs]") +   theme(legend.position = "none", axis.title.y = element_blank()) +
  geom_hline(yintercept=0,size=0.5,linetype="dotted") +
  geom_hline(yintercept=activeangle,size=0.5,color="red")+
  geom_hline(yintercept=startangle, linetype="dashed",size=0.5,color="red") +
  geom_line(aes(x=(as.numeric(rownames(rotation4))+24)/1000,
                y=rollmean(rowMeans(rotation4), 200, fill = NA)), size=.5) +
  labs(color="trajectories") + ylim(ylimits) +  ggtitle("Rep4") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10), 
        plot.title = element_text(size=10,hjust = 0.5))

meanrot5 <- ggplot() +
  xlab("Time [µs]") +   theme(legend.position = "none", axis.title.y = element_blank()) +
  geom_hline(yintercept=0,size=0.5,linetype="dotted") +
  geom_hline(yintercept=activeangle,size=0.5,color="red")+
  geom_hline(yintercept=startangle, linetype="dashed",size=0.5,color="red") +
  geom_line(aes(x=(as.numeric(rownames(rotation5))+22)/1000,
                y=rollmean(rowMeans(rotation5), 200, fill = NA)), size=.5) +
  labs(color="trajectories") + ylim(ylimits) +  ggtitle("Rep5") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10), 
        plot.title = element_text(size=10,hjust = 0.5))


meanrot <- ggplot() +
  xlab("time in ns") + ylab("angle") + 
  geom_line(aes(x=as.numeric(rownames(rotation1))+24,
                y=rollmean(rowMeans(rotation1), 40, fill = NA)), size=.5, color=rainbow(5)[1]) +
  geom_line(aes(x=as.numeric(rownames(rotation2))+20,
                y=rollmean(rowMeans(rotation2), 40, fill = NA)), size=.5, color=rainbow(5)[2]) +
  geom_line(aes(x=as.numeric(rownames(rotation3))+17,
                y=rollmean(rowMeans(rotation3), 40, fill = NA)), size=.5, color=rainbow(5)[3]) +
  geom_line(aes(x=as.numeric(rownames(rotation4))+24,
                y=rollmean(rowMeans(rotation4), 40, fill = NA)), size=.5, color=rainbow(5)[4]) +
  geom_line(aes(x=as.numeric(rownames(rotation5))+22,
                y=rollmean(rowMeans(rotation5), 40, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  ylim(ylimits) +
  geom_hline(yintercept=activeangle, linetype="dashed",size=0.5)+
  geom_hline(yintercept=startangle, linetype="dotted",size=0.5)+
  geom_hline(yintercept=0,size=0.5)+
  annotation_custom(textGrob("active",x=0.94,y=0.85,hjust=0,gp=gpar(fontsize=8))) + 
  annotation_custom(textGrob("start",x=0.955,y=0.74,hjust=0,gp=gpar(fontsize=8))) +   
  annotation_custom(textGrob("inactive",x=0.92,y=0.38,hjust=0,gp=gpar(fontsize=8))) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))
meanrot

newplot <- grid.arrange(meanrot,legend,ncol=2,widths=c(0.9,0.15))

overview_rotation <- grid.arrange(meanrot1,meanrot2,meanrot3,meanrot4,meanrot5,ncol=5, left="\u03B")
ggsave("20210125_rotation_angle.png", plot=newplot, device=png(), 
       path="~/Desktop/figures/", width=16, height = 4, units = "cm", dpi=300)
dev.off()
ggsave("20210125_rotation_angle_single.png", plot=overview_rotation, device=png(), 
       path="~/Desktop/figures/", width=16, height = 4, units = "cm", dpi=300)
dev.off()


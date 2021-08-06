library(ggplot2)
library(gridExtra)
library(grid)
library(zoo)

#####bulk-reading####
setwd("~/Desktop/run002/analyses_complete/distance")        #supply path
temp <- list.files(pattern = "*.xvg", recursive = TRUE) #list files in subfolders
dist_calc <- lapply(temp,function(i){
  read.csv(i,header = FALSE,skip = 24,sep = "")         #read with skipping gromacs output header
})
names(dist_calc) <- basename(list.files(pattern = "*.xvg", recursive = TRUE)) #strip parent folder names
print(names(dist_calc))

#plots center of geometry distance between TMD and finger loop Calpha in angstrom over the course of the trajectory

update_geom_defaults("line",list(size=0.5))
rep1 <- ggplot(dist_calc$distavg_rep1_fingerloop_TMD.xvg, 
               aes(dist_calc$distavg_rep1_fingerloop_TMD.xvg$V1/1000000,
                   dist_calc$distavg_rep1_fingerloop_TMD.xvg$V2*10)) + 
  annotate("segment", x=0.75, xend=0.75, y=40, yend=38, 
           color="red",size=0.5,arrow=arrow(length = unit(0.1,"cm"))) +
  annotate("segment", x=2, xend=2, y=40, yend=38, 
           color="red",size=0.5,arrow=arrow(length = unit(0.1,"cm"))) +
  xlab("Time [µs]") + ylim(26,40) +
  geom_line(aes(y=rollmean(dist_calc$distavg_rep1_fingerloop_TMD.xvg$V2*10, 200, fill = NA))) + ggtitle("Rep1") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  theme(axis.title = element_text(size=10),axis.text = element_text(size=8), 
        plot.title = element_text(size=10,hjust = 0.5))


rep2 <- ggplot(dist_calc$distavg_rep2_fingerloop_TMD.xvg, 
               aes(dist_calc$distavg_rep2_fingerloop_TMD.xvg$V1/1000000,
                   dist_calc$distavg_rep2_fingerloop_TMD.xvg$V2*10)) +
  xlab("Time [µs]") + ylim(26,40) +
  geom_line(aes(y=rollmean(dist_calc$distavg_rep2_fingerloop_TMD.xvg$V2*10, 200, fill = NA))) + ggtitle("Rep2") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  theme(axis.title = element_text(size=10),axis.text = element_text(size=8), 
        plot.title = element_text(size=10,hjust = 0.5))



rep3 <- ggplot(dist_calc$distavg_rep3_fingerloop_TMD.xvg, 
               aes(dist_calc$distavg_rep3_fingerloop_TMD.xvg$V1/1000000,
                   dist_calc$distavg_rep3_fingerloop_TMD.xvg$V2*10)) + 
  xlab("Time [µs]") + ylim(26,40) +
  geom_line(aes(y=rollmean(dist_calc$distavg_rep3_fingerloop_TMD.xvg$V2*10, 200, fill = NA))) + ggtitle("Rep3") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  theme(axis.title = element_text(size=10),axis.text = element_text(size=8), 
        plot.title = element_text(size=10,hjust = 0.5))



rep4 <- ggplot(dist_calc$distavg_rep4_fingerloop_TMD.xvg, 
               aes(dist_calc$distavg_rep4_fingerloop_TMD.xvg$V1/1000000,
                   dist_calc$distavg_rep4_fingerloop_TMD.xvg$V2*10)) + 
  xlab("Time [µs]") + ylim(26,40) +
  geom_line(aes(y=rollmean(dist_calc$distavg_rep4_fingerloop_TMD.xvg$V2*10, 200, fill = NA))) + ggtitle("Rep4") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  theme(axis.title = element_text(size=10),axis.text = element_text(size=8), 
        plot.title = element_text(size=10,hjust = 0.5))



rep5 <- ggplot(dist_calc$distavg_rep5_fingerloop_TMD.xvg, 
               aes(dist_calc$distavg_rep5_fingerloop_TMD.xvg$V1/1000000,
                   dist_calc$distavg_rep5_fingerloop_TMD.xvg$V2*10)) + 
  xlab("Time [µs]") + ylim(26,40) +
  geom_line(aes(y=rollmean(dist_calc$distavg_rep5_fingerloop_TMD.xvg$V2*10, 200, fill = NA))) + ggtitle("Rep5") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  theme(axis.title = element_text(size=10),axis.text = element_text(size=8), 
        plot.title = element_text(size=10,hjust = 0.5))



overview_dist <- grid.arrange(rep1,rep2,rep3,rep4,rep5,
                             ncol=5, left="Distance [\u00C5]")

ggsave("20200928_fingerloop_distance.png", plot=overview_dist, device=png(), 
       path="~/Desktop/figures/", width=16, height = 4, units = "cm", dpi=300)
dev.off()
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
setwd("~/Desktop/run002/analyses_complete/rmsd")        #supply path here
temp <- list.files(pattern = "*.xvg", recursive = TRUE) #list files in subfolders
rmsd_calc <- lapply(temp,function(i){
  read.csv(i,header = FALSE,skip = 18,sep = "")         #read with skipping gromacs output header
})
names(rmsd_calc) <- basename(list.files(pattern = "*.xvg", recursive = TRUE)) #strip parent folder names
expnames <- names(rmsd_calc)[1:length(temp)]
replicats <- 5
expnum <- length(rmsd_calc)/replicats

index <- 0
for (j in 1:replicats) {
  for (i in 1:expnum) {
    rmsd_calc[[index+i]] <- cbind(rmsd_calc[[index+i]],run=rep(as.character(j),length(rmsd_calc[[index+i]]$V1)))  #add feature run number for tidy plotting
  }
  index <- index + expnum
}

index <- 0
for (j in 1:replicats) {
  for (i in 1:expnum) {
    rmsd_calc[[index+i]] <- cbind(rmsd_calc[[index+i]],exp=rep(as.character(i),length(rmsd_calc[[index+i]]$V1)))  #add feature exp number for tidy plotting
  }
  index <- index + expnum
}

rmsd_calc <- do.call("rbind", rmsd_calc)
expnames
############whole_complex###############################

repl_l <- ggplot(rmsd_calc[which(rmsd_calc$exp==1),],aes(x = V1, y = V2*10)) +
  xlab("") + ylab("RMSD [\u00C5] ---  \u03B2arr-1") + ylim(0,10) +
  geom_line(aes(color=run), size=0.5) + scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
legend <- get_legend(repl_l)
# 
# ntr_rmsd <- ggplot(rmsd_calc[which(rmsd_calc$exp==13),],aes(x = V1, y = V2*10)) + ylim(0,10) +
#   xlab("Time [ns]") + ylab("RMSD [\u00C5] --- NTS1R") +
#   geom_line(aes(color=run)) + scale_color_manual(values=rainbow(5)) +
#   labs(color="trajectories") +
#   theme(axis.text=element_text(size=10), axis.title=element_text(size=14))
# 
# overview <- grid.arrange(arr_rmsd,ntr_rmsd,nrow=2,top="C\u03B1-RMSD of NTS1R-\u03B2Arr1 complex for 5 x 3 µs simulations")
# overview

dev.off()
update_geom_defaults("line",list(size=.5))
arr_rmsd <- ggplot() +
  xlab("") + ylab("RMSD [\u00C5] - \u03B2arr-1") + ylim(0,10) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==1),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==1),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==1),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==1),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==1),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==1),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==1),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==1),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==1),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==1),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(rmsd_calc$run, values=rainbow(5)) +
  labs(color="trajectories") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))

ntr_rmsd <- ggplot() + ylim(0,10) +
  xlab("Time [ns]") + ylab("RMSD [\u00C5] - NTS1R") +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==13),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==15),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==13),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==15),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==13),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==15),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==13),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==15),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==13),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==15),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(as.factor(c(1:5)),values=rainbow(5)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))

overview <- grid.arrange(arrangeGrob(arr_rmsd,ntr_rmsd,nrow=2),legend,ncol=2, widths=c(0.88,0.12))
#                         top="C\u03B1-RMSD of \u03B2arr-1-NTS1R complex for 5 x 3 µs simulations"

ggsave("20201022_RMSD_arr_ntr_test.png", plot=overview, device=png(), 
       path="~/Desktop/figures/", width=16, height = 8, units = "cm", dpi=300)
overview <- grid.arrange(arrangeGrob(arr_rmsd,ntr_rmsd,nrow=2),legend,ncol=2, widths=c(0.88,0.12))
dev.off()




###########arrestin plots###########################################
arrylim = 14
ndom_rmsd <- ggplot() +
  xlab("") + ylab("N-domain") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==14),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==14),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==14),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==14),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==14),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==14),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==14),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==14),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==14),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==14),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) + 
  theme(legend.title = element_text(size=10),legend.text = element_text(size=8), 
        legend.spacing.y = unit(0.5,"cm")) + 
  guides(shape=guide_legend(override.aes = list(size=0.3))) + 
#  theme(legend.position = "bottom") +  
  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#ndom_rmsd
legendfeatures <- get_legend(ndom_rmsd)
ndom_rmsd <- ndom_rmsd + theme(legend.position = "none")

cdom_rmsd <- ggplot() +
  xlab("") + ylab("C-domain") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==4),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==4),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==4),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==4),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==4),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==4),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==4),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==4),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==4),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==4),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#cdom_rmsd

sheetI_rmsd <- ggplot() +
  xlab("") + ylab("\u03B2-sheet I") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==17),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==17),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==17),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==17),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==17),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==17),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==17),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==17),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==17),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==17),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#sheetI_rmsd

fingerL_rmsd <- ggplot() +
  xlab("") + ylab("finger loop") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==7),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==7),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==7),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==7),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==7),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==7),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==7),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==7),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==7),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==7),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories")+ scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#fingerL_rmsd

middleL_rmsd <- ggplot() +
  xlab("") + ylab("middle loop") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==13),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==13),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==13),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==13),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==13),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==13),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==13),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==13),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==13),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==13),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#middleL_rmsd

onesixtyL_rmsd <- ggplot() +
  xlab("") + ylab("160-loop") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==16),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==16),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==16),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==16),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==16),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==16),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==16),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==16),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==16),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==16),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#onesixtyL_rmsd

gateL_rmsd <- ggplot() +
  xlab("") + ylab("gate loop") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==8),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==8),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==8),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==8),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==8),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==8),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==8),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==8),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==8),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==8),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#gateL_rmsd

cedge_rmsd <- ggplot() +
  xlab("") + ylab("C-edge") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==5),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==5),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==5),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==5),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==5),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==5),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==5),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==5),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==5),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==5),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#cedge_rmsd

cloop_rmsd <- ggplot() +
  xlab("") + ylab("C-loop") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==6),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==6),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==6),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==6),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==6),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==6),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==6),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==6),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==6),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==6),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#cloop_rmsd

backloop_rmsd <- ggplot() +
  xlab("") + ylab("back loop") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==3),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==3),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==3),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==3),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==3),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==3),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==3),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==3),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==3),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==3),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#backloop_rmsd

helixI_rmsd <- ggplot() +
  xlab("") + ylab("\u03B1-helix I") + ylim(0,arrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==2),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==2),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==2),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==2),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==2),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==2),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==2),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==2),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==2),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==2),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
#helixI_rmsd

overview_arr <- grid.arrange(ndom_rmsd,cdom_rmsd,sheetI_rmsd,fingerL_rmsd,
                             helixI_rmsd,middleL_rmsd,onesixtyL_rmsd,
                             cloop_rmsd,gateL_rmsd, backloop_rmsd,cedge_rmsd,legendfeatures,
                             
                             ncol=3,top="C\u03B1-RMSD of \u03B2Arr1 structural features", 
                             )
overview_arr


#### NTS1R #####
ntrylim = 14
icl1_rmsd <- ggplot() +
  xlab("") + ylab("ICL1") + ylim(0,ntrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==10),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==10),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==10),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==10),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==10),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==10),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==10),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==10),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==10),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==10),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + xlab("") +  scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))


icl2_rmsd <- ggplot() +
  xlab("") + ylab("ICL2") + ylim(0,ntrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==11),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==11),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==11),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==11),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==11),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==11),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==11),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==11),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==11),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==11),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + xlab("") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))


icl3_rmsd <- ggplot() +
  xlab("") + ylab("ICL3") + ylim(0,ntrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==12),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==12),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==12),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==12),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==12),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==12),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==12),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==12),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==12),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==12),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") + xlab("") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))


h8_rmsd <- ggplot() +
  xlab("") + ylab("Helix 8") + ylim(0,ntrylim) +
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==9),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==1 & rmsd_calc$exp==9),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[1]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==9),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==2 & rmsd_calc$exp==9),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[2]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==9),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==3 & rmsd_calc$exp==9),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[3]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==9),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==4 & rmsd_calc$exp==9),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[4]) + 
  geom_line(aes(x=rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==9),]$V1,
                y=rollmean(rmsd_calc[which(rmsd_calc$run==5 & rmsd_calc$exp==9),]$V2*10, 200, fill = NA)), size=.5, color=rainbow(5)[5]) + 
  scale_color_manual(values=rainbow(5)) +
  labs(color="trajectories") +  scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) +
  theme(legend.position = "none")+  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))


overview_ntr <- grid.arrange(icl1_rmsd,icl2_rmsd,icl3_rmsd,h8_rmsd,
                             nrow=4,top="C\u03B1-RMSD of NTS1R structural features", 
                             bottom="Time [ns]")
overview_ntr

emptyt <- textGrob("")

overview_both <- grid.arrange(arrangeGrob(icl1_rmsd,icl2_rmsd,icl3_rmsd,h8_rmsd,nrow=4,top = "NTS1R"),
                              emptyt,
                              arrangeGrob(ndom_rmsd,cdom_rmsd,sheetI_rmsd,fingerL_rmsd,
                                          helixI_rmsd,middleL_rmsd,onesixtyL_rmsd,
                                          cloop_rmsd,gateL_rmsd, backloop_rmsd,
                                          cedge_rmsd,legendfeatures,
                                          nrow=4,top = "\u03B2arr-1"),ncol=3, widths=c(0.23,0.04,0.73))
ggsave("20201022_RMSD_arr_ntr_features_fullpage.png", plot=overview_both, device=png(), 
       path="~/Desktop/figures/", width=24.7, height = 14, units = "cm", dpi=300)
dev.off()

##########DUMP############
# old code for timetrace without smoothing
# arrylim = 17
# ndom_rmsd <- ggplot(rmsd_calc[which(rmsd_calc$exp==11),],aes(x = V1, y = V2*10)) +
#   xlab("") + ylab("N-domain") + ylim(0,arrylim) +
#   geom_line(aes(color=run)) + scale_color_manual(values=rainbow(5)) +
#   labs(color="trajectories") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
#   theme(axis.text=element_text(size=10), axis.title=element_text(size=14)) + 
#   theme(legend.position = "none") +  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
# #ndom_rmsd




#combined_arr[[i]] <<- rbind(eval(parse(text = paste("arr_parts_runs$","rmsd1_",as.character(i),".xvg",sep = ""))),
#                            eval(parse(text = paste("arr_parts_runs$","rmsd2_",as.character(i),".xvg",sep = ""))),
#                            eval(parse(text = paste("arr_parts_runs$","rmsd3_",as.character(i),".xvg",sep = ""))), 
#                            eval(parse(text = paste("arr_parts_runs$","rmsd4_",as.character(i),".xvg",sep = ""))),
#                            eval(parse(text = paste("arr_parts_runs$","rmsd5_",as.character(i),".xvg",sep = ""))), )
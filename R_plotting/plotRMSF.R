library(ggplot2)
library(gridExtra)
library(grid)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#####bulk-reading####
setwd("~/Desktop/run002/analyses_complete/rmsf")        #supply path
temp <- list.files(pattern = "rmsf_rep*", recursive = TRUE) #list files in subfolders
rmsf_calc <- lapply(temp,function(i){
  read.csv(i,header = FALSE,skip = 17,sep = "")         #read with skipping gromacs output header
})
names(rmsf_calc) <- basename(list.files(pattern = "rmsf_rep*", recursive = TRUE)) #strip parent folder names

rmsf_arr <- data.frame(cbind(rmsf_calc$rmsf_rep1_arr.xvg$V2,
                             rmsf_calc$rmsf_rep2_arr.xvg$V2,
                             rmsf_calc$rmsf_rep3_arr.xvg$V2,
                             rmsf_calc$rmsf_rep4_arr.xvg$V2,
                             rmsf_calc$rmsf_rep5_arr.xvg$V2), 
                       row.names = rmsf_calc$rmsf_rep1_arr.xvg$V1,stringsAsFactors = FALSE)
arrxval <- as.numeric(rmsf_calc$rmsf_rep1_arr.xvg$V1)
arryval <- matrix(rowMeans(rmsf_arr)*10, c(length(rmsf_calc$rmsf_rep1_arr.xvg$V1),1)) 
arryvalsd <- matrix(apply(rmsf_arr,1,sd)*10, c(length(rmsf_calc$rmsf_rep1_arr.xvg$V1),1))

rmsf_arr_converged <- data.frame(cbind(rmsf_calc$rmsf_rep2_arr.xvg$V2,
                                       rmsf_calc$rmsf_rep3_arr.xvg$V2,
                                       rmsf_calc$rmsf_rep4_arr.xvg$V2,
                                       rmsf_calc$rmsf_rep5_arr.xvg$V2), 
                                 row.names = rmsf_calc$rmsf_rep1_arr.xvg$V1,stringsAsFactors = FALSE)
arryval_converged <- matrix(rowMeans(rmsf_arr_converged)*10, c(length(rmsf_calc$rmsf_rep1_arr.xvg$V1),1)) 

arryval_one <- matrix(rmsf_calc$rmsf_rep1_arr.xvg$V2*10, 
                      c(length(rmsf_calc$rmsf_rep1_arr.xvg$V1),1))




rmsf_ntr <- data.frame(cbind(rmsf_calc$rmsf_rep1_ntr.xvg$V2,
                             rmsf_calc$rmsf_rep2_ntr.xvg$V2,
                             rmsf_calc$rmsf_rep3_ntr.xvg$V2,
                             rmsf_calc$rmsf_rep4_ntr.xvg$V2,
                             rmsf_calc$rmsf_rep5_ntr.xvg$V2),
                       row.names = rmsf_calc$rmsf_rep1_ntr.xvg$V1,stringsAsFactors = FALSE)
ntrxval <- c(50:382)#as.numeric(rmsf_calc$rmsf_rep1_ntr.xvg$V1)
ntryval <- matrix(rowMeans(rmsf_ntr)*10, c(length(rmsf_calc$rmsf_rep1_ntr.xvg$V1),1))
ntryvalsd <- matrix(apply(rmsf_ntr,1,sd)*10, c(length(rmsf_calc$rmsf_rep1_ntr.xvg$V1),1))

####helices####
helices <- read.csv("~/PycharmProjects/untitled/gromacs_tools/6up7_helices.csv",header = FALSE,sep = ",")
sheets <- read.csv("~/PycharmProjects/untitled/gromacs_tools/6up7_sheets.csv",header = FALSE,sep = ",")

helices_arr <- helices[which(as.character(helices$V1) == "B"),]
helices_ntr <- helices[which(as.character(helices$V1) == "R"),]
sheets_arr <- sheets[which(as.character(sheets$V1) == "B"),]
sheets_ntr <- sheets[which(as.character(sheets$V1) == "R"),]

#arrestin elements
helarr <- data.frame(xmin=helices_arr$V2,
                     xmax=helices_arr$V3,
                     ymin=rep(0,length(helices_arr$V1)),
                     ymax=rep(Inf,length(helices_arr$V1)),
                     fill = rep("\u03B1-helix",length(helices_arr$V1)))

shearr <- data.frame(xmin=sheets_arr$V2,
                     xmax=sheets_arr$V3,
                     ymin=rep(0,length(sheets_arr$V1)),
                     ymax=rep(Inf,length(sheets_arr$V1)),
                     fill = rep("\u03B2-sheet",length(sheets_arr$V1)))
#receptor elements
helntr <- data.frame(xmin=helices_ntr$V2,
                     xmax=helices_ntr$V3,
                     ymin=rep(0,length(helices_ntr$V1)),
                     ymax=rep(Inf,length(helices_ntr$V1)),
                     fill = rep("\u03B1-helix",length(helices_ntr$V1)))

shentr <- data.frame(xmin=sheets_ntr$V2,
                     xmax=sheets_ntr$V3,
                     ymin=rep(0,length(sheets_ntr$V1)),
                     ymax=rep(Inf,length(sheets_ntr$V1)),
                     fill = rep("\u03B2-sheet",length(sheets_ntr$V1)))


#arrestin plot
update_geom_defaults("line", list(size=.5))
arrestin_rmsf <- ggplot() +  
  xlab(paste("\u03B2arr-1 position")) + ylab("RMSF [\u00C5]") +   ylim(0,12.5) +
  geom_rect(data=helarr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shearr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_ribbon(aes(x=arrxval,y=arryval,ymin=arryval[,1]-arryvalsd[,1], ymax=arryval[,1]+arryvalsd[,1],group=1),
              fill="grey60", size=.5)+
  geom_line(aes(x=arrxval,y=arryval[,1],group=1)) + scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),
        legend.title = element_text(size = 10), legend.position = "top", 
        legend.text = element_text(size=8)) 
legend_sec <- get_legend(arrestin_rmsf)
arrestin_rmsf <- arrestin_rmsf + theme(legend.position = "none")
legend_sec <- arrangeGrob(legend_sec,top = textGrob("a", 
                                                          x=unit(0,"npc"), 
                                                          y = unit(1,"npc"),
                                                          just=c("left","top")))


# 
# arrestin_plot_four <- ggplot() +  
#   xlab(paste("\u03B2arr-1 position")) + ylab("RMSF [\u00C5]") +   ylim(0,10) +
#   geom_rect(data=helarr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
#   geom_rect(data=shearr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
#   geom_line(aes(x=arrxval,y=arryval_converged[,1],group=1)) + scale_fill_discrete(name = "secondary structure") +
#   theme(axis.text = element_text(size=8), axis.title = element_text(size=10),
#         legend.title = element_text(size = 10), legend.position = "bottom", 
#         legend.text = element_text(size=8) )
# arrestin_plot_four
# ggsave("20201021_RMSFarr_ntr_four.png", plot=arrestin_plot_four, device=png(), 
#        path="~/Desktop/figures/", width=16, height = 8, units = "cm", dpi=300)
# #top="Average RMSF per residue for 5 x 3 µs simulations"
# dev.off()
# 
# arrestin_plot_one <- ggplot() +  
#   xlab(paste("\u03B2arr-1 position")) + ylab("RMSF [\u00C5]") +   ylim(0,10) +
#   geom_rect(data=helarr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
#   geom_rect(data=shearr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
#   geom_line(aes(x=arrxval,y=arryval_one[,1],group=1)) + scale_fill_discrete(name = "secondary structure") +
#   theme(axis.text = element_text(size=8), axis.title = element_text(size=10),
#         legend.title = element_text(size = 10), legend.position = "bottom", 
#         legend.text = element_text(size=8) )
# arrestin_plot_four
# ggsave("20201021_RMSFarr_ntr_one.png", plot=arrestin_plot_one, device=png(), 
#        path="~/Desktop/figures/", width=16, height = 8, units = "cm", dpi=300)
# #top="Average RMSF per residue for 5 x 3 µs simulations"
# dev.off()


#ntr plot
nts1r_rmsf <- ggplot() +  
  xlab("NTS1R position") + ylab("RMSF [\u00C5]") +   ylim(0,12.5) +
  geom_rect(data=helntr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shentr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_ribbon(aes(x=ntrxval,
                  y=c(ntryval[1:224,1],rep(NA,12),ntryval[225:length(ntryval),1],group=1),
                  ymin=c(ntryval[1:224,1],rep(NA,12),ntryval[225:length(ntryval),1],group=1)-
                    c(ntryvalsd[1:224,1],rep(NA,12),ntryvalsd[225:length(ntryvalsd),1],group=1), 
                  ymax=c(ntryval[1:224,1],rep(NA,12),ntryval[225:length(ntryval),1],group=1)+
                    c(ntryvalsd[1:224,1],rep(NA,12),ntryvalsd[225:length(ntryvalsd),1],group=1)), 
              fill="grey60", size=.5)+
  geom_line(aes(x=ntrxval,y=c(ntryval[1:224,1],rep(NA,12),ntryval[225:length(ntryval),1],group=1))) + 
  scale_fill_discrete(name = "secondary\nstructures") + 
  scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) + 
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),legend.position = "none")


rmsf_overview <- grid.arrange(legend_sec,arrangeGrob(arrestin_rmsf,nts1r_rmsf,nrow=2),nrow=2,heights=c(0.1,0.9)) 
ggsave("20201021_RMSFarr_ntr_wholetrj_sd.png", plot=rmsf_overview, device=png(), 
       path="~/Desktop/figures/", width=16, height = 8, units = "cm", dpi=300)
#top="Average RMSF per residue for 5 x 3 µs simulations"
dev.off()

####raw plots####
update_geom_defaults("line", list(size=.5))
nts1r_rmsf_1 <- ggplot() +  
  xlab("") + ylab("") +   ylim(0,12.5) +
  geom_rect(data=helntr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shentr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_line(aes(x=ntrxval,y=c(rmsf_ntr$X1[1:224]*10,rep(NA,12),rmsf_ntr$X1[225:length(rmsf_ntr$X1)]*10,group=1))) + 
  scale_fill_discrete(name = "secondary\nstructures") + 
  scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) + 
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),legend.position = "none")

nts1r_rmsf_2 <- ggplot() +  
  xlab("") + ylab("") +   ylim(0,12.5) +
  geom_rect(data=helntr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shentr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_line(aes(x=ntrxval,y=c(rmsf_ntr$X2[1:224]*10,rep(NA,12),rmsf_ntr$X2[225:length(rmsf_ntr$X2)]*10,group=1))) + 
  scale_fill_discrete(name = "secondary\nstructures") + 
  scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) + 
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),legend.position = "none")

nts1r_rmsf_3 <- ggplot() +  
  xlab("") + ylab("") +   ylim(0,12.5) +
  geom_rect(data=helntr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shentr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_line(aes(x=ntrxval,y=c(rmsf_ntr$X3[1:224]*10,rep(NA,12),rmsf_ntr$X3[225:length(rmsf_ntr$X3)]*10,group=1))) + 
  scale_fill_discrete(name = "secondary\nstructures") + 
  scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) + 
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),legend.position = "none")

nts1r_rmsf_4 <- ggplot() +  
  xlab("") + ylab("") +   ylim(0,12.5) +
  geom_rect(data=helntr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shentr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_line(aes(x=ntrxval,y=c(rmsf_ntr$X4[1:224]*10,rep(NA,12),rmsf_ntr$X4[225:length(rmsf_ntr$X4)]*10,group=1))) + 
  scale_fill_discrete(name = "secondary\nstructures") + 
  scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) + 
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),legend.position = "none")

nts1r_rmsf_5 <- ggplot() +  
  xlab("NTS1R position") + ylab("") +   ylim(0,12.5) +
  geom_rect(data=helntr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shentr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_line(aes(x=ntrxval,y=c(rmsf_ntr$X5[1:224]*10,rep(NA,12),rmsf_ntr$X5[225:length(rmsf_ntr$X5)]*10,group=1))) + 
  scale_fill_discrete(name = "secondary\nstructures") + 
  scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) + 
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),legend.position = "none")


arrestin_rmsf_1 <- ggplot() +  
  xlab(paste("")) + ylab("RMSF [\u00C5] - Rep1") +   ylim(0,12.5) +
  geom_rect(data=helarr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shearr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_line(aes(x=arrxval,y=rmsf_arr$X1*10,group=1)) + scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),
        legend.title = element_text(size = 10), legend.position = "top", 
        legend.text = element_text(size=8))
legend_sec <- get_legend(arrestin_rmsf_1)
arrestin_rmsf_1 <- arrestin_rmsf_1 + theme(legend.position = "none")

arrestin_rmsf_2 <- ggplot() +  
  xlab(paste("")) + ylab("RMSF [\u00C5] - Rep2") +   ylim(0,12.5) +
  geom_rect(data=helarr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shearr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_line(aes(x=arrxval,y=rmsf_arr$X2*10,group=1)) + scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),
        legend.title = element_text(size = 10), legend.position = "none")

arrestin_rmsf_3 <- ggplot() +  
  xlab(paste("")) + ylab("RMSF [\u00C5] - Rep3") +   ylim(0,12.5) +
  geom_rect(data=helarr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shearr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_line(aes(x=arrxval,y=rmsf_arr$X3*10,group=1)) + scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),
        legend.title = element_text(size = 10), legend.position = "none")

arrestin_rmsf_4 <- ggplot() +  
  xlab(paste("")) + ylab("RMSF [\u00C5] - Rep4") +   ylim(0,12.5) +
  geom_rect(data=helarr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shearr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_line(aes(x=arrxval,y=rmsf_arr$X4*10,group=1)) + scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),
        legend.title = element_text(size = 10), legend.position = "none")

arrestin_rmsf_5 <- ggplot() +  
  xlab(paste("\u03B2arr-1 position")) + ylab("RMSF [\u00C5] - Rep5") +   ylim(0,12.5) +
  geom_rect(data=helarr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_rect(data=shearr, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  geom_line(aes(x=arrxval,y=rmsf_arr$X5*10,group=1)) + scale_fill_manual(name = "secondary structure", values=c("tan1","paleturquoise")) +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10),
        legend.title = element_text(size = 10), legend.position = "none")

arrestins <- arrangeGrob(arrestin_rmsf_1,arrestin_rmsf_2,arrestin_rmsf_3,arrestin_rmsf_4,arrestin_rmsf_5,ncol=1)
arrestins <- arrangeGrob(arrestins,top = textGrob("a",
                                                        x=unit(0.08,"npc"), 
                                                        y = unit(1,"npc"),
                                                        just=c("left","top")))
receptors <- arrangeGrob(nts1r_rmsf_1,nts1r_rmsf_2,nts1r_rmsf_3,nts1r_rmsf_4,nts1r_rmsf_5,ncol=1)
receptors <- arrangeGrob(receptors,top = textGrob("b",
                                                  x=unit(0.08,"npc"), 
                                                  y = unit(1,"npc"),
                                                  just=c("left","top")))
rmsf_raw <- grid.arrange(arrangeGrob(arrestins,
                                     receptors,ncol=2),
                         legend_sec,
                         ncol=1,heights=c(0.95,0.05)) 
ggsave("20201021_RMSFarr_ntr_raw.png", plot=rmsf_raw, device=png(), 
       path="~/Desktop/figures/", width=16, height = 20, units = "cm", dpi=300)
#top="Average RMSF per residue for 5 x 3 µs simulations"
dev.off()
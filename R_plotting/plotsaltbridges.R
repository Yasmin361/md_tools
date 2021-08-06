library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(tidyr)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

rep1 <- read.csv("~/Desktop/run002/analyses_complete/VMD_saltbridges/rep1_all_saltbridges.tml_saltbridges_inter.list", sep = "\t", header = FALSE,
                   stringsAsFactors=FALSE)
rep2 <- read.csv("~/Desktop/run002/analyses_complete/VMD_saltbridges/rep2_all_saltbridges.tml_saltbridges_inter.list", sep = "\t", header = FALSE,
                 stringsAsFactors=FALSE)
rep3 <- read.csv("~/Desktop/run002/analyses_complete/VMD_saltbridges/rep3_all_saltbridges.tml_saltbridges_inter.list", sep = "\t", header = FALSE,
                 stringsAsFactors=FALSE)
rep4 <- read.csv("~/Desktop/run002/analyses_complete/VMD_saltbridges/rep4_all_saltbridges.tml_saltbridges_inter.list", sep = "\t", header = FALSE,
                 stringsAsFactors=FALSE)
rep5 <- read.csv("~/Desktop/run002/analyses_complete/VMD_saltbridges/rep5_all_saltbridges.tml_saltbridges_inter.list", sep = "\t", header = FALSE,
                 stringsAsFactors=FALSE)

reps <- list(data.frame(rep1[-c(2,3)]),
             data.frame(rep2[-c(2,3)]),
             data.frame(rep3[-c(2,3)]),
             data.frame(rep4[-c(2,3)]),
             data.frame(rep5[-c(2,3)]))
longr <- list()
plots <- list()

# for (i in c(1:5)) {
#   rownames(reps[i]) <- reps[i]$V1
#   print(length(reps[i]))
#   print(ncol(reps[i]))
#   colnames(reps[i]) <- c("pairnames",seq(1,length(reps[i])-1))
#   arrres <- as.numeric(c(sub("...","",sub("--.*","",reps[i]$pairnames))))
#   ntrres <- as.numeric(c(sub("...","",sub(".*--","",reps[i]$pairnames))))
#   
#   longr[i] <- pivot_longer(reps[i], cols = !pairnames, names_to = "times", values_to = "saltbridge")
#   longr[i] <- cbind(arrres=rep(arrres,each=length(reps[i])-1),longr[i])
#   longr[i] <- cbind(ntrres=rep(ntrres,each=length(reps[i])-1),longr[i])
#   plots[i] <- ggplot(longr[i]) + geom_tile(aes(x=as.numeric(times), 
#                                                y=reorder(pairnames,-arrres), 
#                                                fill=factor(as.numeric(saltbridge)*ntrres))) + 
#     ylab("barr1-NTSR1R pairs") + xlab("Time [ns]") + main(paste("Rep ",i,sep = "")) + 
#     theme(legend.position = "none") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
#     scale_fill_manual(name= "Salt Bridge", values=c("white",viridis(length(ntrres)))) +
#     theme(axis.text=element_text(size=10), axis.title=element_text(size=14)) + 
#     theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
# }

rep1 <- data.frame(rep1)
rownames(rep1) <- rep1$V1
rep1 <- subset(rep1,select = -c(2,3))
colnames(rep1) <- c("pairnames",seq(1,length(rep1)-1))
arrres <- as.numeric(c(sub("...","",sub("--.*","",rep1$pairnames))))
ntrres <- as.numeric(c(sub("...","",sub(".*--","",rep1$pairnames))))

r1 <- pivot_longer(rep1, cols = !pairnames, names_to = "times", values_to = "saltbridge")
r1 <- cbind(arrres=rep(arrres,each=length(rep1)-1),r1)
r1 <- cbind(ntrres=rep(ntrres,each=length(rep1)-1),r1)
plot1 <- ggplot(r1) + geom_tile(aes(x=as.numeric(times), y=reorder(pairnames,-arrres), fill=factor(saltbridge))) + 
  ylab("Rep 1") + xlab("Time [ns]") + theme(legend.position = "none") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  scale_fill_manual(name= "Salt Bridge", values=c("white","red")) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=14)) + 
  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))


rep2 <- data.frame(rep2)
rownames(rep2) <- rep2$V1
rep2 <- subset(rep2,select = -c(2,3))
colnames(rep2) <- c("pairnames",seq(1,length(rep2)-1))
arrres <- as.numeric(c(sub("...","",sub("--.*","",rep2$pairnames))))
ntrres <- as.numeric(c(sub("...","",sub(".*--","",rep2$pairnames))))

r2 <- pivot_longer(rep2, cols = !pairnames, names_to = "times", values_to = "saltbridge")
r2 <- cbind(arrres=rep(arrres,each=length(rep2)-1),r2)
r2 <- cbind(ntrres=rep(ntrres,each=length(rep2)-1),r2)
plot2 <- ggplot(r2) + geom_tile(aes(x=as.numeric(times), y=reorder(pairnames,-arrres), fill=factor(saltbridge))) + 
  ylab("Rep 2") + xlab("Time [ns]") + theme(legend.position = "none") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  scale_fill_manual(name= "Salt Bridge", values=c("white","red")) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=14)) + 
  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))


rep3 <- data.frame(rep3)
rownames(rep3) <- rep3$V1
rep3 <- subset(rep3,select = -c(2,3))
colnames(rep3) <- c("pairnames",seq(1,length(rep3)-1))
arrres <- as.numeric(c(sub("...","",sub("--.*","",rep3$pairnames))))
ntrres <- as.numeric(c(sub("...","",sub(".*--","",rep3$pairnames))))

r3 <- pivot_longer(rep3, cols = !pairnames, names_to = "times", values_to = "saltbridge")
r3 <- cbind(arrres=rep(arrres,each=length(rep3)-1),r3)
r3 <- cbind(ntrres=rep(ntrres,each=length(rep3)-1),r3)
plot3 <- ggplot(r3) + geom_tile(aes(x=as.numeric(times), y=reorder(pairnames,-arrres),  fill=factor(saltbridge))) + 
  ylab("Rep 3") + xlab("Time [ns]") + theme(legend.position = "none") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  scale_fill_manual(name= "Salt Bridge", values=c("white","red")) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=14)) + 
  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))


rep4 <- data.frame(rep4)
rownames(rep4) <- rep4$V1
rep4 <- subset(rep4,select = -c(2,3))
colnames(rep4) <- c("pairnames",seq(1,length(rep4)-1))
arrres <- as.numeric(c(sub("...","",sub("--.*","",rep4$pairnames))))
ntrres <- as.numeric(c(sub("...","",sub(".*--","",rep4$pairnames))))

r4 <- pivot_longer(rep4, cols = !pairnames, names_to = "times", values_to = "saltbridge")
r4 <- cbind(arrres=rep(arrres,each=length(rep4)-1),r4)
r4 <- cbind(ntrres=rep(ntrres,each=length(rep4)-1),r4)
plot4 <- ggplot(r4) + geom_tile(aes(x=as.numeric(times), y=reorder(pairnames,-arrres),  fill=factor(saltbridge))) + 
  ylab("Rep 4") + xlab("Time [ns]") + theme(legend.position = "none") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  scale_fill_manual(name= "Salt Bridge", values=c("white","red")) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=14)) + 
  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))


rep5 <- data.frame(rep5)
rownames(rep5) <- rep5$V1
rep5 <- subset(rep5,select = -c(2,3))
colnames(rep5) <- c("pairnames",seq(1,length(rep5)-1))
arrres <- as.numeric(c(sub("...","",sub("--.*","",rep5$pairnames))))
ntrres <- as.numeric(c(sub("...","",sub(".*--","",rep5$pairnames))))

r5 <- pivot_longer(rep5, cols = !pairnames, names_to = "times", values_to = "saltbridge")
r5 <- cbind(arrres=rep(arrres,each=length(rep5)-1),r5)
r5 <- cbind(ntrres=rep(ntrres,each=length(rep5)-1),r5)
plot5 <- ggplot(r5) + geom_tile(aes(x=as.numeric(times), y=reorder(pairnames,-arrres), fill=factor(saltbridge))) + 
  ylab("Rep 5") + xlab("Time [ns]") + scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  scale_fill_manual(name= "Salt Bridge", values=c("white","red")) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=14)) + 
  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))

legend1 <- get_legend(plot5)
plot5 <- plot5  + theme(legend.position = "none")

# lay <- rbind(c(1,1,6),
#              c(2,2,6),
#              c(3,3,6),
#              c(4,4,6),
#              c(5,5,6))

overview_salt <- grid.arrange(plot1,plot2,plot3,plot4,plot5,legend1,ncol=6,
                              top="Intermolecular Salt Bridges")
overview_salt


lysglu <- read.csv("~/Desktop/run002/analyses_complete/VMD_saltbridges/LYS77_GLU313.list", sep = "", header = FALSE,
                   stringsAsFactors=FALSE, na = "0")

lysglu <- data.frame(lysglu)
lysglu[is.na(lysglu)] <- 0
rownames(lysglu) <- c("Rep1","Rep2","Rep3","Rep4","Rep5")
lysglu$V1 <- c("Rep1","Rep2","Rep3","Rep4","Rep5")
lysglu <- subset(lysglu,select = -c(2,3))
colnames(lysglu) <- c("pairnames",seq(1,length(lysglu)-1))
lysg <- pivot_longer(lysglu, cols = !pairnames, names_to = "times", values_to = "saltbridge")
lgplot <- ggplot(lysg) + geom_tile(aes(x=as.numeric(times), y=rev(pairnames), fill=factor(saltbridge))) + 
  ylab("\u03B2arr-1 K77-E313") + xlab("Time [ns]") + 
  scale_y_discrete(labels=rev(c("Rep1","Rep2","Rep3","Rep4","Rep5"))) +
  scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +  
  scale_fill_manual(name= "Salt Bridge", values=c("white","red")) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=14))+
  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
lgplot



####old plots 
# rep2 <- data.frame(rep2)
# rownames(rep2) <- rep2$V1
# rep2 <- subset(rep2,select = -c(2,3))
# colnames(rep2) <- c("pairnames",seq(1,length(rep2)-1))
# r2 <- pivot_longer(rep2, cols = !pairnames, names_to = "times", values_to = "saltbridge")
# plot2 <- ggplot(r2) + geom_tile(aes(x=as.numeric(times), y=pairnames, fill=factor(saltbridge))) + 
#   ylab("Rep 2") + xlab("") + theme(legend.position = "none") + 
#   scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
#   scale_fill_manual(name= "Salt Bridge", values=c("white","red"))+
#   theme(axis.text=element_text(size=10), axis.title=element_text(size=14)) +
#   theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
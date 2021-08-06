library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(tidyr)
library(tibble)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

rep1 <- read.table("~/Desktop/run002/analyses_complete/MD-task_crosscorr/crosscorrelation_rep1_500ns.txt", 
                   quote="\"", comment.char="", stringsAsFactors=FALSE)
residues <- read.table("~/Desktop/run002/analyses_complete/MD-task_crosscorr/names.list", 
                       quote="\"", sep = "\t", comment.char="", stringsAsFactors=FALSE)
rownames(residues) <- c(1:(nrow(residues)))
colnames(rep1) <- row.names(residues)
rownames(rep1) <- row.names(residues)
rep1_nol <- rep1[-seq(349,354),-seq(349,354)]
rep1_nol <- data.frame(melt(rep1_nol))
reslabels <- residues$V2[c(seq(63,343,70),seq(429,639,70))]
rep1_nol["res"] <- rep(1:668,668)


cbind(rep1[471],residues$V2)

rep2 <- read.table("~/Desktop/run002/analyses_complete/MD-task_crosscorr/crosscorrelation_rep2_500ns.txt", 
                   quote="\"", comment.char="", stringsAsFactors=FALSE)
colnames(rep2) <- row.names(residues)
rownames(rep2) <- row.names(residues)
rep2_nol <- rep2[-seq(349,354),-seq(349,354)]
rep2_nol <- data.frame(melt(rep2_nol))
rep2_nol["res"] <- rep(1:668,668)

rep3 <- read.table("~/Desktop/run002/analyses_complete/MD-task_crosscorr/crosscorrelation_rep3_500ns.txt", 
                   quote="\"", comment.char="", stringsAsFactors=FALSE)
colnames(rep3) <- row.names(residues)
rownames(rep3) <- row.names(residues)
rep3_nol <- rep3[-seq(349,354),-seq(349,354)]
rep3_nol <- data.frame(melt(rep3_nol))
rep3_nol["res"] <- rep(1:668,668)

rep4 <- read.table("~/Desktop/run002/analyses_complete/MD-task_crosscorr/crosscorrelation_rep4_500ns.txt", 
                   quote="\"", comment.char="", stringsAsFactors=FALSE)
colnames(rep4) <- row.names(residues)
rownames(rep4) <- row.names(residues)
rep4_nol <- rep4[-seq(349,354),-seq(349,354)]
rep4_nol <- data.frame(melt(rep4_nol))
rep4_nol["res"] <- rep(1:668,668)


rep5 <- read.table("~/Desktop/run002/analyses_complete/MD-task_crosscorr/crosscorrelation_rep5_500ns.txt", 
                   quote="\"", comment.char="", stringsAsFactors=FALSE)
colnames(rep5) <- row.names(residues)
rownames(rep5) <- row.names(residues)
rep5_nol <- rep5[-seq(349,354),-seq(349,354)]
rep5_nol <- data.frame(melt(rep5_nol))
rep5_nol["res"] <- rep(1:668,668)

####plots####
update_geom_defaults("line", list(size=.5))
plot1 <- ggplot(rep1_nol) + geom_tile(aes(x=variable, y=res, fill=value)) +
  scale_fill_gradient2(low ="cyan",mid = "white",high = "red",midpoint = 0,labels=c(-1,-0.5,0,0.5,1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_hline(yintercept =  348, color="black") +
  geom_vline(xintercept =  348, color="black") +
  geom_hline(yintercept =  572, color="black", linetype="dashed") +
  geom_vline(xintercept =  572, color="black", linetype="dashed") +
  scale_x_discrete(breaks=c(seq(63,343,70),seq(429,639,70)),labels=reslabels) +
  scale_y_discrete(breaks=c(seq(63,343,70),seq(429,639,70)),labels=reslabels) + 
  theme(legend.text = element_text(size=8), legend.title = element_text(size=10))+
  labs(fill="correlation coefficient") +
  xlab("Rep 1") + ylab("")
leg <- get_legend(plot1)
plot1 <- plot1 + coord_fixed() + theme(legend.position = "none") + theme(axis.text=element_text(size=8),
                                                                         axis.title=element_text(size=10))


plot2 <- ggplot(rep2_nol) + geom_tile(aes(x=variable, y=res, fill=value)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_gradient2(low ="cyan",mid = "white",high = "red",midpoint = 0) +
  geom_hline(yintercept =  348, color="black") +
  geom_vline(xintercept =  348, color="black") +
  geom_hline(yintercept =  572, color="black", linetype="dashed") +
  geom_vline(xintercept =  572, color="black", linetype="dashed") +
  scale_x_discrete(breaks=c(seq(63,343,70),seq(429,639,70)),labels=reslabels) +
  scale_y_discrete(breaks=c(seq(63,343,70),seq(429,639,70)),labels=reslabels) + 
  xlab("Rep 2") + ylab("")
plot2 <- plot2 + coord_fixed() + theme(legend.position = "none") + theme(axis.text=element_text(size=8),
                                                                         axis.title=element_text(size=10))

plot3 <- ggplot(rep3_nol) + geom_tile(aes(x=variable, y=res, fill=value)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_gradient2(low ="cyan",mid = "white",high = "red",midpoint = 0) +
  geom_hline(yintercept =  348, color="black") +
  geom_vline(xintercept =  348, color="black") +
  geom_hline(yintercept =  572, color="black", linetype="dashed") +
  geom_vline(xintercept =  572, color="black", linetype="dashed") +
  scale_x_discrete(breaks=c(seq(63,343,70),seq(429,639,70)),labels=reslabels) +
  scale_y_discrete(breaks=c(seq(63,343,70),seq(429,639,70)),labels=reslabels) + 
  xlab("Rep 3")  +  ylab("")
plot3 <- plot3 + coord_fixed() + theme(legend.position = "none") + theme(axis.text=element_text(size=8),
                                                                          axis.title=element_text(size=10))

plot4 <- ggplot(rep4_nol) + geom_tile(aes(x=variable, y=res, fill=value)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_gradient2(low ="cyan",mid = "white",high = "red",midpoint = 0) +
  geom_hline(yintercept =  348, color="black") +
  geom_vline(xintercept =  348, color="black") +
  geom_hline(yintercept =  572, color="black", linetype="dashed") +
  geom_vline(xintercept =  572, color="black", linetype="dashed") +
  scale_x_discrete(breaks=c(seq(63,343,70),seq(429,639,70)),labels=reslabels) +
  scale_y_discrete(breaks=c(seq(63,343,70),seq(429,639,70)),labels=reslabels) + 
  xlab("Rep 4") + ylab("ßarr-1 and NTS1R residues")
plot4 <- plot4 + coord_fixed() + theme(legend.position = "none") + theme(axis.text=element_text(size=8),
                                                                         axis.title=element_text(size=10))

plot5 <- ggplot(rep5_nol) + geom_tile(aes(x=variable, y=res, fill=value)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_gradient2(low ="cyan",mid = "white",high = "red",midpoint = 0) +
  geom_hline(yintercept =  348, color="black") +
  geom_vline(xintercept =  348, color="black") +
  geom_hline(yintercept =  572, color="black", linetype="dashed") +
  geom_vline(xintercept =  572, color="black", linetype="dashed") +
  scale_x_discrete(breaks=c(seq(63,343,70),seq(429,639,70)),labels=reslabels) +
  scale_y_discrete(breaks=c(seq(63,343,70),seq(429,639,70)),labels=reslabels) + 
  xlab("Rep 5") + labs(fill = "Correlation")  +  ylab("")
plot5 <- plot5 + coord_fixed() + theme(legend.position = "none") + theme(axis.text=element_text(size=8),
                                                                         axis.title=element_text(size=10))

overview_cc <- grid.arrange(plot1,plot2,plot3,plot4,plot5,leg, nrow=3, ncol=2)
overview_cc

ggsave("20210128_crosscorrelation_new.png", plot=overview_cc, device=png(), 
       path="~/Desktop/figures/", width=16, height = 20, units = "cm", dpi=300)
dev.off()

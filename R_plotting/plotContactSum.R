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
setwd("~/Desktop/run002/analyses_complete/mdciao")        #supply path here
arr_contacts <- read.csv("sum_of_replicas_arr_contacts.list", header = FALSE,sep = "\t", stringsAsFactors=FALSE)
ntr_contacts <- read.csv("sum_of_replicas_ntr_contacts.list", header = FALSE,sep = "\t", stringsAsFactors=FALSE)
pair_contacts <- read.csv("sum_of_replicas_pairs.list", header = FALSE,sep = "\t", stringsAsFactors=FALSE)

arr <- data.frame(description = arr_contacts$V1,
                  resseq =  as.numeric(sub("[A-Z]","", sub("@.*", "", arr_contacts$V1))),
                  freq=arr_contacts$V2)
arr_plot <- ggplot(arr, aes(x=reorder(description,resseq), y=freq)) +
  geom_bar(stat = "identity") + ylab("Cumulative contact frequency for all simulations") +
  xlab("barr1 residues") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
arr_plot

ntr <- data.frame(description = c(ntr_contacts$V1),
                  freq=c(ntr_contacts$V2))
ntr_plot <- ggplot(ntr, aes(x=reorder(description, -freq), y=freq)) +
  geom_bar(stat = "identity") + ylab("Cumulative contact frequency for all simulations") +
  xlab("NTS1R residues") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ntr_plot

pairs <- data.frame(description = as.character(pair_contacts$V1),
                    shortdescription = sub("@.*","",sub("@frag0","", pair_contacts$V1)),
                    arrseq = as.numeric(sub("[A-Z]","", sub("@.*", "", pair_contacts$V1))),
                    ntrseq = as.numeric(sub("[A-Z]","",sub("@.*","",sub(".*-", "", pair_contacts$V1)))),
                    freq=as.numeric(pair_contacts$V2))

####plots####

####consurf####
setwd("~/Desktop/run002/analyses_complete/consurf")        #supply path here
ntrconsurf <- read.csv("6up7-R/consurf.grades", header = FALSE,skip=15,sep = "\t", stringsAsFactors=FALSE)
ntrconsurf$V3 <- sub("\\s*-",NA,ntrconsurf$V3)
ntrconsurf <- ntrconsurf[complete.cases(ntrconsurf$V3),]
ntrconsurf$V3 <- as.numeric(sub("[A-Z][A-Z][A-Z]","",sub(":.*","",ntrconsurf$V3)))
ntrconsurf$V6 <- as.numeric(sub("[\\*]","",ntrconsurf$V6))
ntrconsurf <- ntrconsurf[complete.cases(ntrconsurf$V6),]
ntrconservation <- c()
for (i in c(1:length(pairs$ntrseq))) {
  ntrconservation[i] <- ntrconsurf[which(ntrconsurf$V3 == pairs[order(pairs$ntrseq),]$ntrseq[i]), ]$V6
}

setwd("~/Desktop/run002/analyses_complete/consurf")        #supply path here
arrconsurf <- read.csv("6up7-B/consurf.grades", header = FALSE,skip=15,sep = "\t", stringsAsFactors=FALSE)
arrconsurf$V3 <- sub("\\s*-",NA,arrconsurf$V3)
arrconsurf <- arrconsurf[complete.cases(arrconsurf$V3),]
arrconsurf$V3 <- as.numeric(sub("[A-Z][A-Z][A-Z]","",sub(":.*","",arrconsurf$V3)))
arrconsurf$V6 <- as.numeric(sub("[\\*]","",arrconsurf$V6))
arrconsurf <- arrconsurf[complete.cases(arrconsurf$V6),]
arrconservation <- c()
for (i in c(1:length(pairs$ntrseq))) {
  arrconservation[i] <- arrconsurf[which(arrconsurf$V3 == pairs[order(pairs$arrseq),]$arrseq[i]), ]$V6
} 



pairs <- cbind(pairs[order(pairs$ntrseq),],ntrconservation,arrconservation)



conservation_palette <- c(brewer.pal(9, "RdBu"))

# order dataframe by one or more columns
pairs <- pairs[order(pairs[,4]),]

####plots####

pair_plot_finger <- ggplot(subset(pairs,subset = arrseq < 79), aes(x=reorder(shortdescription, arrseq), y=freq)) +
  geom_bar(stat = "identity", aes(fill=ntrseq)) + ylab("CCF") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) + labs(fill="NTS1R residues") +
  xlab(paste("\u03B2arr-1 finger loop (\u03A3 ", sum(subset(pairs,subset = arrseq < 79)$freq),")",sep="")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_viridis_c() + theme(legend.position = "none")  # theme(legend.position=c(0.25,0.9),legend.direction = "horizontal")
pair_plot_finger

ggsave("20201201_cumulative_pairsums_fl.png", plot=pair_plot_finger, device=png(), 
       path="~/Desktop/figures/", width=12, height = 4, units = "cm", dpi=300)
dev.off()

pair_plot_finger_cons <- ggplot(subset(pairs,subset = arrseq < 79), aes(x=reorder(shortdescription, arrseq), y=freq)) +
  geom_bar(stat = "identity", aes(fill=ntrconservation)) + ylab("CCF") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) + labs(fill="NTS1R residues") +
  xlab("Intermolecular contacts of finger loop") + 
  scale_fill_viridis_c(option = "C") + theme(legend.position = "none")  # theme(legend.position=c(0.25,0.9),legend.direction = "horizontal")
pair_plot_finger_cons



pair_plot_other <- ggplot(subset(pairs,subset = arrseq > 79), aes(x=reorder(shortdescription, arrseq), y=freq)) +
  geom_bar(stat = "identity", aes(fill=ntrseq)) + ylab("CCF") + 
  scale_fill_viridis_c() + theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) + 
  labs(fill="NTS1R residues") + ylim(c(0,10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(legend.position = "none") +# theme(legend.position=c(0.25,0.8),legend.direction = "horizontal") + 
  geom_vline(xintercept = 7.5, linetype="dashed", color = "black") +
  geom_vline(xintercept = 11.5, linetype="dashed", color = "black") +
  geom_vline(xintercept = 25.5, linetype="dashed", color = "black") +
  geom_vline(xintercept = 29.5, linetype="dashed", color = "black") +
  xlab(paste("other \u03B2arr-1 loops (\u03A3 ", sum(subset(pairs,subset = arrseq > 78)$freq),")",sep = "")) +
  annotation_custom(textGrob(paste("middle loop\n(\u03A3 ", sum(subset(pairs,subset = arrseq > 126 & arrseq < 144)$freq),")",sep=""),x=0.02,y=0.85,hjust=0,gp=gpar(fontsize=8))) + 
  annotation_custom(textGrob(paste("160-loop\n(\u03A3 ", sum(subset(pairs,subset = arrseq > 150 & arrseq < 162)$freq),")",sep=""),x=0.22,y=0.85,hjust=0,gp=gpar(fontsize=8))) + 
  annotation_custom(textGrob(paste("C-loop\n(\u03A3 ", sum(subset(pairs,subset = arrseq > 239 & arrseq < 250)$freq),")",sep=""),x=0.34,y=0.85,hjust=0,gp=gpar(fontsize=8))) +
  annotation_custom(textGrob(paste("lariat loop\n(\u03A3 ", sum(subset(pairs,subset = arrseq > 270 & arrseq < 300)$freq),")",sep=""),x=0.74,y=0.85,hjust=0,gp=gpar(fontsize=8))) +
  annotation_custom(textGrob(paste("back loop\n(\u03A3 ", sum(subset(pairs,subset = arrseq > 302 & arrseq < 318)$freq),")",sep = ""),x=0.865,y=0.85,hjust=0,gp=gpar(fontsize=8))) 
#  xlab(paste("middle loop (\u03A3 ", sum(subset(pairs,subset = arrseq > 126 & arrseq < 144)$freq),")",
#             "160-loop (\u03A3 ", sum(subset(pairs,subset = arrseq > 150 & arrseq < 162)$freq),")",
#             "C-loop (\u03A3 ", sum(subset(pairs,subset = arrseq > 239 & arrseq < 250)$freq),")",
#             "lariat loop (\u03A3 ", sum(subset(pairs,subset = arrseq > 270 & arrseq < 300)$freq),")",
#             "back loop (\u03A3 ", sum(subset(pairs,subset = arrseq > 302 & arrseq < 318)$freq),")",sep = "")
#  )
pair_plot_other

ggsave("20201201_cumulative_pairsums_other.png", plot=pair_plot_other, device=png(), 
       path="~/Desktop/figures/", width=16, height = 4, units = "cm", dpi=300)
dev.off()

pair_plot_other_cons <- ggplot(subset(pairs,subset = arrseq > 79), aes(x=reorder(shortdescription, arrseq), y=freq)) +
  geom_bar(stat = "identity", aes(fill=ntrconservation)) + ylab("Cumulative contact frequency for all simulations") + 
  scale_fill_viridis_c(option = "C") + theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) + labs(fill="NTS1R residues") +
  xlab("Intermolecular contacts of arrestin loops") +  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(legend.position = "none") +# theme(legend.position=c(0.25,0.8),legend.direction = "horizontal") + 
  geom_vline(xintercept = 7.5, linetype="dotted", color = "black") +
  geom_vline(xintercept = 11.5, linetype="dotted", color = "black") +
  geom_vline(xintercept = 25.5, linetype="dotted", color = "black") +
  geom_vline(xintercept = 29.5, linetype="dotted", color = "black") +
  annotation_custom(textGrob("middle loop",x=0.05,y=0.95,hjust=0,gp=gpar(fontsize=10))) + 
  annotation_custom(textGrob("160-loop",x=0.23,y=0.95,hjust=0,gp=gpar(fontsize=10))) + 
  annotation_custom(textGrob("C-loop",x=0.5,y=0.95,hjust=0,gp=gpar(fontsize=10))) +
  annotation_custom(textGrob("lariat loop",x=0.75,y=0.95,hjust=0,gp=gpar(fontsize=10))) +
  annotation_custom(textGrob("back loop",x=0.9,y=0.95,hjust=0,gp=gpar(fontsize=10))) 
pair_plot_other_cons

# pair_plot_subset <- ggplot(subset(pairs,freq>=1), aes(x=reorder(residues, -freq), y=freq)) +
#   geom_bar(stat = "identity") + ylab("Cumulative contact frequency for all simulations") +
#   xlab("Intermolecular contact") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# pair_plot_subset


####plots per NTS1R element#####
limitres <- c(81,91,92,98,99,112,166,172,173,182,183,193,248,270,271,294,295,318,361,367,368,382)
lengthsele <- c()
for (i in c(1:11)) {
  lengthsele[i] <- limitres[2*i] - limitres[2*i-1] + 1
}
elements <- c("TM1","ICL1","TM2","TM3","ICL2","TM4","TM5","ICL3","TM6","TM7","H8")

for (i in c(1:11)) {
  plotname <- paste("plot",i,sep="")
  plot <- ggplot(subset(pairs,subset = ntrseq <= limitres[2*i] & ntrseq >= limitres[2*i-1]), 
                 aes(x=reorder(shortdescription, arrseq), y=freq)) +
    geom_bar(stat = "identity", aes(fill=arrseq)) + ylab("CCF") +
    theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) + labs(fill="ßarr-1 residues") +
    xlab(elements[i]) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylim(0,8) +
    scale_fill_viridis_c(limits=c(min(pairs$arrseq),max(pairs$arrseq)))  # theme(legend.position=c(0.25,0.9),legend.direction = "horizontal")
  arrlegend <- get_legend(plot)
  plot <- plot + theme(legend.position = "none")
  assign(plotname,plot)
}
viewplots <- grid.arrange(plot1,plot2,plot4,plot5,plot7,plot8,plot9,plot10,plot11,
                          ncol=3,nrow=3)

# top=textGrob("Intermolecular contacts of \u03B2arr-1 residues with NTS1R structural elements",gp=gpar(fontsize=14)))
viewplots

####plots per NTS1R element cons#####
limitres <- c(81,91,92,98,99,112,166,172,173,182,183,193,248,270,271,294,295,318,361,367,368,382)

elements <- c("TM1","ICL1","TM2","TM3","ICL2","TM4","TM5","ICL3","TM6","TM7","H8")

for (i in c(1:11)) {
  plotname <- paste("plot",i,sep="")
  plot <- ggplot(subset(pairs,subset = ntrseq <= limitres[2*i] & ntrseq >= limitres[2*i-1]), 
                 aes(x=reorder(shortdescription, arrseq), y=freq)) +
    geom_bar(stat = "identity", aes(fill=ntrconservation)) + ylab("CCF") +
    theme(axis.text=element_text(size=8), axis.title=element_text(size=10)) + labs(fill="ßarr-1 residues") +
    xlab(paste(elements[i]," (\u03A3 ",sum(subset(pairs,subset = ntrseq <= limitres[2*i] & ntrseq >= limitres[2*i-1])$freq),")",sep="")) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylim(0,8.5) +
    scale_fill_viridis_c(option = "C") +  # theme(legend.position=c(0.25,0.9),legend.direction = "horizontal")
    theme(legend.position = "none")
  assign(plotname,plot)
}

df <- data.frame(
  dfx=as.factor(c(1:9)),
  dfy=rep(1,9),
  colord=c(1:9)
)

plasma <- ggplot(df) + geom_tile(aes(x=dfx,y=dfy,fill=colord)) + 
  scale_fill_viridis_c(option="C") + theme(legend.position = "none") +
  ylab("") + xlab("ConSurf Score of\ncontact in NTS1R") + theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank()) +
  coord_equal() + theme(axis.title = element_text(size=10), axis.text = element_text(size=8)) +
  annotation_custom(textGrob("variable",x=0.05,y=0.5,hjust=0,gp=gpar(fontsize=8,col="white"))) +
  annotation_custom(textGrob("conserved",x=0.65,y=0.5,hjust=0,gp=gpar(fontsize=8)))
  plasma



viewplots_consurf <- grid.arrange(
  arrangeGrob(plot1,plot2,plot4,nrow = 1,widths = c(5,18,17)),
  arrangeGrob(plot5,plot7,plot8,nrow=1, widths= c(21,9,25)),
  arrangeGrob(plot9,plot10,plot11,plasma, nrow=1,widths=c(6,3,9,4)),
  nrow=3)

# top=textGrob("Intermolecular contacts of \u03B2arr-1 residues with NTS1R structural elements",gp=gpar(fontsize=14)))
viewplots_consurf


ggsave("20210201_cfc_by_receptor_consurf.png", plot=viewplots_consurf, device=png(), 
       path="~/Desktop/figures/", width=24.7, height = 14, units = "cm", dpi=300)
dev.off()



####ligand####
setwd("~/Desktop/run002/analyses_complete/ligand_mdc")        #supply path here
ligand_contact <- read.csv("sum_of_replicas_pairs.list", header = FALSE,sep = "\t", stringsAsFactors=FALSE)
ligand <- data.frame(description = as.character(ligand_contact$V1),
                    shortdescription = sub("@.*","",sub("@frag0","", ligand_contact$V1)),
                    ligseq = as.numeric(sub("[A-Z]","", sub("@.*", "", ligand_contact$V1))),
                    ntrseq = as.numeric(sub("[A-Z]","",sub("@.*","",sub(".*-", "", ligand_contact$V1)))),
                    freq=as.numeric(ligand_contact$V2))

pair_plot_ligand <- ggplot(ligand, aes(x=reorder(shortdescription, ligseq), y=freq)) +
  geom_bar(stat = "identity", aes(fill=ntrseq)) + ylab("Cumulative contact frequency for all simulations") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + labs(fill="NTS1R residues") +
  xlab("Intermolecular contacts of NTS(8-13)") + 
  scale_fill_viridis_c() + theme(legend.position=c(0.1,0.8))
pair_plot_ligand

limitres2 <- c(50,58,59,91,99,130,131,136,137,172,183,207,208,228,229,270,295,329,330,334,335,367)
lengthsele2 <- c()
for (i in c(1:11)) {
  lengthsele2[i] <- limitres2[2*i] - limitres2[2*i-1] + 1
}
elements <- c("NTER","TM1","TM2","ECL1","TM3","TM4","ECL2","TM5","TM6","ECL3","TM7")

for (i in c(1:11)) {
  plotname <- paste("plotl",i,sep="")
  plot <- ggplot(subset(ligand,subset = ntrseq <= limitres2[2*i] & ntrseq >= limitres2[2*i-1]),
                 aes(x=reorder(shortdescription, ligseq), y=freq)) +
    geom_bar(stat = "identity", aes(fill=ligseq)) + ylab("Cumulative contact frequency") +
    theme(axis.text=element_text(size=10), axis.title=element_text(size=14)) + labs(fill="ligand residues") +
    xlab(paste(elements[i]," (\u03A3 ",sum(subset(ligand,subset = ntrseq <= limitres2[2*i] & ntrseq >= limitres2[2*i-1])$freq),")",sep="")) + ylim(0,22.5) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_fill_viridis_c(limits=c(min(ligand$ligseq),max(ligand$ligseq)))  # theme(legend.position=c(0.25,0.9),legend.direction = "horizontal")
  liglegend <- get_legend(plot)
  plot <- plot + theme(legend.position = "none")
  assign(plotname,plot)
}
ligplots <- grid.arrange(plotl1,plotl3,plotl4,plotl5,plotl6,plotl7,plotl8,plotl9,plotl10,plotl11,
                          ncol=5,nrow=2,
                          top=textGrob("Intermolecular contacts of NTS(8-18) residues with NTS1R structural elements",
                                       gp=gpar(fontsize=14)))
ligplots
plotl3

#####dump#####

# arr_plot <- ggplot(arr, aes(x=reorder(residues, -freq), y=freq)) +
#   geom_bar(stat = "identity") + ylab("Cumulative contact frequency for all simulations") +
#   xlab("barr1 residues") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# arr_plot
# 
# ntr <- data.frame(residues = c(ntr_contacts$V1),
#                   freq=c(ntr_contacts$V2))
# ntr_plot <- ggplot(ntr, aes(x=reorder(residues, -freq), y=freq)) +
#   geom_bar(stat = "identity") + ylab("Cumulative contact frequency for all simulations") +
#   xlab("NTS1R residues") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ntr_plot
# 
# pairs <- data.frame(residues = c(pair_contacts$V1),
#                   freq=c(pair_contacts$V2))
# pair_plot <- ggplot(pairs, aes(x=reorder(residues, -freq), y=freq)) +
#   geom_bar(stat = "identity") + ylab("Cumulative contact frequency for all simulations") +
#   xlab("Intermolecular contact") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# pair_plot
# 
# pair_plot_subset <- ggplot(subset(pairs,freq>=1), aes(x=reorder(residues, -freq), y=freq)) +
#   geom_bar(stat = "identity") + ylab("Cumulative contact frequency for all simulations") +
#   xlab("Intermolecular contact") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# pair_plot_subset
#c(rep("finger loop",14),rep("middle loop",5),rep("C-loop",13),rep("lariat loop",2),rep("back loop",3))

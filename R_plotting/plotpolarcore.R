library(ggplot2)
library(gridExtra)
library(grid)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

get_cfc<-function(contacttable){
  trjlength <- length(contacttable)
  ctcs <- length(contacttable[which(contacttable <=4)])
  return(ctcs/trjlength)
}

#####bulk-reading####
setwd("~/Desktop/run002/analyses_complete/mdc_polarcore")        #supply path here
temp <- list.files(pattern = "*ARG169@βarr1\\.step7_production_ext_rep._complete_prot_fit\\.dat", recursive = TRUE) #list files in subfolders
neighborhood <- lapply(temp,function(i){
  read.csv(i,header = FALSE,sep = "",stringsAsFactors = FALSE)         
})

for (i in c(1:5)) {
  neighborhood[[i]][1,1:(ncol(neighborhood[[i]])-1)] <- neighborhood[[i]][1,2:ncol(neighborhood[[i]])]
  neighborhood[[i]] <- neighborhood[[i]][,1:(ncol(neighborhood[[i]])-1)]
  names(neighborhood[[i]]) <- neighborhood[[i]][1,]
  neighborhood[[i]] <- neighborhood[[i]][2:nrow(neighborhood[[i]]),]
  neighborhood[[i]] <- data.frame(neighborhood[[i]])
  neighborhood[[i]] <- lapply(neighborhood[[i]], function(x) as.numeric(x))
}

names(neighborhood) <- basename(list.files(pattern = "*ARG169@βarr1\\.step7_production_ext_rep._complete_prot_fit\\.dat", recursive = TRUE)) #strip parent folder names

#D26
d26_1 <- get_cfc(neighborhood$`rep1.ARG169@βarr1.step7_production_ext_rep1_complete_prot_fit.dat`$D26.βarr1.R169.βarr1.Ang)
d26_2 <- get_cfc(neighborhood$`rep2.ARG169@βarr1.step7_production_ext_rep2_complete_prot_fit.dat`$D26.βarr1.R169.βarr1.Ang)
d26_3 <- get_cfc(neighborhood$`rep3.ARG169@βarr1.step7_production_ext_rep3_complete_prot_fit.dat`$D26.βarr1.R169.βarr1.Ang)
d26_4 <- get_cfc(neighborhood$`rep4.ARG169@βarr1.step7_production_ext_rep4_complete_prot_fit.dat`$D26.βarr1.R169.βarr1.Ang)
d26_5 <- get_cfc(neighborhood$`rep5.ARG169@βarr1.step7_production_ext_rep5_complete_prot_fit.dat`$D26.βarr1.R169.βarr1.Ang)
#D290
d290_1 <- get_cfc(neighborhood$`rep1.ARG169@βarr1.step7_production_ext_rep1_complete_prot_fit.dat`$R169.βarr1.D290.βarr1.Ang)
d290_2 <- get_cfc(neighborhood$`rep2.ARG169@βarr1.step7_production_ext_rep2_complete_prot_fit.dat`$R169.βarr1.D290.βarr1.Ang)
d290_3 <- get_cfc(neighborhood$`rep3.ARG169@βarr1.step7_production_ext_rep3_complete_prot_fit.dat`$R169.βarr1.D290.βarr1.Ang)
d290_4 <- get_cfc(neighborhood$`rep4.ARG169@βarr1.step7_production_ext_rep4_complete_prot_fit.dat`$R169.βarr1.D290.βarr1.Ang)
d290_5 <- get_cfc(neighborhood$`rep5.ARG169@βarr1.step7_production_ext_rep5_complete_prot_fit.dat`$R169.βarr1.D290.βarr1.Ang)
#D297
d297_1 <- get_cfc(neighborhood$`rep1.ARG169@βarr1.step7_production_ext_rep1_complete_prot_fit.dat`$R169.βarr1.D297.βarr1.Ang) 
d297_2 <- 0
d297_3 <- 0
d297_4 <- 0
d297_5 <- get_cfc(neighborhood$`rep5.ARG169@βarr1.step7_production_ext_rep5_complete_prot_fit.dat`$R169.βarr1.D297.βarr1.Ang) 

polarcore <- data.frame(rep1=c(d26_1,d290_1,d297_1),
                        rep2=c(d26_2,d290_2,d297_2),rep3=c(d26_3,d290_3,d297_3),
                        rep4=c(d26_4,d290_4,d297_4),rep5=c(d26_5,d290_5,d297_5), residues = c("D26","D290","D297"))

bar1 <- ggplot(data=polar, aes(x=residues, y=rep1*100)) + geom_bar(stat ="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   xlab("") + 
  ylab("Rep1") + theme(axis.title = element_text(size=10),axis.text = element_text(size=8))
bar2 <- ggplot(data=polar, aes(x=residues, y=rep2*100)) + geom_bar(stat ="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + 
  ylab("Rep2") + theme(axis.title = element_text(size=10),axis.text = element_text(size=8))
bar3 <- ggplot(data=polar, aes(x=residues, y=rep3*100)) + geom_bar(stat ="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + 
  ylab("Rep3") + theme(axis.title = element_text(size=10),axis.text = element_text(size=8))
bar4 <- ggplot(data=polar, aes(x=residues, y=rep4*100)) + geom_bar(stat ="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + 
  ylab("Rep4") + theme(axis.title = element_text(size=10),axis.text = element_text(size=8))
bar5 <- ggplot(data=polar, aes(x=residues, y=rep5*100)) + geom_bar(stat ="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + 
  ylab("Rep5") + theme(axis.title = element_text(size=10),axis.text = element_text(size=8))

overview_polarcore <- grid.arrange(bar1,bar2,bar3,bar4,bar5,ncol=5, left="CF(%)")
ggsave("20210128_polarcore_ctcs.png", plot=overview_polarcore, device=png(), 
       path="~/Desktop/figures/", width=16, height = 4, units = "cm", dpi=300)
dev.off()
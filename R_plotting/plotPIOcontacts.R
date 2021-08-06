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

#####bulk-reading####
setwd("~/Desktop/run002/analyses_complete/VMD_PIO_contacts/cutoff")        #supply path here
temp <- list.files(pattern = "*.list", recursive = TRUE) #list files in subfolders
pio_ctcs <- lapply(temp,function(i){
  read.csv(i,header = FALSE,sep = "\t",stringsAsFactors = FALSE)         #read with skipping gromacs output header
})
names(pio_ctcs) <- basename(list.files(pattern = "*.list", recursive = TRUE)) #strip parent folder names
expnames <- names(pio_ctcs)
replicats <- 5
expnum <- length(pio_ctcs)/replicats

index <- 0
plots <- list()
for (j in 1:replicats) {
  for (i in 1:expnum) {
    pio_ctcs[[index+i]] <- subset(pio_ctcs[[index+i]],select = -c(2)) #prepare for getting long format with pivot_longer
    colnames(pio_ctcs[[index+i]]) <- c("res",seq(1,length(pio_ctcs[[index+i]])-1))
    rownames(pio_ctcs[[index+i]]) <- pio_ctcs[[index+i]]$res
    pio_long <- pivot_longer(pio_ctcs[[index+i]], cols = !res, names_to = "times", values_to = "contacts") 
    plots[[index+i]] <- ggplot(pio_long) + geom_raster(aes(x=as.numeric(times), y=res, fill=factor(contacts))) + 
      theme(legend.position = "none", ) + 
      scale_fill_manual(name= "Contact", values=c("white","red")) + 
      scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
      theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt")) +
      theme(axis.title = element_blank(), axis.text = element_text(size=7))
    }   
    index <- index + expnum
}

ctc_num <- lapply(pio_ctcs, nrow)
ctc_num <- unlist(ctc_num, use.names = FALSE)

overview_salt <- grid.arrange(arrangeGrob(plots[[1]],plots[[3]],plots[[5]],plots[[7]],plots[[9]], 
                                          nrow = 5, heights = ctc_num[seq(1,9,2)], top="\u03B2arr1" , bottom = "Time [ns]"),
                              arrangeGrob(plots[[2]],plots[[4]],plots[[6]],plots[[8]],plots[[10]], 
                                          nrow = 5, heights = ctc_num[seq(2,10,2)], top="NTS1R", bottom = "Time [ns]"),
                              ncol=2,top="")

ggsave("20201123_PIO_contacts.png", plot=overview_salt, device=png(), 
       path="~/Desktop/figures/", width=16, height = 24.7, units = "cm", dpi=300)
dev.off()
# overview_salt <- grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]], 
#                               plots[[6]],plots[[7]],plots[[8]],plots[[9]],plots[[10]], 
#                               nrow = 5, ncol = 2,top="Intermolecular Salt Bridges to PIP2")


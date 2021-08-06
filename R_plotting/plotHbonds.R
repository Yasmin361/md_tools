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
setwd("~/Desktop/run002/analyses_complete/VMD_H_bonds")        #supply path here
arg <- list.files(pattern = "*ARG169_any_intra_arr.list", recursive = TRUE) #list files in subfolders
arg_bonds <- lapply(arg,function(i){
  read.csv(i,header = FALSE,sep = "\t",stringsAsFactors = FALSE)         #read with skipping gromacs output header
})
names(arg_bonds) <- basename(list.files(pattern = "*ARG169_any_intra_arr.list", recursive = TRUE)) #strip parent folder names
replicats <- 5
expnum <- length(arg_bonds)/replicats

index <- 0
plots <- list()
for (j in 1:replicats) {
  for (i in 1:expnum) {
    arg_bonds[[index+i]] <- subset(arg_bonds[[index+i]],select = -c(3,4))
    colnames(arg_bonds[[index+i]]) <- c("res1", "res2",seq(1,length(arg_bonds[[index+i]])-2))
    rownames(arg_bonds[[index+i]]) <- paste(as.character(arg_bonds[[index+i]]$res1), as.character(arg_bonds[[index+i]]$res2))
    arg_long <- pivot_longer(arg_bonds[[index+i]], cols = !c(res1,res2), names_to = "times", values_to = "contacts")
    plots[[index+i]] <- ggplot(arg_long) + geom_tile(aes(x=as.numeric(times), y=res1, fill=factor(contacts))) + 
      ylab("") + xlab("") + theme(legend.position = "none") + 
      scale_fill_manual(name= "Contact", values=c("white","red"),expand = c(0,0)) + 
      scale_x_continuous(limits = c(0,3000), expand = c(0,0)) + 
      theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))
  }   
  index <- index + expnum
}

overview_arg <- grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],nrow=1,ncol=5,left="ARG169 contacts",bottom="Time [ns]")
overview_arg

overview_salt <- grid.arrange(arrangeGrob(plots[[1]],plots[[3]],plots[[5]],plots[[7]],plots[[9]], 
                                          nrow = 5, top="\u03B2arr1" , bottom = "Time [ns]"),
                              arrangeGrob(plots[[2]],plots[[4]],plots[[6]],plots[[8]],plots[[10]], 
                                          nrow = 5, top="NTS1R", bottom = "Time [ns]"),
                              ncol=2,top="Intermolecular Salt Bridges to PIP2")

# overview_salt <- grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]], 
#                               plots[[6]],plots[[7]],plots[[8]],plots[[9]],plots[[10]], 
#                               nrow = 5, ncol = 2,top="Intermolecular Salt Bridges to PIP2")
overview_salt

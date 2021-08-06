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

get_cfc<-function(contacttable){
  trjlength <- length(contacttable)
  ctcs <- length(contacttable[which(contacttable <=4)])
  return(ctcs/trjlength)
}

#####bulk-reading####
setwd("~/Desktop/run002/analyses_complete/mdciao/traces")        #supply path here
temp <- list.files(pattern = "*fit\\.dat", recursive = TRUE) #list files in subfolders
traces <- lapply(temp,function(i){
  read.csv(i,header = FALSE,sep = "",stringsAsFactors = FALSE)         
})

for (i in c(1:6)) {
  traces[[i]][1,1:(ncol(traces[[i]])-1)] <- traces[[i]][1,2:ncol(traces[[i]])]
  traces[[i]] <- traces[[i]][,1:(ncol(traces[[i]])-1)]
  names(traces[[i]]) <- traces[[i]][1,]
  traces[[i]] <- traces[[i]][2:nrow(traces[[i]]),]
  traces[[i]] <- data.frame(traces[[i]])
  traces[[i]] <- lapply(traces[[i]], function(x) as.numeric(x))
}
names(traces) <- c("rep0","rep1","rep2","rep3","rep4","rep5") #to use n=250 for replica 1.

update_geom_defaults("line",list(size=.5))
colors <- c(rep("red",5))
plot_timetrace <- function(tracetable,replica,contact){
  contact_label <- unlist(strsplit(as.character(contact),".",fixed = TRUE))
  replicaname <- paste("rep",replica,sep="")
  traceplot <- ggplot() +
    xlab(paste("Time [µs] (", 
               round(get_cfc(get(contact,get(replicaname,tracetable))),2)*100,
               "%)",sep="")) + # for single traces
    ylab("")+ ylim(c(0,10)) + # for grid
#    ylab(paste(contact_label[1],"-",contact_label[3],"(", contact_label[4],".",contact_label[5],") [\u00C5]",sep="")) + ylim(0,10) +
    geom_hline(yintercept = 4, linetype="dashed") +
    geom_line(aes(x=get("time.ns",get(replicaname,tracetable))/1000,
                  y=rollmean(get(contact,get(replicaname,tracetable)), 40, fill = NA)), 
              size=.5, color=colors[replica]) + #alternative , color=rainbow(5)[replica]
    theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
  return(traceplot)
}

trace1 <- plot_timetrace(traces,1,"V70.βarr1.G301.6.33.Ang")
trace3 <- plot_timetrace(traces,1,"F75.βarr1.M266.5.68.Ang")
trace4 <- plot_timetrace(traces,1,"L71.βarr1.L263.5.65.Ang")


trace_other <- grid.arrange(
  arrangeGrob(  textGrob("Rep2", just="centre"),
                textGrob("Rep3", just="centre"),  
                textGrob("Rep4", just="centre"), 
                textGrob("Rep5", just="centre"), nrow=1),
  arrangeGrob( plot_timetrace(traces,2,"V70.βarr1.G301.6.33.Ang"),
               plot_timetrace(traces,3,"V70.βarr1.G301.6.33.Ang"),
               plot_timetrace(traces,4,"V70.βarr1.G301.6.33.Ang"),
               plot_timetrace(traces,5,"V70.βarr1.G301.6.33.Ang"),ncol=4,
               left = textGrob("V70 - G301",gp=gpar(fontsize=10),rot=90)),
  arrangeGrob(  plot_timetrace(traces,2,"F75.βarr1.M266.5.68.Ang"),
                plot_timetrace(traces,3,"F75.βarr1.M266.5.68.Ang"),
                plot_timetrace(traces,4,"F75.βarr1.M266.5.68.Ang"),
                plot_timetrace(traces,5,"F75.βarr1.M266.5.68.Ang"),ncol=4,
                left = textGrob("F75 - M266",gp=gpar(fontsize=10),rot=90)),
  arrangeGrob(  plot_timetrace(traces,2,"L71.βarr1.L263.5.65.Ang"),
                plot_timetrace(traces,3,"L71.βarr1.L263.5.65.Ang"),
                plot_timetrace(traces,4,"L71.βarr1.L263.5.65.Ang"),
                plot_timetrace(traces,5,"L71.βarr1.L263.5.65.Ang"),ncol=4,
                left = textGrob("L71 - L263",gp=gpar(fontsize=10),rot=90)),
  nrow=4, heights=c(0.08,0.4,0.4,0.4))


trace_other

trace1_plots <- grid.arrange(trace1,trace3,trace4,ncol=3)
trace1_plots

ggsave("20210201_disrupted_rep1.png", plot=trace1_plots, device=png(), 
       path="~/Desktop/figures/", width=16, height = 4, units = "cm", dpi=300)
dev.off()

ggsave("20210201_disrupted_other.png", plot=trace_other, device=png(), 
       path="~/Desktop/figures/", width=16, height = 8, units = "cm", dpi=300)
dev.off()

grid.arrange(
plot_timetrace(traces,1,"N245.βarr1.H172.3.56.Ang"),
plot_timetrace(traces,2,"N245.βarr1.H172.3.56.Ang"),
plot_timetrace(traces,3,"N245.βarr1.H172.3.56.Ang"),
plot_timetrace(traces,4,"N245.βarr1.H172.3.56.Ang"),
plot_timetrace(traces,5,"N245.βarr1.H172.3.56.Ang"),
nrow=1
)

#"R285.βarr1.Q95.NTS1R.Ang"  "L71.βarr1.I170.3.54.Ang"  
#[4] "K284.βarr1.L94.NTS1R.Ang"  "K250.βarr1.S96.NTS1R.Ang"  "Q248.βarr1.Q98.NTS1R.Ang" 
#[7]  "N245.βarr1.H172.3.56.Ang"

#contacts
E66.βarr1.K91.1.60.Ang #meh

F75.βarr1.M266.5.68.Ang #good
V70.βarr1.G301.6.33.Ang #good
#L71.βarr1.G301.6.33.Ang #good
#L71.βarr1.R166.3.50.Ang #good

L71.βarr1.L263.5.65.Ang #second move
 #second move

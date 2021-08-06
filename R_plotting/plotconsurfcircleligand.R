library(circlize)
library(dplyr)
library(viridis)
library(RColorBrewer)

### Make data
####ligand####
setwd("~/Desktop/run002/analyses_complete/ligand_mdc")        #supply path here
ligand_contact <- read.csv("sum_of_replicas_pairs.list", header = FALSE,sep = "\t", stringsAsFactors=FALSE)
ligand <- data.frame(description = as.character(ligand_contact$V1),
                     shortdescription = sub("@.*","",sub("@frag0","", ligand_contact$V1)),
                     ligseq = as.numeric(sub("[A-Z]","", sub("@.*", "", ligand_contact$V1))),
                     ntrseq = as.numeric(sub("[A-Z]","",sub("@.*","",sub(".*-", "", ligand_contact$V1)))),
                     freq=as.numeric(ligand_contact$V2))

###consurf###
setwd("~/Desktop/run002/analyses_complete/consurf")        #supply path here
ntrconsurf <- read.csv("6up7-R/consurf.grades", header = FALSE,skip=15,sep = "\t", stringsAsFactors=FALSE)
ntrconsurf$V3 <- sub("\\s*-",NA,ntrconsurf$V3)
ntrconsurf <- ntrconsurf[complete.cases(ntrconsurf$V3),]
ntrconsurf$V3 <- as.numeric(sub("[A-Z][A-Z][A-Z]","",sub(":.*","",ntrconsurf$V3)))
ntrconsurf$V6 <- as.numeric(sub("[\\*]","",ntrconsurf$V6))
ntrconsurf <- ntrconsurf[complete.cases(ntrconsurf$V6),]
ntrconservation <- c()
for (i in c(1:length(ligand$ntrseq))) {
  ntrconservation[i] <- ntrconsurf[which(ntrconsurf$V3 == ligand[order(ligand$ntrseq),]$ntrseq[i]), ]$V6
} 
ligand <- cbind(ligand[order(ligand$ntrseq),],ntrconservation)

ln <- unique(ligand$ntrseq)
rescon <- c()
for (i in c(1:length(ln))) {
  rescon[i] <- ntrconsurf[which(ntrconsurf$V3 == ln[i]), ]$V6
} 
ln <- c(0,ln)

conservation_palette <- c(brewer.pal(9, "RdBu"))
                   
####calculate segment sum and residue sum####                   
limitres2 <- c(50,58,59,91,99,130,131,136,137,172,183,207,208,228,229,270,295,329,330,334,335,367)
lengthsele2 <- c()
for (i in c(1:11)) {
  lengthsele2[i] <- limitres2[2*i] - limitres2[2*i-1] + 1
}
elements <- c("NTS(8-13)", "NTER","TM1", "TM2","ECL1","TM3","TM4","ECL2","TM5","TM6","ECL3","TM7")
freqsums <- c()
for (i in c(1:11)) {
  freqsums[i] <- sum(subset(ligand,subset = ntrseq <= limitres2[2*i] & ntrseq >= limitres2[2*i-1])$freq)
}
freqsums
segmentsums <- c(sum(freqsums), freqsums)
segmentsums

ligand_by_ntres <- ligand[order(ligand$ntrseq),]

residuesums <- c(0,sum(freqsums))
for (i in c(2:(length(ln)-1))) {
  residuesums[i+1] <- sum(ligand_by_ntres[which(ligand_by_ntres$ntrseq == ln[i]),]$freq) + residuesums[i]
}
residuesums

####dataframe####
m <- data.frame(order = 1:12,
                segment = c("NTS(8-13)", "N-ter", "TM1", "TM2", "ECL1", "TM3", "TM4", "ECL2", "TM5", "TM6", "ECL3", "TM7"),
                V3 = c(1, freqsums),
                V4 = c(freqsums[1], rep(0.001,length(freqsums))),
                V5 = c(freqsums[2], rep(0.001,length(freqsums))),
                V6 = c(freqsums[3], rep(0.001,length(freqsums))),
                V7 = c(freqsums[4], rep(0.001,length(freqsums))),
                V8 = c(freqsums[5], rep(0.001,length(freqsums))),
                V9 = c(freqsums[6], rep(0.001,length(freqsums))),
                V10 = c(freqsums[7], rep(0.001,length(freqsums))),
                V11 = c(freqsums[8], rep(0.001,length(freqsums))),
                V12 = c(freqsums[9], rep(0.001,length(freqsums))),
                V13 = c(freqsums[10], rep(0.001,length(freqsums))),
                V14 = c(freqsums[11], rep(0.001,length(freqsums))),
                vircol = c("gray",viridis(11)),
                stringsAsFactors = FALSE)
df1 <- m[, c(1,2,15)]
m <- m[,-c(1,2,15)]
m <- as.matrix(m[,c(1:12)])
dimnames(m) <- list(orig = df1$segment, dest = df1$segment)
#Sort order of data.frame and matrix for plotting in circos
df1 <- arrange(df1, order)
df1$segment <- factor(df1$segment, levels = df1$segment)
m <- m[levels(df1$segment),levels(df1$segment)]


### Define ranges of circos sectors and their colors (both of the sectors and the links)
df1$xmin <- 0
df1$xmax <- rowSums(m) #+ colSums(m)
n <- nrow(df1)
df1$rcol<-df1$vircol
df1$lcol<-df1$vircol

#### Plot sectors (outer part) ####
par(mar=rep(0,4))
circos.clear()

### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), gap.degree =0,start.degree = 90)

### Sector details
circos.initialize(factors = df1$segment, xlim = cbind(df1$xmin, df1$xmax))

### Plot sectors
circos.trackPlotRegion(ylim = c(0, 1), factors = df1$segment, track.height=0.1,
                       #panel.fun for each sector
                       panel.fun = function(x, y) {
                         #select details of current sector
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         
                         #text direction (dd) and adjusmtents (aa)
                         theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                         dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                         aa = c(1, 0.5)
                         if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                         
                         #plot segment labels
                         circos.text(x=mean(xlim), y=1.7, labels=name, facing = dd, cex=0.6,  adj = aa)
                         
                         #plot main sector in upper track half
                         circos.rect(xleft=xlim[1], ybottom=0.5, xright=xlim[2], ytop=ylim[2], 
                                     col = df1$rcol[i], border=df1$rcol[i])

                         #plot axis
                         circos.axis(labels.cex=0.6, direction = "outside",  
                                     minor.ticks=1)
                         

                       })

set.current.cell(sector.index = "NTS(8-13)", track.index = 1)
#plot residue-wise contribution in lower track half
circos.rect(xleft=residuesums[1], ybottom=0, xright=residuesums[2], ytop=0.5,
            col = "grey", border = "grey")
for (i in c(2:length(ln))) {
  circos.rect(xleft=residuesums[i], ybottom=0, xright=residuesums[i+1], ytop=0.5,
              col = conservation_palette[rescon[i-1]], border = conservation_palette[rescon[i-1]])
}


# ### Plot links (inner part)
# ### Add sum values to df1, marking the x-position of the first links
# ### out (sum1) and in (sum2). Updated for further links in loop below.
# #df1$sum1 <- colSums(m)
# df1$sum2 <- numeric(n)
# 
# ### Create a data.frame of the flow matrix sorted by flow size, to allow largest flow plotted first
# df2 <- cbind(as.data.frame(m),orig=rownames(m),  stringsAsFactors=FALSE)
# df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
#                timevar="dest", time=rownames(m),  v.names = "m")
# df2 <- arrange(df2,desc(m))
# 
# ### Keep only the largest flows to avoid clutter
# df2 <- subset(df2, m > quantile(m,0.5))
# 
# ### Plot links
# for(k in 1:nrow(df2)){
#   #i,j reference of flow matrix
#   i<-match(df2$orig[k],df1$segment)
#   j<-match(df2$dest[k],df1$segment)
#   
#   #plot link
#   circos.link(sector.index1=df1$segment[i], point1=c(df1$sum1[i], df1$sum1[i] + abs(m[i, j])),
#               sector.index2=df1$segment[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(m[i, j])),
#               col = df1$lcol[i])
#   
#   #update sum1 and sum2 for use when plotting the next link
# #  df1$sum1[i] = df1$sum1[i] + abs(m[i, j])
#   df1$sum2[j] = df1$sum2[j] + abs(m[i, j])
# }
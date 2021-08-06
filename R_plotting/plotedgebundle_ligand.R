library(circlize)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(ggraph)

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
ntrconsurf <- read.csv("4jqi-A/consurf.grades", header = FALSE,skip=15,sep = "\t", stringsAsFactors=FALSE)
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
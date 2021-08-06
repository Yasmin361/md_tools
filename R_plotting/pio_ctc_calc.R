
pairsum <- data.frame(description = as.character(sum_of_replicas_pairs$V3),
                     shortdescription = sub(".*.@","",sub("@frag0","", sum_of_replicas_pairs$V3)),
                     freq=as.numeric(sum_of_replicas_pairs$V4))
sum(pairsum[which(pairsum$shortdescription != "βarr1"),]$freq)
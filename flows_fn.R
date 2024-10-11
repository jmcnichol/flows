# flow estimate function 

flows.est <- function(clust.list,zones.seq,zones.sample,burn){
infect.prob.list <- list()
for (i in 1:length(clust.list)) {
  infect.prob.list[[i]] <- computeMatWIW(clust.list[[i]], burnin = burn)
}

trans.probs.df <- list()
for (i in 1:length(infect.prob.list)){
  trans.probs.df[[i]]<-reshape2::melt(as.matrix(infect.prob.list[[i]]), varnames = c('Infector', 'Infectee'), na.rm = TRUE)
}
trans.probs.df <- do.call(rbind,trans.probs.df)

# keep only most probable infector for each infectee
infecteeList<-as.character(unique(trans.probs.df$Infectee))
remove<-numeric()
for (i in 1:length(infecteeList)){
  cols<-which(trans.probs.df$Infectee==infecteeList[i])
  if (length(cols)>1){
    remove<-c(remove,cols[which(trans.probs.df$value[cols]!=max(trans.probs.df$value[cols]))])
  }
}
trans.probs.df<-trans.probs.df[-remove,]

trans.probs.df <- trans.probs.df[-which(trans.probs.df$value==0),]

# assign groups to infector and infectees
infector.group <- infectee.group <- rep(NA,length(trans.probs.df$Infector))
for (i in 1:nrow(trans.probs.df)){
  if (sum(zones.seq[,1] == trans.probs.df$Infectee[i]) > 0 & sum(zones.seq[,1] == trans.probs.df$Infector[i]) > 0) {
    infector.group[i] <- zones.seq[which(zones.seq[,1] == trans.probs.df$Infector[i]),2]
    infectee.group[i] <- zones.seq[which(zones.seq[,1] == trans.probs.df$Infectee[i]),2]
  }
}

group.df <- na.omit(data.frame(Infector=infector.group,Infectee=infectee.group))


A=group.df$Infector;B=group.df$Infectee

Am = matrix(0, nrow=length(unique(zones.sample[,1])), ncol=length(unique(zones.sample[,1])))
for( i in seq_along(A)) {
  Am[A[i],B[i]] <- Am[A[i],B[i]] + 1
}
group.count = Am

top <- pi <- matrix(NA,nrow(group.count),ncol(group.count))
for (a in 1:nrow(group.count)){
  for (b in 1:ncol(group.count)){
    top[a,b] <- group.count[a,b]/(sample.prop[a]*sample.prop[b])
  }
}

for (a in 1:nrow(group.count)){
  for (b in 1:ncol(group.count)){
    pi[a,b] <- top[a,b]/sum(top[-a,-b])
  }
}
return(pi)
}
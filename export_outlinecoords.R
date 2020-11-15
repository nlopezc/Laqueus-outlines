#How to export individual outline coordinates to .txt files
  #Dataset 1
for(i in seq_along(1:40)) {
  write.table(LaqueusCTOutlines.al$coo[[i]], file = paste(i,".txt"), row.names=F, col.names=F)
}
  #Dataset 2--extant specimens
for(i in seq_along(1:99)) {
  write.table(LaqueusOutlines$coo[[i]], file = paste(i,".txt"), row.names=F, col.names=F)
}
  #Dataset 2--fossils only
for(i in seq_along(1:16)) {
  write.table(laqueusfossils.out$coo[[i]], file = paste(i,".txt"), row.names=F, col.names=F)
}


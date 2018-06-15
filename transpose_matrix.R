rm(list=ls())

data.dir <- "/home/pablo/genetica/NeuralNetNeumo/data/"
nas.file <- "repucri16052018_filtro.NAs.txt"
nas.t.file <- "repucri16052018_filtro.NAs.t.txt"
setwd(data.dir)

nas.matrix <- read.table(nas.file, sep="\t", header=T)
write.table(t(nas.matrix), file=nas.t.file,sep="\t", quote=F)

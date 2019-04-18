lsCol<-c()

for (name in 1:nrow(ColFactBact)){
if (ColFactBact[name,1]=="Betaproteobacteria"){
  lsCol<-c(lsCol,"#1f78b4")
}else if (ColFactBact[name,1]=="Bacilli"){
  lsCol<-c(lsCol,"#e31a1c")
}else if (ColFactBact[name,1]=="Actinobacteria"){
  lsCol<-c(lsCol,"#ff7f00")
}else if (ColFactBact[name,1]=="Gammaproteobacteria"){
  lsCol<-c(lsCol,"#6a3d9a")
}else if (ColFactBact[name,1]=="Alphaproteobacteria"){
  lsCol<-c(lsCol,"#b15928")
}else if (ColFactBact[name,1]=="Sphingobacteriia"){
  lsCol<-c(lsCol,"#a6cee3")
}else if (ColFactBact[name,1]=="Flavobacteriia"){
  lsCol<-c(lsCol,"#b2df8a")
}else if (ColFactBact[name,1]=="Unclassified"){
  lsCol<-c(lsCol,"#bfbfbf")
}else if (ColFactBact[name,1]=="Acidobacteria_Gp4"){
  lsCol<-c(lsCol,"#7a7a7a")
}else if (ColFactBact[name,1]=="Verrucomicrobiae"){
  lsCol<-c(lsCol,"#33a02c")
}else if (ColFactBact[name,1]=="Cytophagia"){
  lsCol<-c(lsCol,"#fb9a99")
}else if (ColFactBact[name,1]=="Thermomicrobia"){
  lsCol<-c(lsCol,"#fdbf6f")
}else if (ColFactBact[name,1]=="Spartobacteria"){
  lsCol<-c(lsCol,"#cab2d6")
}else if (ColFactBact[name,1]=="Deltaproteobacteria"){
  lsCol<-c(lsCol,"#ffff99")
}else {lsCol<-c(lsCol,"FALSE")}
}

ColFactBact<-cbind(ColFactBact,lsCol)

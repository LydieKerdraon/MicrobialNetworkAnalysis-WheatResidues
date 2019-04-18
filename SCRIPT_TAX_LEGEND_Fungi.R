lsCol<-c()

for (name in 1:nrow(ColFactFung)){
if (ColFactFung[name,1]=="Sordariomycetes"){
  lsCol<-c(lsCol,"#1f78b4")
}else if (ColFactFung[name,1]=="Dothideomycetes"){
  lsCol<-c(lsCol,"#e31a1c")
}else if (ColFactFung[name,1]=="Pezizomycotina_cls_Incertae_sedis"){
  lsCol<-c(lsCol,"#ff7f00")
}else if (ColFactFung[name,1]=="Leotiomycetes"){
  lsCol<-c(lsCol,"#6a3d9a")
}else if (ColFactFung[name,1]=="Agaricomycetes"){
  lsCol<-c(lsCol,"#b15928")
}else if (ColFactFung[name,1]=="Eurotiomycetes"){
  lsCol<-c(lsCol,"#a6cee3")
}else if (ColFactFung[name,1]=="Pezizomycetes"){
  lsCol<-c(lsCol,"#b2df8a")
}else if (ColFactFung[name,1]=="Unclassified"){
  lsCol<-c(lsCol,"#bfbfbf")
}else if (ColFactFung[name,1]=="Tremellomycetes"){
  lsCol<-c(lsCol,"#33a02c")
}else if (ColFactFung[name,1]=="unclassified_Chytridiomycota"){
  lsCol<-c(lsCol,"#bfbfbf")
}else if (ColFactFung[name,1]=="Orbiliomycetes"){
  lsCol<-c(lsCol,"#fb9a99")
}else if (ColFactFung[name,1]=="Chytridiomycetes"){
  lsCol<-c(lsCol,"#fdbf6f")
}else if (ColFactFung[name,1]=="Agaricostilbomycetes"){
  lsCol<-c(lsCol,"#cab2d6")
}else if (ColFactFung[name,1]=="Cystobasidiomycetes"){
  lsCol<-c(lsCol,"#ffff99")
}else {lsCol<-c(lsCol,"FALSE")}
}

ColFactFung<-cbind(ColFactFung,lsCol)

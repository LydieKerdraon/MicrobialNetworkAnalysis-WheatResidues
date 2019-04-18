
phyloFungAll<-phyloRunFungifilt
phyloBactAll<-phyloRunBactfilt
taxa_names(phyloFungAll)
taxa_names(phyloBactAll)


library(ggplot2)
library(igraph)
library(Matrix)
library(SpiecEasi)
library(phyloseq)
#Define Path

path<-"./Network/"
path2<-paste0(path,"Network_withLabels/")
path3<-paste0(path,"ModulesNetwork/")
path4<-paste0(path,"Tax_network/")
path5<-paste0(path,"Betamat_results/")
path6<-paste0(path,"betweeness/")
dir.create(file.path(path), showWarnings = FALSE)
dir.create(file.path(path2), showWarnings = FALSE)
dir.create(file.path(path3), showWarnings = FALSE)
dir.create(file.path(path4), showWarnings = FALSE)
dir.create(file.path(path5), showWarnings = FALSE)
dir.create(file.path(path6), showWarnings = FALSE)


x<-data.frame()
COND<-list('2016-07-Jul.','2016-10-Above ground','2016-12-Above ground','2017-02-Above ground','2016-10-Soil Contact','2016-12-Soil Contact','2017-02-Soil Contact','2017-10-Above ground','2017-12-Above ground','2018-02-Above ground','2017-07-Jul.','2017-10-Soil Contact','2017-12-Soil Contact','2018-02-Soil Contact')
for (i in COND){
  
## Condition
  
Bact<-subset_samples(phyloBactAll, Condspiec == i)
Fung<-subset_samples(phyloFungAll, Condspiec == i)

condition <- function(x) { x > 0 } 
taxafilt <- genefilter_sample(Bact, condition, 1)
Bact<-prune_taxa(taxafilt, Bact)
taxafilt <- genefilter_sample(Fung, condition, 1)
Fung<-prune_taxa(taxafilt, Fung)

taxafilt <- genefilter_sample(Bact, condition, 6)
Bact.occ8<-prune_taxa(taxafilt, Bact)
taxafilt <- genefilter_sample(Fung, condition, 6)
Fung.occ8<-prune_taxa(taxafilt, Fung)

sample_names(Fung.occ8)
sample_names(Bact.occ8)

set.seed(1)

## Network

spiec <- spiec.easi(list(Bact.occ8, Fung.occ8), method='mb', nlambda=100,
                    lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))

dtypex <- c(rep(1,ntaxa(Bact.occ8)), rep(2,ntaxa(Fung.occ8)))
dtype <- c(taxa_names(Bact.occ8), taxa_names(Fung.occ8))

list(name=dtype)

ig<-adj2igraph(symBeta(getOptBeta(spiec)))
igRefit<-adj2igraph(getRefit(spiec))
CoordDelete <- layout_nicely(igRefit)

nodedegree<-(degree(ig)==0)
names(nodedegree)<-dtype
NodeInteracting<-table(nodedegree)["FALSE"]
NodeNonInteracting<-table(nodedegree)["TRUE"]


##Network with weight

#pdf(paste0(path2,i,".pdf"), width = 30,height = 30)

set.seed(1)
plot(delete.vertices((adj2igraph(getRefit(spiec), vertex.attr=list(name=dtype,color=dtypex,size=9,label.cex=4),
                                edge.attr=list(width=13,color=(ifelse(E(ig)$weight > 0.0000, 'green', 'red'))))),names(which(nodedegree == "TRUE"))),)+title(i,cex.main = 4)


#dev.off()



######Modularite et Modules
mod.net<-as(spiec$refit$stars,"matrix")
colnames(mod.net) <- rownames(mod.net) <- dtype

testadj<-graph.adjacency(mod.net, mode='undirected', add.rownames = TRUE)
ig2<-adj2igraph(symBeta(getOptBeta(spiec)))

modules =cluster_fast_greedy(testadj)
modules
modularityNote<-modularity(modules)


V(testadj)$color=modules$membership

#####networks with Modules

#pdf(paste0(path3,i,"_modules.pdf"), width = 30,height = 30)
plot(testadj, col = modules, vertex.size = 9, vertex.label.cex=4, vertex.label = dtype, layout=layout.auto)
#dev.off()


#---------Network with colors------#
##.define Shape / colors 

dtypeShape <- c(rep("circle",ntaxa(Bact.occ8)), rep("square",ntaxa(Fung.occ8)))

dtype


ColBact<-c(tax_table(Bact.occ8)[,"Class"])
ColFung<-c(tax_table(Fung.occ8)[,"Class"])

ColFactBact <- transform( ColBact,Col=as.numeric(factor(ColBact)))
ColFactFung <- transform( ColFung,Col=as.numeric(factor(ColFung)))

source("SCRIPT_TAX_LEGEND_Fungi.R")
source("SCRIPT_TAX_LEGEND_Bact.R")

my_color=c(as.vector(ColFactBact[,3]),as.vector(ColFactFung[,3]))
names(my_color)<-c(ColBact,ColFung)



##Make Network

pdf(paste0(path4,i,"_Colored.pdf"), width = 30,height = 30)

set.seed(1)
plot(delete.vertices((adj2igraph(getRefit(spiec), vertex.attr=list(name=dtype,color=my_color, label.cex=0.2,label.color=my_color,shape=dtypeShape,size=7),
                                 edge.attr=list(width=13,color=(ifelse(E(ig)$weight > 0.0000, 'green', 'red'))))),names(which(nodedegree == "TRUE"))),)+title(i,cex.main = 4)     

dev.off()


WeightNeg<-table(E(ig)$weight>0)["FALSE"]
WeightPos<-table(E(ig)$weight>0)["TRUE"]


betaMat=as.matrix(symBeta(getOptBeta(spiec)))
betaMat[betaMat==0] <- NA
betaMat[lower.tri(betaMat, diag = TRUE)] <- NA

colnames(betaMat) <- rownames(betaMat) <- dtype
testbetamat <- na.omit(data.frame(as.table(betaMat)))
BB<-as.data.frame(dtype)
BB2<-cbind(BB,"","")
colnames(BB2)<-colnames(testbetamat)
testbetamat2<-rbind(testbetamat,BB2)
#write.csv(testbetamat2,paste(path5,i,"_testbetamat.csv"))

taxBact<-cbind(as(tax_table(Bact.occ8),"matrix"),"","")
taxFungi<-as(tax_table(Fung.occ8),"matrix")
taxAll<-rbind(taxBact,taxFungi)
testsuite<-(degree(ig))
names(testsuite)<-dtype
taxAll_withDegree<- merge(as.data.frame(taxAll), as.data.frame(testsuite), by='row.names', all=TRUE)
datamodules<-data.frame(modules$names, modules$membership)

#write.csv(datamodules,paste0(path3,i,"_modules.csv"))
#write.csv(taxAll_withDegree,paste0(path,i,"_tax.csv"))

bet<-betweenness(ig, v = V(ig), directed = FALSE, weights = NA,
                 nobigint = TRUE, normalized = FALSE)

names(bet) <-  dtype

deg<-degree(ig)
names(deg)<-dtype
BetTax<-cbind(as.data.frame(deg),as.data.frame(bet),taxAll)

assign(paste("PlotBetweeness_",i, sep = ""),ggplot(BetTax, aes_string(x="deg", y="bet"))+geom_point(aes_string(shape="Kingdom", color="Class"),size=2)+theme_bw()+scale_colour_manual(values=my_color)+scale_shape_manual(values=c(16, 15))+ggtitle(i)+ theme(legend.position="none")+xlim(0,13)+ylim(0,2600)
)

write.csv(BetTax, paste0(path6,i,"_Bet.csv"))
x<-rbind(x,data.frame(i, ntaxa(Bact),ntaxa(Fung),ntaxa(Bact.occ8),ntaxa(Fung.occ8), getStability(spiec),sum(getRefit(spiec))/2,NodeInteracting,NodeNonInteracting,modularity(modules),length(modules),diameter(ig, directed = F,weights = NA),
                      names(bet)[which.max(bet)], mean_distance(ig, directed = FALSE, unconnected = TRUE),WeightNeg,WeightPos,mean(bet),mean(deg),sum(deg)/NodeInteracting))
}

write.csv(x, paste(path,"donnees_reseaux.csv"))

library(cowplot)
plot_grid(`PlotBetweeness_2016-07-Jul.`,`PlotBetweeness_2016-10-Above ground`,`PlotBetweeness_2016-10-Soil Contact`,
          `PlotBetweeness_2016-12-Above ground`,`PlotBetweeness_2016-12-Soil Contact`,
          `PlotBetweeness_2017-02-Above ground`,`PlotBetweeness_2017-02-Soil Contact`,
          `PlotBetweeness_2017-07-Jul.`,`PlotBetweeness_2017-10-Above ground`,`PlotBetweeness_2017-10-Soil Contact`,
          `PlotBetweeness_2017-10-Above ground`,`PlotBetweeness_2017-10-Soil Contact`,
          `PlotBetweeness_2017-12-Above ground`,`PlotBetweeness_2017-12-Soil Contact`,
          `PlotBetweeness_2018-02-Above ground`,`PlotBetweeness_2018-02-Soil Contact`, ncol=3)




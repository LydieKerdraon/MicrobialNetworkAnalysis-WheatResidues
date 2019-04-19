#Data used for these analyses are given in the R environment 16S_Data.RData


## Packages required + custom functions
Packages <- c("phyloseq", "data.table", "ggplot2", "plyr","dplyr","reshape2","grid",
              "gridExtra","scales","dplyr", "ggpubr","vegan","multcompView","rcompanion","betapart")

lapply(Packages, library, character.only = TRUE)


########################
# Dataset constitution #
########################

phyloRunBact  <- phyloseq(otu_table(t(seqTabBact), taxa_are_rows=TRUE),
                         sample_data(data.tableBact),
                         tax_table(taxBact))
taxa_names(phyloRunBact) <- paste0("B", seq(ntaxa(phyloRunBact)))


######---- Suppression Mock ----######
phyloMock<-subset_samples(phyloRunBact, Statut == "mock")
phyloRunBact<-subset_samples(phyloRunBact, Statut != "mock")

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phyloRunBact, condition, 1)
phyloRunBact<-prune_taxa(taxaToKeep, phyloRunBact)



##
phyloRunBactfilt <- transform_sample_counts(phyloRunBact,function(x) ifelse(x>=0.003*sum(x),x,0))

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phyloRunBactfilt, condition, 1)
phyloRunBactfilt<-prune_taxa(taxaToKeep, phyloRunBactfilt)

table(tax_table(phyloRunBactfilt)[, "Phylum"])
phyloRunBactfilt <- subset_taxa(phyloRunBactfilt, !Phylum %in% c("", "Cyanobacteria/Chloroplast", "Unclassified"))
phyloRunBactfilt<-subset_samples(phyloRunBactfilt, Time != "2016-07")
phyloBactnorm <- transform_sample_counts(phyloRunBactfilt, function(x) round(x/sum(x) *100000 ))




###########################
####--- Alpha - div ---####
###########################

measures=c("Observed", "Shannon")
Rich<-estimate_richness(phyloBactnorm,measures = measures)
Rich2<-cbind(as(sample_data(phyloBactnorm),"matrix"),Rich)
Rich2$Month.ent <- factor(Rich2$Month.ent, levels = c("Jul.", "Oct.", "Dec.", "Feb."))

my_comparisons <- list( c("Jul.","Oct."),c("Oct." , "Dec."),c("Dec." , "Feb."))

## Soil contact
RichSoil<-subset(Rich2, Soil != "Above ground")

ShannonSC<-ggplot(RichSoil, aes(x=Month.ent, y=Shannon,fill=Inoc)) +
  stat_compare_means(aes(group = Inoc), label = "p.signif",label.y=4.5)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y=4.8)+ # Add pairwise comparisons p-value
  geom_boxplot(lwd=0.15, width = 0.5)+scale_fill_manual(values=c("#FFCC66","#99CCFF"))+
  ylab(paste("Shannon index (Bacterial diversity)")) + xlab("") + ylim(0,5)+facet_grid(Saison~.)+theme_bw()
ShannonSC

ObservedSC<-ggplot(RichSoil, aes(x=Month.ent, y=Observed ,fill=Inoc)) +
  stat_compare_means(aes(group = Inoc), label = "p.signif",label.y=100)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y=103)+ # Add pairwise comparisons p-value
  geom_boxplot(lwd=0.15, width = 0.5)+scale_fill_manual(values=c("#FFCC66","#99CCFF"))+
  ylab(paste("Number of ASV (Bacterial diversity)")) + xlab("")+ylim(0,105)+facet_grid(Saison~.)+theme_bw()
ObservedSC

## Above ground
RichAG<-subset(Rich2, Soil != "Soil Contact")

ShannonAG<-ggplot(RichAG, aes(x=Month.ent, y=Shannon,fill=Inoc)) +
  stat_compare_means(aes(group = Inoc), label = "p.signif",label.y=4.5)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y=4.8)+ # Add pairwise comparisons p-value
  geom_boxplot(lwd=0.15, width = 0.5)+scale_fill_manual(values=c("#FFCC66","#99CCFF"))+
  ylab(paste("Shannon index (Bacterial diversity)")) + xlab("") + ylim(0,5)+facet_grid(Saison~.)+theme_bw()
ShannonAG

ObservedAG<-ggplot(RichAG, aes(x=Month.ent, y=Observed ,fill=Inoc)) +
  stat_compare_means(aes(group = Inoc), label = "p.signif",label.y=100)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y=103)+ # Add pairwise comparisons p-value
  geom_boxplot(lwd=0.15, width = 0.5)+scale_fill_manual(values=c("#FFCC66","#99CCFF"))+
  ylab(paste("Number of ASV (Bacterial diversity)")) + xlab("")+ylim(0,105)+facet_grid(Saison~.)+theme_bw()
ObservedAG

##########################
####--- Beta - div ---####
##########################

#Clustering
require(ape)
sampleType <- get_variable(phyloBactnorm, "Condition.Time")
palette <- hue_pal()(length(levels(sampleType)))
tipColor = col_factor(palette, levels = levels(sampleType))(sampleType)
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
d <- distance(phyloBactnorm, method="bray", type="samples")
hell.hclust     <- as.phylo(hclust(d, method="average"))
hell.hclust$tip.label<-as(sample_data(phyloBactnorm)$Condition.Time,"vector")
plot(hell.hclust, edge.lty = 1,
     tip.color =tipColor,  show.tip.label = TRUE,direction = "downwards")

#MDS

Bact.Jul<-subset_samples(phyloBactnorm, Month.ent == "Jul.")
OTU.ordBact<-ordinate(Bact.Jul,"MDS","bray")
MdsJul=plot_ordination(Bact.Jul,OTU.ordBact,type="sample",color="Month.Inoc", shape="Saison",axes = c(1,2))+
  scale_shape_manual(values = c(17)) + geom_point(size=3)+
  theme(panel.spacing = unit(1, "lines"))+ 
  scale_color_manual(values = c("#009900","#66ff66"))
MdsJul


Bact.withoutJul<-subset_samples(phyloBactnorm, Month.ent != "Jul.")
OTU.ordBact<-ordinate(Bact.withoutJul,"MDS","bray")
MdsBact=plot_ordination(Bact.withoutJul,OTU.ordBact,type="sample",color="Month.Inoc", shape="Saison",axes = c(1,2))+
  scale_shape_manual(values = c(19,17)) + geom_point(size=3)+
  theme(panel.spacing = unit(1, "lines"))+ 
  scale_color_manual(values = c("#cc0000","#ff9999","#0033cc","#99b3ff","#0d0d0d","#b3b3b3"))
MdsBact+facet_grid(~Soil)



#Permanova
set.seed(1)
Bact.withoutJuly<-subset_samples(phyloBactnorm, Soil != "Jul.")
adonis2(distance(Bact.withoutJuly, "bray") ~ Soil +Month.ent +Saison+ Inoc  ,by="margin", data = as(sample_data(Bact.withoutJuly), "data.frame"))

set.seed(1)
Bact.withJuly<-subset_samples(phyloBactnorm, Soil == "Jul.")
adonis2(distance(Bact.withJuly, "bray") ~  Inoc ,by="margin", data = as(sample_data(Bact.withJuly), "data.frame"))


#########################
####---  Heatmap  ---####
#########################

Bactglom<-tax_glom(phyloBactnorm,"Genus")
test   <- subset_taxa(Bactglom, Genus !="Unclassified")
TopNOTUs <- names(sort(taxa_sums(test), TRUE)[1:50])
ent10   <- prune_taxa(TopNOTUs, Bactglom)

nameX<- rev(tax_table(ent10)[,"Genus"][order(tax_table(ent10)[,"Genus"]),])

sample_data(ent10)$Soil<-factor(sample_data(ent10)$Soil, levels = c("Jul.","Above ground","Soil Contact"))

plot_heatmap(ent10,  method=NULL, taxa.label = "Genus", taxa.order = taxa_names(nameX), low="#F9F8F8", high="#000000", na.value = "#FFFFFF")+
  facet_grid(~Saison+Soil+Month.Inoc, scales="free")+
  theme(panel.spacing = unit(0, "lines"))+
  theme(axis.text.y = element_text(face="italic",size=8))


####################
## Abundance plot ##
####################

GlomClass<-tax_glom(phyloBactnorm,"Class")
condition<-function(x){x>0}
taxaToKeep<-genefilter_sample(GlomClass,condition,1)
GlomClass<-prune_taxa(taxaToKeep,GlomClass)
GlomClass
sample_data(GlomClass)$Month.ent <- factor(sample_data(GlomClass)$Month.ent, levels = c("Jul.","Oct.", "Dec.", "Feb."))


Genus<-c("Betaproteobacteria",
         'Bacilli',
         'Actinobacteria',
         'Gammaproteobacteria',
         'Alphaproteobacteria',
         'Sphingobacteriia',
         'Flavobacteriia',
         "Verrucomicrobiae",
         'Cytophagia')



my_comparisons <- list( c("Oct." , "Dec."),c("Dec." , "Feb."))

plot_list = list()
library(ggpubr)
for (i in Genus){
  GlomClass2<-subset_taxa(GlomClass,Class==i)
  GlomClass2<-subset_samples(GlomClass2,Soil !="Jul.")
  dat<-psmelt(GlomClass2)#createdataframe
  
  p <- ggplot(dat,aes(x=Month.ent,y=Abundance, fill= Soil))+
    geom_boxplot(outlier.alpha = 0, alpha=0.25)+
    stat_compare_means(aes(group = Soil), method='wilcox.test',label = "p.signif",label.y=7000)+
    # stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y=10000)+
    scale_x_discrete(name=NULL)+scale_y_continuous(name=NULL)+ ggtitle(i)+theme(plot.title = element_text(size = 12))+
    theme(strip.text.x = element_text(size = 12), legend.position = "none")+
    theme(axis.text = element_text(size = 12))+
    theme(panel.spacing = unit(0.3, "lines"), strip.background = element_rect(colour="black", fill="#E1E1E1", 
                                                                              size=1.5, linetype="solid"))+
    facet_grid(.~Saison)
  
  plot_list[[i]]=p
}
cowplot::plot_grid(plotlist = plot_list)



###################
## Mock Analysis ##
###################
Mocknorm <- transform_sample_counts(phyloMock, function(x) round(x/sum(x) *100000 ))

TopNOTUs <- names(sort(taxa_sums(Mocknorm), TRUE)[1:9])
ent10   <- prune_taxa(TopNOTUs, Mocknorm)

plot_heatmap(ent10, method = NULL, taxa.label = "Genus", taxa.order = c('B1','B120','B449','B26','B477',
                                                                        'B14','B285','B29','B15'))+facet_grid(~Sample, scales="free")

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(Mocknorm, condition, 1)
Mocknorm<-prune_taxa(taxaToKeep, Mocknorm)
plot_list = list()
for (Name in sample_names(Mocknorm)){ 
  phyloEch<- subset_samples(Mocknorm, sample_names(Mocknorm) == Name)
  
  tdt = data.table(tax_table(phyloEch),
                   TotalCounts = taxa_sums(phyloEch),
                   OTU = taxa_names(phyloEch))
  
  taxcumsum = tdt[, .N, by = TotalCounts]
  setkey(taxcumsum, TotalCounts)
  taxcumsum[, CumSum := cumsum(N)]
  # Define the plot
  p=pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
    geom_point() +
    xlab("Filtering threshold, minimum total counts (log2)") +
    ylab("ASV filtered") +
    ggtitle(Name) + geom_vline(xintercept = 300, color="red",size=1.5)+ scale_x_continuous(trans='log2', limits = c(1,38000))+theme_bw()
  
  plot_list[[Name]]=p
}
cowplot::plot_grid(plotlist = plot_list, ncol = 1, nrow=2)

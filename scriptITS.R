#Data used for these analyses are given in the R environment ITS_Data.RData


## Packages required + custom functions
Packages <- c("phyloseq", "data.table", "ggplot2", "plyr","dplyr","reshape2","grid",
              "gridExtra","scales","dplyr", "ggpubr","vegan","multcompView","rcompanion","betapart")

lapply(Packages, library, character.only = TRUE)


########################
# Dataset constitution #
########################

data.tableFungi$Month <- factor(data.tableFungi$Month, levels = c("7", "10", "12", "2"))
data.tableFungi$Month.ent <- factor(data.tableFungi$Month.ent, levels = c("Jul.", "Oct.", "Dec.", "Feb."))
data.tableFungi$Month.Inoc <- factor(data.tableFungi$Month.Inoc, levels = c("Jul.-I", "Jul.-NI","Oct.-I", "Oct.-NI", "Dec.-I","Dec.-NI", "Feb.-I","Feb.-NI"))
data.tableFungi$Position <- factor(data.tableFungi$Position, levels = c("1", "2", "3", "4","5","6","7","8","9","10","11","12","13","14","15"))
taxFungi <- replace(taxFungi, is.na(taxFungi), "Unclassified")

##
phyloRunFungi <- phyloseq(otu_table(t(seqTabFungi), taxa_are_rows=TRUE),
                          sample_data(data.tableFungi),
                          tax_table(taxFungi))

taxa_names(phyloRunFungi) <- paste0("F", seq(ntaxa(phyloRunFungi)))
tax_table(phyloRunFungi)[,1] <- gsub("k__", "", tax_table(phyloRunFungi)[,1]);tax_table(phyloRunFungi)[,2] <- gsub("p__", "", tax_table(phyloRunFungi)[,2]);tax_table(phyloRunFungi)[,3] <- gsub("c__", "", tax_table(phyloRunFungi)[,3]);tax_table(phyloRunFungi)[,4] <- gsub("o__", "", tax_table(phyloRunFungi)[,4]);tax_table(phyloRunFungi)[,5] <- gsub("f__", "", tax_table(phyloRunFungi)[,5]);tax_table(phyloRunFungi)[,6] <- gsub("g__", "", tax_table(phyloRunFungi)[,6])


##
phyloMock<-subset_samples(phyloRunFungi, Statut == "mock") #see Mock Analysis section
phyloRunFungi<-subset_samples(phyloRunFungi, Statut != "mock")

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phyloRunFungi, condition, 1)
phyloRunFungi<-prune_taxa(taxaToKeep, phyloRunFungi)


##
phyloRunFungifilt <- transform_sample_counts(phyloRunFungi,function(x) ifelse(x>=0.003*sum(x),x,0)) #see Mock Analysis section

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phyloRunFungifilt, condition, 1)
phyloRunFungifilt<-prune_taxa(taxaToKeep, phyloRunFungifilt)

table(tax_table(phyloRunFungifilt)[, "Phylum"])
phyloRunFungifilt <- subset_taxa(phyloRunFungifilt, !Phylum %in% c("", "Cyanobacteria/Chloroplast", "Unclassified"))
taxa_names(phyloRunFungifilt)
phyloFungnorm <- transform_sample_counts(phyloRunFungifilt, function(x) round(x/sum(x) *100000 ))


###########################
####--- Alpha - div ---####
###########################

measures=c("Observed", "Shannon")
Rich<-estimate_richness(phyloFungnorm,measures = measures)
Rich2<-cbind(as(sample_data(phyloFungnorm),"matrix"),Rich)
Rich2$Month.ent <- factor(Rich2$Month.ent, levels = c("Jul.", "Oct.", "Dec.", "Feb."))

my_comparisons <- list( c("Jul.","Oct."),c("Oct." , "Dec."),c("Dec." , "Feb."))

## Soil contact
RichSoil<-subset(Rich2, Soil != "Above ground")

ShannonSC<-ggplot(RichSoil, aes(x=Month.ent, y=Shannon,fill=Inoc)) +
  stat_compare_means(aes(group = Inoc), label = "p.signif",label.y=3.6)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y=4)+ # Add pairwise comparisons p-value
  geom_boxplot(lwd=0.15, width = 0.5)+scale_fill_manual(values=c("#FFCC66","#99CCFF"))+
  ylab(paste("Shannon index (Fungal diversity)")) + xlab("") + ylim(0,4.3)+facet_grid(Saison~.)+theme_bw()
ShannonSC

ObservedSC<-ggplot(RichSoil, aes(x=Month.ent, y=Observed ,fill=Inoc)) +
  stat_compare_means(aes(group = Inoc), label = "p.signif",label.y=53)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y=58)+ # Add pairwise comparisons p-value
  geom_boxplot(lwd=0.15, width = 0.5)+scale_fill_manual(values=c("#FFCC66","#99CCFF"))+
  ylab(paste("Number of ASV (Fungal diversity)")) + xlab("")+ylim(0,61)+facet_grid(Saison~.)+theme_bw()
ObservedSC

## Above ground
RichAG<-subset(Rich2, Soil != "Soil Contact")

ShannonAG<-ggplot(RichAG, aes(x=Month.ent, y=Shannon,fill=Inoc)) +
  stat_compare_means(aes(group = Inoc), label = "p.signif",label.y=3.6)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y=4)+ # Add pairwise comparisons p-value
  geom_boxplot(lwd=0.15, width = 0.5)+scale_fill_manual(values=c("#FFCC66","#99CCFF"))+
  ylab(paste("Shannon index (Fungal diversity)")) + xlab("") + ylim(0,4.3)+facet_grid(Saison~.)+theme_bw()
ShannonAG

ObservedAG<-ggplot(RichAG, aes(x=Month.ent, y=Observed ,fill=Inoc)) +
  stat_compare_means(aes(group = Inoc), label = "p.signif",label.y=53)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y=58)+ # Add pairwise comparisons p-value
  geom_boxplot(lwd=0.15, width = 0.5)+scale_fill_manual(values=c("#FFCC66","#99CCFF"))+
  ylab(paste("Number of ASV (Fungal diversity)")) + xlab("")+ylim(0,61)+facet_grid(Saison~.)+theme_bw()
ObservedAG

##########################
####--- Beta - div ---####
##########################

#Clustering
require(ape)
sampleType <- get_variable(phyloFungnorm, "Condition.Time")
palette <- hue_pal()(length(levels(sampleType)))
tipColor = col_factor(palette, levels = levels(sampleType))(sampleType)
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
d <- distance(phyloFungnorm, method="bray", type="samples")
hell.hclust     <- as.phylo(hclust(d, method="average"))
hell.hclust$tip.label<-as(sample_data(phyloFungnorm)$Condition.Time,"vector")
plot(hell.hclust, edge.lty = 1,
     tip.color =tipColor,  show.tip.label = TRUE,direction = "downwards")

#MDS
sample_data(phyloFungnorm)$Soil<- factor(sample_data(phyloFungnorm)$Soil, levels = c("Jul.", "Above ground", "Soil Contact"))
OTU.ordITS<-ordinate(phyloFungnorm,"MDS","bray")
MdsITS=plot_ordination(phyloFungnorm,OTU.ordITS,type="sample",color="Month.Inoc", shape="Saison",axes=1:3)+
  scale_shape_manual(values = c(16,17)) + geom_point(size=3)+
  theme(panel.spacing = unit(1, "lines"))+ facet_grid(.~Soil)+
  scale_color_manual(values = c("#009900","#66ff66","#cc0000","#ff9999","#0033cc","#99b3ff","#0d0d0d","#b3b3b3"))
MdsITS


#Permanova
set.seed(1)
Fung.withoutJuly<-subset_samples(phyloFungnorm, Soil != "Jul.")
adonis2(distance(Fung.withoutJuly, "bray") ~ Soil +Month.ent +Saison+ Inoc  ,by="margin", data = as(sample_data(Fung.withoutJuly), "data.frame"))

set.seed(1)
Fung.withJuly<-subset_samples(phyloFungnorm, Soil == "Jul.")
adonis2(distance(Fung.withJuly, "bray") ~ Saison+ Inoc ,by="margin", data = as(sample_data(Fung.withJuly), "data.frame"))


#########################
####---  Heatmap  ---####
#########################

Funglom<-tax_glom(phyloFungnorm,"Genus")
test   <- subset_taxa(Funglom, Genus !="Unclassified")
test   <- subset_taxa(test, Genus !="unclassified_Phaeosphaeriaceae")
TopNOTUs <- names(sort(taxa_sums(test), TRUE)[1:30])
ent10   <- prune_taxa(TopNOTUs, Funglom)


nameX<- rev(tax_table(ent10)[,"Genus"][order(tax_table(ent10)[,"Genus"]),])
sample_data(ent10)$Soil<-factor(sample_data(ent10)$Soil, levels = c("Jul.","Above ground","Soil Contact"))


plot_heatmap(ent10,  method=NULL, taxa.label = "Genus", taxa.order = taxa_names(nameX), low="#F9F8F8", high="#000000", na.value = "#FFFFFF")+
  facet_grid(~Saison+Soil+Month.Inoc, scales="free")+
  theme(panel.spacing = unit(0, "lines"))+
  theme(axis.text.y = element_text(face="italic",size=10))


####################
## Abundance plot ##
####################

## Zymo
GlomG<-tax_glom(phyloFungnorm,"Genus")
require(microbiome)

Zymo <- transform_sample_counts(GlomG, function(x) round(x/sum(x) *100 ))

sample_data(Zymo)$Soil<-factor(sample_data(Zymo)$Soil, levels =c ("Jul.","Above ground","Soil Contact"))


boxplot_abundance(Zymo, x="Month.ent", y="F11", violin = FALSE, na.rm = FALSE,show.points = FALSE)+
  aes(color=Inoc)+scale_colour_manual(values=c("#db0416","#0fb507"))+labs(x=NULL, y="Percentage of reads")+
  facet_grid(Saison~Soil, scales="free_x",space = "free")+theme_bw()+ 
  theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12),
        axis.title.y = element_text(size=12),strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))+stat_compare_means(aes(group = Inoc), label = "p.signif",label.y=40)


##Classes
GlomClass<-tax_glom(phyloFungnorm,"Class")

condition<-function(x){x>0}
taxaToKeep<-genefilter_sample(GlomClass,condition,1)
GlomClass<-prune_taxa(taxaToKeep,GlomClass)
sample_data(GlomClass)$Month.ent <- factor(sample_data(GlomClass)$Month.ent, levels = c("Jul.","Oct.", "Dec.", "Feb."))

classes<-c('Dothideomycetes','Pezizomycetes','Agaricomycetes','Sordariomycetes','Pezizomycotina_cls_Incertae_sedis',
         'Leotiomycetes','Orbiliomycetes')

plot_list = list()
library(ggpubr)
for (i in classes){
  Boxpp2<-subset_taxa(GlomClass,Class==i)
  Boxpp2<-subset_samples(Boxpp2,Soil !="Jul.")
  dat<-psmelt(Boxpp2)#createdataframe
  
  p <- ggplot(dat,aes(x=Month.ent,y=Abundance, fill= Soil))+
    geom_boxplot(outlier.alpha = 0, alpha=0.25)+
    stat_compare_means(aes(group = Soil), label = "p.signif",label.y=7000)+
    scale_x_discrete(name=NULL)+scale_y_continuous(name=NULL)+ ggtitle(i)+ theme(plot.title = element_text(size = 12))+
    theme(strip.text.x = element_text(size = 12), legend.position = "bottom")+
    theme(axis.text = element_text(size = 12))+
    theme(panel.spacing = unit(0.3, "lines"), strip.background = element_rect(colour="black", fill="#E1E1E1", 
                                                                              size=1.5, linetype="solid"))+
    facet_grid(.~Saison)
  
  plot_list[[i]]=p
}

cowplot::plot_grid(plotlist = plot_list)



####################
## Mock Analysis ##
####################

sample_data(phyloMock)
taxa_names(phyloMock) <- paste0("MF", seq(ntaxa(phyloMock)))


Mocknorm <- transform_sample_counts(phyloMock, function(x) round(x/sum(x) *100000 ))

condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(Mocknorm, condition, 1)
Mocknorm<-prune_taxa(taxaToKeep, Mocknorm)
sample_names(Mocknorm)<-c("Mock Run1","Mock Run3")
sample_data(Mocknorm)


#Filtering Threshold
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


#Heatmap
TopNOTUs <- names(sort(taxa_sums(Mocknorm), TRUE)[1:40])
Top40   <- prune_taxa(TopNOTUs, Mocknorm)

plot_heatmap(Top40, method = NULL,taxa.label = "Genus" )+facet_grid(~Sample, scales="free")


##################################
#Step 1: Prep packages and input data

#several packages that need separate installation
# If a package is not installed, it will be installed from CRAN.
# First select the packages of interest
packages <- c("MASS", "GUniFrac", "ggpubr", "randomForest", "ggdendro", "ggrepel", "agricolae", "picante", "pheatmap", "igraph", "rgexf", 
              "ggalluvial", "ggh4x", "FSA", "gridExtra", "aplot", "tidytree", "microeco")
# Now check or install
for(x in packages){
  if(!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
}

#load the packages
library(MASS) #linear discriminant analysis
library(GUniFrac) #UniFrac distance matrix
library(ggpubr) #some plotting functions
library(randomForest) #random forest analysis
library(ggdendro) #plotting clustering dendrogram
library(ggrepel) #reduce the text overlap in the plot
library(agricolae) #multiple comparisons in anova
library(picante) #Faith's phylogenetic alpha diversity
library(pheatmap) #correlation heatmap with clustering dendrogram 
library(tidytree) #plot the taxonomic tree
library(igraph) #network related operations
library(rgexf) #save network with gexf style
library(ggalluvial) #alluvial plot
library(ggh4x) #multiple facets
library(FSA) #Dunn's Kruskal-Wallis Multiple Comparisons
library(gridExtra) #merge plots
library(aplot) #merge plots
library(magrittr) #set of operators which promote semantics 
library(ggplot2) #plotting graphs
install.packages("readxl")
library(readxl)
install.packages("microeco")
library(microeco)
library(dplyr)

install.packages("forcats")
library(forcats) #This package makes organization much easier, by allowing easier factor relevelling. 

if (!requireNamespace("BiocManager", quietly = TRUE, force = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("microbiome")
install.packages("file2meco", repos = BiocManager::repositories())

install.packages("rcompanion")

library(file2meco)
library(phyloseq)
library(microbiome)
library(rcompanion)

#Import the data

#Diets, n = 2 per treatment
dietsampleinfo <- read.csv('DietMetadataSeYe.csv', sep=',', header=TRUE, row.names=1, check.names=F)
dietasvtable <- read.csv('DietASV.csv', sep=',', header=TRUE, row.names=1, check.names=F)
dietasvtable <- as.matrix(dietasvtable)
dietasvtable <- t(dietasvtable)
dietasvtable <- as.data.frame(dietasvtable)
diettaxatable <- read.csv('DietTaxa.csv', sep=',', header=TRUE, row.names=1, check.names=F)

#Feces from the intestine, n = 9 per treatment
gutsampleinfo <- read.csv('MetadataSeYe.csv', sep=',', header=TRUE, row.names=1, check.names=F)
gutsampleinfo
gutasvtable <- read.csv('ASV.csv', sep=',', header=TRUE, row.names=1, check.names=F)
gutasvtable <- as.matrix(gutasvtable)
gutasvtable <- t(gutasvtable)
gutasvtable <- as.data.frame(gutasvtable)
guttaxatable <- read.csv('Taxa.csv', sep=',', header=TRUE, row.names=1, check.names=F)

#env_table_16S <- read.csv('Growth-data.csv', sep=',', header=TRUE, row.names=1, check.names=F)

#Class command allows us to determin e the type of our inputted data. 
class(sampleinfo)
class(asvtable)
class(taxatable)
#class(env_table_16S)

#Knowing that all three are data.frames, we can make a microtable. All downstream analysis can be applied to the diet and gut tables. 
mydatadiet <- microtable$new(sample_table = dietsampleinfo, otu_table = dietasvtable, tax_table = diettaxatable)
mydatagut <- microtable$new(sample_table = gutsampleinfo, otu_table = gutasvtable, tax_table = guttaxatable)


#Now, we clean out dataset
#Command: filter_pollution takes out non-bacterial DNA sources (ie, mitochondria, chloroplasts)
mydatagut$filter_pollution(taxa = c("mitochondria", "chloroplast", "eukaryota", "bacteria_uncl_uncl_uncl_uncl_uncl_uncl", "bacteria_uncl_uncl_uncl_uncl_uncl"))
mydatadiet$filter_pollution(taxa = c("mitochondria", "chloroplast", "eukaryota", "bacteria_uncl_uncl_uncl_uncl_uncl_uncl", "bacteria_uncl_uncl_uncl_uncl_uncl"))

#Command: tidy_dataset trims the dataset further to unify the taxonomic information
mydata$tidy_dataset()

write.csv(dataset, file = "my_otu_table.csv", row.names = TRUE) # Set to TRUE if you want row names (sample IDs)


#Tidy taxonomy

#Next step is to check the sequencing depth of our samples - we base our decision on
#smallest depth, so that the sequence depth does not impact analysis
#Command: sample_sums checks the sequence depth of each sample
mydata$sample_sums() %>% range

#Now we can view our data: 
#Now all samples have an equal sequencing depth 

print(mydata)

All downstream analysis can be applied to the diet and gut tables. 

# from microtable to phyloseq object
dietphyseq <- meco2phyloseq(mydatadiet)
dietphyseq
gutphyseq <- meco2phyloseq(mydatagut)
gutphyseq


samplenames()
##################################
#Step 2: Create a folder for taxa abundance

#Command: cal_abund will calculate the taxa abundance of our data
mydatagut$cal_abund()
#The result is stored in object$taxa_abund ...
class(mydatagut$taxa_abund)
#Classified as a 'list' datatype
#From here, we can check relative abundance at taxonomic levels 
#Phylum 
mydatagut$taxa_abund$Species

#Saves the relative abundance for each phylogenetic level in a folder 'taxa_abund'
mydatagut$save_abund(dirpath = "taxa_abund")

##################################
#Step 3: Calculating alpha diversity 

#Command: cal_alphadiv will provide overall diversity for each sample
#with no statistics. 
#Use PD = (Faith's phylogenetic diversity), T = true if added phylogenetic tree to dataset
#Otherwise, like here, F = false is what we use. 
mydatagut$cal_alphadiv(PD = FALSE)
#The result is stored in object$alpha_diversity ...
class(mydatagut$alpha_diversity)
#Stored as a dataframe data type

#Saves the alpha diversity for each treatment in a folder called 'alpha_diversity'
mydatagut$save_alphadiv(dirpath = "alpha_diversity")

#Now, we create a trans_alpha object - two return data.frames will be given, 
#data_alpha: used to follow differential tests and plots
#data_stat

mydatagut$sample_table

t1 <- trans_alpha$new(dataset = mydatagut, group = "Diet")
t1$data_stat[1:5,] #Shows 5 of the 18 lines
t1$data_stat #Shows all lines

#To test differences among two groups, we will use a non-parametric t-test Wilcoxon - Rank Sum Tests
t1$cal_diff(method = "KW")
#The result is stored in object$res_diff ...

t1$res_diff #this shows the overview of all alpha diversity statistics
#now you can plot a specific index that is significant to see which group is higher
par(mfrow=c(2,2))

t1$plot_alpha(measure = "Simpson")
t1$plot_alpha(measure = "InvSimpson")
t1$plot_alpha(measure = "Shannon")
t1$plot_alpha(measure = "Chao1")
t1$plot_alpha(measure = "Pielou")
t1$plot_alpha(measure = "Fisher")
t1$plot_alpha(measure = "Observed")
t1$plot_alpha(measure = "Coverage")

tab <-microbiome::alpha(physeq, index = "all")
table(head(tab))


p.shannon <- boxplot_alpha(physeq, 
                           index = "shannon",
                           x_var = "Diet")
p.shannon <- p.shannon + theme_minimal() +  
  labs(x="\nTreatment Group", y="Shannon diversity\n") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.shannon

p.chao <- boxplot_alpha(physeq, 
                        index = "chao1",
                        x_var = "Diet")
p.chao <- p.chao + theme_minimal() + 
  labs(x="\nTreatment Group", y="Chao1 richness\n") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.chao

p.obs <- boxplot_alpha(physeq, 
                       index = "observed",
                       x_var = "Diet")
p.obs <- p.obs + theme_minimal() + 
  labs(x="\nTreatment Group", y="Observed ASVs\n") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.obs

grid.arrange(p.obs, p.shannon, p.chao, ncol=1)   


##################################
#Step 4: Calculating beta diversity 

#Command: cal_betadiv will provide overall diversity for each sample
#with no statistics. 
#unifrac = FALSE means to not calculate the unifrac method, because again, 
#there is no phylogenetic data/fasta in the dataset.

mydatadiet$cal_betadiv(unifrac = FALSE)
mydatagut$cal_betadiv(unifrac = FALSE)

#As seen, the result is stored in a file called object$beta_diversity
#We can call upon it using:
mydatagut$beta_diversity
class(mydatagut$beta_diversity)
#Classified as a list data type.
#We can now save the data to a directory. 
mydatagut$save_betadiv(dirpath = "beta_diversity")

#In order to understand the structure of the data, we create a trans_beta object
#measuring the parameter can invoke the distance matrix in dataset$beta_diversity

t1 <- trans_beta$new(dataset = mydatagut, group = "Diet", measure = "bray")
t1$cal_manova(manova_all = FALSE) #To see pairwise, set manova_all = FALSE
t1$res_manova

t2 <- trans_beta$new(dataset = mydatagut, group = "Diet", measure = "bray")
t2$cal_manova(manova_all = FALSE)
t2$res_manova


# Assuming t1 is your trans_beta object
# Define your own colors
my_colors <- c(Day0="#e3342f", InCON="#f6993f", InY1="#ffed4a", InY2="#38c172", InY3="#4dc0b5", OrCON="#3490dc", OrY1="#6574cd", OrY2="#9561e2",OrY3="#f66d9b")
my_colorsdiet <- c(InCON="#f6993f", InY1="#ffed4a", InY2="#38c172", InY3="#4dc0b5", OrCON="#3490dc", OrY1="#6574cd", OrY2="#9561e2",OrY3="#f66d9b")

#PCoA plot for Bray-Curtis.
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
plot.gut <- t1$plot_ordination(plot_type = c("point", "chull"), plot_color = "Diet", color_values = my_colors)
plot.diet 
plot.gut


grid.arrange(plot.diet, plot.gut, ncol=2)   


#PCoA plot for Jaccard
t2$cal_ordination(method = "PCoA")
class(t2$res_ordination)
# plot the PCoA result with confidence ellipse
t2$plot_ordination(plot_color = "Diet", plot_type = c("point", "chull"))

##################################
#Step 5: Determine what taxonomic composition is present in these fish. 

#The easiest way to present microbiome data is as a percentage, referred to as
#the relative abundance for the sample composition. 
#In short, we cannot know exact concentrations of these microbes - but we how
#much there is relative to each other. 

#Here, we'll use the microbiome package to plot microbial composition. 
#Papers typically report the phylum and genus levels. 

#First, we'll report thr 10 phyla with highest abundance in the dataset. 

p1_plot <- trans_abund$new(dataset = mydatagut, taxrank = "Phylum", ntaxa = 10)
p1_plot$plot_bar(others_color = "grey70", facet = "Diet", xtext_keep = FALSE, legend_text_italic = FALSE)


#Using groupmean parameter, we can attain the groupmean barplot.
#What this does is combine each of the individual samples into one mean barplot. 

p2_plot <- trans_abund$new(dataset = mydatagut, taxrank = "Phylum", ntaxa = 10, groupmean = "Diet")
phy2_plot <- p2_plot$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
phy2_plot + theme_classic() + theme(axis.title.y = element_text(size = 18))

#Can also be done for all taxanomic ranks, just change taxrank = "Family", "Class", etc. as needed. Typically phylum + genus are used in publications.

#Genera plots
g1_plot <- trans_abund$new(dataset = mydatagut, taxrank = "Genus", ntaxa = 20)
g1_plot$plot_bar(others_color = "grey70", facet = "Diet", xtext_keep = FALSE, legend_text_italic = FALSE, guide_legend(ncol = 1))

#Groupmean
g2_plot <- trans_abund$new(dataset = mydatagut, taxrank = "Genus", ntaxa = 20, groupmean = "Diet", )
gen2_plot <- g2_plot$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
gen2_plot + theme_classic() + theme(axis.title.y = element_text(size = 18)) + guides(color = guide_legend(ncol = 1))

##################################
#Step 6: Determine which specific microbes are enriched in gut microbiome when
#fed different diets. This is done via a differential abundance approach. 

#We will test statistically, which microbial taxa are significantly different 
#in between the treatments. 

#Metastat analysis at the Genus level, which depends on the permutations and t-test.
#Used for comparisons of taxonomic abundance between any paired groups, 
#regardless of taxonomic level. 

#For Phylum: 
p1 <- trans_diff$new(dataset = dataset_rarefied, method = "metastat", group
                     = "Diet", taxa_level = "Phylum")
select <- c("*","**","***","****") #we select for significance
select
p1$res_diff %<>% subset(Significance %in% select)
p1$plot_diff_abund(add_sig = T, add_sig_label = "Significance")
#No significant difference in relative abundance between treatments at Phylum level. 

#For Genus:  
g1 <- trans_diff$new(dataset = dataset_rarefied, method = "metastat", group
                     = "Diet", taxa_level = "Genus")
select <- c("*","**","***","****") #we select for significance
select
g1$res_diff %<>% subset(Significance %in% select)
g1$plot_diff_abund(add_sig = T, add_sig_label = "Significance")


##################################
#Step 7: Correlation of microbiome data and the performance of host. 
#Often done using redundancy analysis as a statistical approach, predicting 
#how the microbial community is affected by certain parameters. 

#First, we'll use add_data to incorporate the growth data to our frame. 
t1 <- trans_env$new(dataset = dataset_rarefied, add_data = env_table_16S)

#Then, we'll utilize Wilcoxon Rank Sum Test as an example to check if the data
#differs between the two dietary groups. 
t1$cal_diff(group = "Diet", method = "wilcox")
head(t1$res_diff)
t1 <- trans_env$new(dataset = dataset_rarefied, add_data = env_table_16S)

#Use Bray-Curtis distance for dbRDA
t1$cal_ordination(method = "RDA", use_measure = "bray")
#t1$res_rda is the result list stored in the object
t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1)
#t1$res_rda_trans is the transformed result for plotting
t1$plot_ordination(plot_color = "Diet", plot_type = c("point", "ellipse"))

#It shows if your data are significantly dissimilar
t1$cal_ordination_anova()
t1$res_ordination_axis
#It shows you if any of the growth correlates with the differences
t1$cal_ordination_envfit()
t1$res_ordination_envfit

##################################
#Step 8: Heat Plots and Other Plots


library(microbiome) # Load libraries
library(phyloseq)
library(dplyr)
library(reshape2)
library(knitr)

p <- microbiome::transform(physeq, "compositional") %>% 
  plot_composition(plot.type = "heatmap",
                   sample.sort = "neatmap", 
                   otu.sort = "neatmap") +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size = 9, hjust=1),
        legend.text = element_text(size = 8)) +
  ylab("Samples")
print(p)



# show 15 taxa at Class level
t1 <- trans_abund$new(dataset = dataset_rarefied, taxrank = "Genus", ntaxa = 15)
t1$plot_box(group = "Diet", xtext_angle = 30)

t2 <- trans_abund$new(dataset = dataset_rarefied, taxrank = "Genus", ntaxa = 15)
t2$plot_box(group = "Diet", xtext_angle = 30)

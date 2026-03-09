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


if (!requireNamespace("BiocManager", quietly = TRUE, force = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("microbiome")
install.packages("file2meco", repos = BiocManager::repositories())

install.packages("rcompanion")

library(file2meco)
library(phyloseq)
library(microbiome)
library(rcompanion)

install.packages("forcats")
library(forcats)

#Import the data
sampleinfo <- read.csv('D0MetadataNO128B106B.csv', sep=',', header=TRUE, row.names=1, check.names=F)
sampleinfo$Diet <- fct_relevel(sampleinfo$Diet, "Day0", "CON", "BAX", "MAO", "MAP", "SAX")
levels(sampleinfo$Diet)


asvtable <- read.csv('D0ASVNO128B106B.csv', sep=',', header=TRUE, row.names=1, check.names=F)
asvtable <- as.matrix(asvtable)
asvtable <- t(asvtable)
asvtable <- as.data.frame(asvtable)
taxatable <- read.csv('D0Taxa.csv', sep=',', header=TRUE, row.names=1, check.names=F)

 #env_table_16S <- read.csv('Growth-data.csv', sep=',', header=TRUE, row.names=1, check.names=F)

#Class command allows us to determin e the type of our inputted data. 
class(sampleinfo)
class(asvtable)
class(taxatable)
#class(env_table_16S)

#Knowing that all three are data.frames, we can make a microtable. 
mydata <- microtable$new(sample_table = sampleinfo, otu_table = asvtable, tax_table = taxatable)

#Now, we clean out dataset
#Command: filter_pollution takes out non-bacterial DNA sources (ie, mitochondria, chloroplasts)
mydata$filter_pollution(taxa = c("mitochondria", "chloroplast", "eukaryota", "bacteria_uncl_uncl_uncl_uncl_uncl_uncl", "bacteria_uncl_uncl_uncl_uncl_uncl"))
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

# from microtable to phyloseq object
data("mydata")
physeq <- meco2phyloseq(mydata)
physeq



##################################
#Step 2: Create a folder for taxa abundance

#Command: cal_abund will calculate the taxa abundance of our data
mydata$cal_abund()
#The result is stored in object$taxa_abund ...
class(mydata$taxa_abund)
#Classified as a 'list' datatype
#From here, we can check relative abundance at taxonomic levels 
#Phylum 
mydata$taxa_abund$Phylum

#Saves the relative abundance for each phylogenetic level in a folder 'taxa_abund'
mydata$save_abund(dirpath = "taxa_abund")

##################################
#Step 3: Calculating alpha diversity 

#Command: cal_alphadiv will provide overall diversity for each sample
#with no statistics. 
#Use PD = (Faith's phylogenetic diversity), T = true if added phylogenetic tree to dataset
#Otherwise, like here, F = false is what we use. 
mydata$cal_alphadiv(PD = FALSE)
#The result is stored in object$alpha_diversity ...
class(mydata$alpha_diversity)

#Stored as a dataframe data type

#Saves the alpha diversity for each treatment in a folder called 'alpha_diversity'
mydata$save_alphadiv(dirpath = "alpha_diversity")

#Now, we create a trans_alpha object - two return data.frames will be given, 
#data_alpha: used to follow differential tests and plots
#data_stat

mydata$sample_table

t1 <- trans_alpha$new(dataset = mydata, group = "Diet")
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
                           x_var = "Diet",
                           fill.colors = c(CON="#41AB8C", BAX="#DD782C", MAO="#8B86BD", MAP="#E84C9B", SAX="#7EB343"))

p.shannon <- p.shannon + theme_minimal() +  
  labs(x="\nTreatment Group", y="Shannon diversity\n") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.shannon

p.chao <- boxplot_alpha(physeq, 
                        index = "chao1",
                        x_var = "Diet",
                        fill.colors = c(CON="#41AB8C", BAX="#DD782C", MAO="#8B86BD", MAP="#E84C9B", SAX="#7EB343"))

p.chao <- p.chao + theme_minimal() + 
  labs(x="\nTreatment Group", y="Chao1 index\n") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.chao

p.obs <- boxplot_alpha(physeq, 
                       index = "observed",
                       x_var = "Diet",
                       fill.colors = c(CON="#41AB8C", BAX="#DD782C", MAO="#8B86BD", MAP="#E84C9B", SAX="#7EB343"))

p.obs <- p.obs + theme_minimal() + 
  labs(x="\nTreatment Group", y="Observed ASVs\n") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.obs


##################################
#Step 4: Calculating beta diversity 

#Command: cal_betadiv will provide overall diversity for each sample
#with no statistics. 
#unifrac = FALSE means to not calculate the unifrac method, because again, 
#there is no phylogenetic data/fasta in the dataset.

mydata$cal_betadiv(unifrac = FALSE)
#As seen, the result is stored in a file called object$beta_diversity
#We can call upon it using:
mydata$beta_diversity
class(mydata$beta_diversity)
#Classified as a list data type.
#We can now save the data to a directory. 
mydata$save_betadiv(dirpath = "beta_diversity")

#In order to understand the structure of the data, we create a trans_beta object
#measuring the parameter can invoke the distance matrix in dataset$beta_diversity

t1 <- trans_beta$new(dataset = mydata, group = "Diet", measure = "bray")
t1$cal_manova(manova_all = TRUE) #To see pairwise, set manova_all = FALSE
t1$res_manova

t2 <- trans_beta$new(dataset = mydata, group = "Diet", measure = "jaccard")
t2$cal_manova(manova_all = TRUE)
t2$res_manova

my_colors <- c(Day0="#e3342f", CON="#41AB8C", BAX="#DD782C", MAO="#8B86BD", MAP="#E84C9B", SAX="#7EB343")

#PCoA plot for Bray-Curtis.
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
plot.d0 <- t1$plot_ordination(plot_type = c("point", "chull"), plot_color = "Diet", color_values = my_colors)
plot.d0
#PCoA plot for Jaccard
t2$cal_ordination(method = "PCoA")
class(t2$res_ordination)
# plot the PCoA result with confidence ellipse
t2$plot_ordination(plot_color = "Diet", plot_type = c("point", "chull"))

p.diet

betaplots <- grid.arrange(p.diet, plot.d0, ncol=2)

##################################
#Step 5: Determine what taxonomic composition is present in these fish. 

#The easiest way to present microbiome data is as a percentage, referred to as
#the relative abundance for the sample composition. 
#In short, we cannot know exact concentrations of these microbes - but we how
#much there is relative to each other. 

#Here, we'll use the microbiome package to plot microbial composition. 
#Papers typically report the phylum and genus levels. 

#First, we'll report thr 10 phyla with highest abundance in the dataset. 

p1_plot <- trans_abund$new(dataset = mydata, taxrank = "Phylum", ntaxa = 10)
p1_plot$plot_bar(others_color = "grey70", facet = "Diet", xtext_keep = FALSE, legend_text_italic = FALSE)

p2_plot <- trans_abund$new(dataset = dataset_rarefied, taxrank = "Phylum", ntaxa = 10)
p2_plot$plot_bar(others_color = "grey70", facet = "SxD", xtext_keep = FALSE, legend_text_italic = FALSE)


#Using groupmean parameter, we can attain the groupmean barplot.
#What this does is combine each of the individual samples into one mean barplot. 

p2_plot <- trans_abund$new(dataset = dataset_rarefied, taxrank = "Phylum", ntaxa = 10, groupmean = "Treatment")
phy2_plot <- p2_plot$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
phy2_plot + theme_classic() + theme(axis.title.y = element_text(size = 18))

#Family plots
f1_plot <- trans_abund$new(dataset = dataset_rarefied, taxrank = "Family", ntaxa = 20)
f1_plot$plot_bar(others_color = "grey70", facet = "NTRMT", xtext_keep = FALSE, legend_text_italic = FALSE)

f2_plot <- trans_abund$new(dataset = dataset_rarefied, taxrank = "Family", ntaxa = 20)
f2_plot$plot_bar(others_color = "grey70", facet = "STRMT", xtext_keep = FALSE, legend_text_italic = FALSE)


#Groupmean
f2_plot <- trans_abund$new(dataset = mydata, taxrank = "Family", ntaxa = 20, groupmean = "STRMT")
fam2_plot <- f2_plot$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
fam2_plot + theme_classic() + theme(axis.title.y = element_text(size = 18))

#Genera plots
g1_plot <- trans_abund$new(dataset = mydata, taxrank = "Genus", ntaxa = 20)
g1_plot$plot_bar(others_color = "grey70", facet = "Diet", xtext_keep = FALSE, legend_text_italic = FALSE)

g2_plot <- trans_abund$new(dataset = dataset_rarefied, taxrank = "Genus", ntaxa = 20)
g2_plot$plot_bar(others_color = "grey70", facet = "SxD", xtext_keep = FALSE, legend_text_italic = FALSE)

#Groupmean
g2_plot <- trans_abund$new(dataset = dataset_rarefied, taxrank = "Genus", ntaxa = 20, groupmean = "Treatment")
gen2_plot <- g2_plot$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
gen2_plot + theme_classic() + theme(axis.title.y = element_text(size = 18))

heatmap <- trans_abund$new(dataset = dataset_rarefied, taxrank = "Species", ntaxa = 10, groupmean = "STRMT")
heatmap_plot <- heatmap$plot_heatmap(
  color_values = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")),
  facet = NULL,
  facet_switch = "y",
  x_axis_name = NULL,
  order_x = NULL,
  withmargin = TRUE,
  plot_numbers = FALSE,
  plot_text_size = 4,
  plot_breaks = NULL,
  margincolor = "white",
  plot_colorscale = "log10",
  min_abundance = 0.01,
  max_abundance = NULL,
  strip_text = 11,
  xtext_keep = TRUE,
  xtext_angle = 0,
  xtext_size = 10,
  ytext_size = 11,
  xtitle_keep = TRUE,
  grid_clean = TRUE,
  legend_title = "% Relative\nAbundance",
  pheatmap = FALSE)
heatmap_plot

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
                     = "STRMT", taxa_level = "Phylum")
select <- c("*","**","***","****") #we select for significance
select
p1$res_diff %<>% subset(Significance %in% select)
p1$plot_diff_abund(add_sig = T, add_sig_label = "Significance")
#No significant difference in relative abundance between treatments at Phylum level. 

#For Family:
f1 <- trans_diff$new(dataset = dataset_rarefied, method = "metastat", group
                     = "STRMT", taxa_level = "Family")
select <- c("*","**","***","****") #we select for significance
select
f1$res_diff %<>% subset(Significance %in% select)
f1$plot_diff_abund(add_sig = T, add_sig_label = "Significance")


#For Genus:  
g1 <- trans_diff$new(dataset = dataset_rarefied, method = "metastat", group
                     = "STRMT", taxa_level = "Genus")
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
t1$plot_box(group = "NTRMT", xtext_angle = 30)

t2 <- trans_abund$new(dataset = dataset_rarefied, taxrank = "Genus", ntaxa = 15)
t2$plot_box(group = "STRMT", xtext_angle = 30)

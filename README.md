Hi internet (and the future me who will use this in order to understand my own code for projects down the line). 

This repository is my slightly messy code base that demonstrates a basic workflow for demultiplexed fastq files off of Illumina NextSeq, sequencing the V3-V5 region of bacteria found in the feces of rainbow trout. 

STEP 1)
Use DADA2.R to take your demultiplexed fastq files and transform them into an ASV table and taxa table. 

STEP 2)
Create a metadata file that is relevant to your data, in the format provided in the .CSV files attached.

STEP 3) 
Input these 3 file types into your AstaxanthinStats.R workflow, which uses microeco and microbiome packages in order to perform statistics on the alpha diversity and beta diversity, as well as make publication ready plots
for both the aforementioned diversity indices and relative abundance stacked barplots. 

STEP 4) 
In order to do individual statistics on the abundances of key bacteria, you will open the taxa_abund file created in your working directory
from STEP 3's workflow. Open the _abund file inside of this at the taxonomic rank you are most interested in (typically phylum, genus, and species)
then take the top 10 to do individual stats on. Copy, transpose into columns and add next to your metadata, in order to properly model using a Kruskal Wallis 
test (one x variable), or even aligned rank transformation (if multiple x variables). Examples of this can be found in X...

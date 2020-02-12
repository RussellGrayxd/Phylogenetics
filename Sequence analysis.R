setwd("")
# install necessary packges 
install.packages("remotes")
remotes::install_github("GuangchuangYu/treeio")
install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("DECIPHER")
install.packages("ape")

library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet("Amanita.txt", format = "fasta")

# look at some of the sequences (optional)
seqs

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
BrowseSeqs(aligned, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned,
                file="Amanita_aligned.fasta")

# read in the aligned data
dna <- read.alignment("Amanita_aligned.fasta", format = "fasta")

# create a distance matrix for the alignment 
D <- dist.alignment(dna, matrix = "similarity")


temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)+
  scale_color_viridis()#darker shades of gray mean a larger distance # you can also make cool color plots but they're much more complicated because they use the image() function

# we can start to see a pattern because the data is ordered by year, 
# but we can't really make any conclusions yet

tre <- nj(D)
class(tre) #all trees created using {ape} package will be of class phylo

tre <- ladderize(tre)

# ~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~Base R plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot(tre, cex = 0.6)
title("similarity in Amanita (ITS)")


# or 
h_cluster <- hclust(D, method = "average", members = NULL) # method = average is used for UPGMA, members can be equal to NULL or a vector with a length of size D
plot(h_cluster, cex = 0.6)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tree Plotting in ggtree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
# you can fan it out 
ggtree(tre, yscale = "NA")+
  geom_tiplab(hjust = -0.3, size=4, align = TRUE)+
  xlim(0,0.5) 

# or whatever this thing does???
ggtree(tre,layout = "daylight")+
  geom_tiplab(hjust = -0.3, size=4, align = TRUE)+
  xlim(0,0.5) 

# plot a basic tree
ggtree(tre) + 
  geom_tiplab(hjust = -0.3, size=4, align = TRUE)+
  xlim(0,0.5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Customize your trees ~~~~~~~~~~~~~~~~~~~~~~~~
# plot using ggtree and highlight clusters
# change the node values for your own data
ggtree(tre) + 
  geom_tiplab(hjust = -0.3, size=4, align = TRUE) + 
  geom_hilight(node=19, fill="purple", alpha = 0.2) + 
  geom_hilight(node=17, fill="dark green", alpha = 0.2) +
  geom_hilight(node=20, fill="gold", alpha = 0.2) +
  xlim(0,0.5) 

# highlight clusters and add a vertical line to group clusters
# change the node values for your own data
ggtree(tre) + 
  geom_tiplab(hjust = -0.3, size=4, align = TRUE) + 
  geom_hilight(node=19, fill="purple", alpha = 0.2) + 
  geom_hilight(node=17, fill="dark green", alpha = 0.2) +
  geom_hilight(node=20, fill="gold", alpha = 0.2) +
  geom_cladelabel(node=19, label=" Cluster 1", 
                  color="purple", offset=.1, barsize = 2,
                  fontsize = 5, align=TRUE, alpha = 0.5) + 
  geom_cladelabel(node=17, label=" Cluster 2", 
                  color="dark green", offset=.1, barsize = 2,
                  fontsize = 5, align=TRUE, alpha = 0.5) + 
  geom_cladelabel(node=20, label=" Cluster 3", 
                  color="gold", offset=.1, barsize = 2,
                  fontsize = 5, align=TRUE, alpha = 0.5) + 
  xlim(0,0.5) 
  


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Plot the allignment with the tree ~~~~~~~~~~~~~~~~~

# lets plot the alignment with the tree, to do this we first have to
# match the names to the tip labels
# set our tree into a new name
tre.new <- tre
# change tip labels to full alignment names
tre.new$tip.label <- aligned@ranges@NAMES

# plot the alignment 
msaplot(p=ggtree(tre.new), fasta="Amanita_aligned.fasta", 
        window=c(150, 175))+
  scale_fill_viridis_d(alpha = 0.8)








#Note: This no longer works due to updates to R

#remove.packages("phytools")
#install.packages("phytools", dependencies = TRUE)
library(phytools)

# Load individual trees
tree1 <- read.tree("1.txt")
tree2 <- read.tree("2.txt")
tree3 <- read.tree("3.txt")

# Combine them into a list
my_trees <- c(tree1, tree2, tree3)
# Generate the combined tree
combined_tree <- phytools::superTree(my_trees, method="mincut")

# Plot the result to see how it looks
plotTree(combined_tree)
# Save as a Newick file
write.tree(combined_tree, file="my_combined_supertree.tre")

# Install and load phangorn
if (!require("phangorn")) install.packages("phangorn")
if (!require("ape")) install.packages("ape")
library(phangorn)
library(ape)

# Load individual trees
tree1 <- read.tree("1.tre")
tree2 <- read.tree("2.tre")
tree3 <- read.tree("3.tre")

# Combine them into a list
my_trees <- c(tree1, tree2, tree3)

# Get a master list of every unique taxon across all trees
all_taxa <- unique(unlist(lapply(my_trees, function(x) x$tip.label)))

# Function to convert each tree into a binary matrix with missing data (?)
# This is the manual version of what 'mrpData' does
make_mrp_part <- function(tr, all_taxa) {
  # Get the splits (partitions) for this specific tree
  bp <- prop.part(tr)
  # Convert splits to a matrix (Rows = Taxa, Cols = Splits)
  mat <- matrix("?", length(all_taxa), length(bp), 
                dimnames = list(all_taxa, NULL))
  
  for(i in 1:length(bp)) {
    present <- tr$tip.label[bp[[i]]]
    absent <- setdiff(tr$tip.label, present)
    mat[present, i] <- "1"
    mat[absent, i] <- "0"
  }
  return(mat)
}

# Combine all trees into one large Supermatrix
full_matrix <- do.call(cbind, lapply(my_trees, make_mrp_part, all_taxa = all_taxa))

# Convert to a format R can use for tree building
mrp_data <- as.phyDat(full_matrix, type = "USER", levels = c("0", "1"))

# Find the best Supertree (The Parsimony Ratchet)
combined_results <- pratchet(mrp_data, trace = 0)

# Final Consensus and Plot
# Experiment with p value as this show how strict the tree is (p=1 is most strict)
final_tree <- consensus(combined_results, p = 0.5)
plot(final_tree, main="Manual MRP Supertree")

# Saving Tree
write.tree(final_tree, file = "combined_supertree.tre")


# Viewing Tree
#Use FigTree to visualize the resulting .tre file
#Root the tree inside of FigTree and choose decending node order
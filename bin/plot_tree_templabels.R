# Load required libraries
library(ggtree)
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)

# Function to plot the phylogenetic tree with temperatures and strain labels, then save to a PNG file
plot_tree_with_temps <- function(treefile, temp_file, output_file) {
  # Load the tree file (returns an S3 phylo object)
  tree <- treeio::read.tree(treefile)
  
  # Load the temperature data
  temp_data <- read_tsv(temp_file)
  
  # Rename the strain column in the temp_data to avoid conflict with the tree's 'label'
  temp_data <- temp_data %>% dplyr::rename(strain_label = strain)
  
  # Extract the tree tip labels
  tree_tip_labels <- tree$tip.label
  
  # Perform partial matching of strain names to tree tip labels
  temp_data$tree_label <- sapply(temp_data$strain_label, function(strain) {
    matching_label <- grep(strain, tree_tip_labels, value = TRUE)
    
    # If a match is found, return the matching tree tip label, otherwise NA
    if (length(matching_label) > 0) {
      return(matching_label[1])  # Take the first match if multiple
    } else {
      return(NA)
    }
  })
  
  # Remove rows where no match was found (optional based on data completeness)
  temp_data <- temp_data %>% filter(!is.na(tree_label))
  
  # Merge the temperature data with the tree's tip labels based on 'tree_label'
  merged_data <- data.frame(label = tree$tip.label) %>%
    left_join(temp_data, by = c("label" = "tree_label"))
  
  # Plot the tree with temperatures mapped to the tips and strain labels added
  p <- ggtree(tree) %<+% merged_data +
    geom_tippoint(aes(color = temp), size = 3) +  # Add color-coded points based on temperature
    geom_tiplab(aes(label = strain_label), align = TRUE, linetype = "dotted", size = 3) +  # Add strain labels
    scale_color_viridis_c(option = "plasma") +  # Color scale for temperature
    theme_tree2()  # Tree styling
  
  # Save the plot to a PNG file
  ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
  
  # Optionally, also display the plot
  print(p)
}

# Example usage:
# plot_tree_with_temps("path/to/cktree.treefile", "path/to/temps.tsv", "tree_plot_with_labels.png")

# Example usage:
plot_tree_with_temps("/home/sheidig/bioenvada/rtrial/CK_0000128b_filtered_NT_checked.afasta.treefile", "/home/sheidig/bioenvada/rtrial/temp_info.tsv", "/home/sheidig/bioenvada/rtrial/tree_pl.png")

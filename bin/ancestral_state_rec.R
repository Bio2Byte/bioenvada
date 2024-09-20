# Load required libraries
library(ggtree)
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(ape)

# Function to label missing internal nodes, perform ancestral state reconstruction, and plot the tree
plot_tree_with_reconstructed_temps <- function(treefile, temp_file, output_file) {
  # Load the tree file (returns an S3 phylo object)
  tree <- treeio::read.tree(treefile)
  
  # Ensure all internal nodes have labels. If any are missing, assign unique labels like NodeA, NodeB, etc.
  if (is.null(tree$node.label) || any(tree$node.label == "")) {
    missing_node_labels <- sum(is.na(tree$node.label) | tree$node.label == "")
    
    # Assign distinct labels like NodeA, NodeB, etc. to missing internal nodes
    tree$node.label[is.na(tree$node.label) | tree$node.label == ""] <- 
      paste0("Node", LETTERS[1:missing_node_labels])
  }

  print(tree$node.label)
  
  # Load the temperature data
  temp_data <- read_tsv(temp_file)
  
  # Rename the strain column in the temp_data to avoid conflict with the tree's 'label'
  temp_data <- temp_data %>% dplyr::rename(strain_label = strain)
  
  # Extract the tree tip labels
  tree_tip_labels <- tree$tip.label
  
  # Perform partial matching of strain names to tree tip labels
  temp_data$tree_label <- sapply(temp_data$strain_label, function(strain) {
    matching_label <- grep(strain, tree_tip_labels, value = TRUE)

    print(matching_label)
    # If a match is found, return the matching tree tip label, otherwise NA
    if (length(matching_label) > 0) {
      return(matching_label[1])  # Take the first match if multiple
    } else {
      return(NA)
    }
  })
  
  print(temp_data)
  # Remove rows where no match was found
  temp_data <- temp_data %>% filter(!is.na(tree_label))
  
  print(temp_data)
  # Create a vector of temperatures for the tree's tips, initialized with NA
  temp_vector <- rep(NA, length(tree$tip.label))
  
  # Assign temperatures to the matching tips
  for (i in seq_along(tree$tip.label)) {
    match_idx <- match(tree$tip.label[i], temp_data$tree_label)
    if (!is.na(match_idx)) {
      temp_vector[i] <- temp_data$temp[match_idx]
    }
  }
  
  # Perform ancestral state reconstruction for the internal nodes using ace()
  # A Brownian motion model is used, but you can change it based on your data
  reconstructed_states <- ace(temp_vector, tree, method = "REML")
  

  print(reconstructed_states)
  # Combine the known tip temperatures and reconstructed internal node temperatures
  # We use the correct numbering for internal nodes from the tree object
  full_temp_data <- data.frame(label = c(tree$tip.label, tree$node.label),
                               temp = c(temp_vector, reconstructed_states$ace))
  

  print(full_temp_data)
  # Create the plot
  p <- ggtree(tree) %<+% full_temp_data +
    geom_tippoint(aes(color = temp), size = 3) +  # Add color-coded points based on known/reconstructed temperature
    geom_nodepoint(aes(color = temp), size = 3) +  # Show temperatures for internal nodes
    geom_tiplab(aes(label = label), size = 3, offset = 0.02 ) +  # Add strain labels with offset
    geom_nodelab(aes(label = label), size = 2, hjust = 1.3, vjust = 1.5) +  # Show labels for internal nodes
    scale_color_viridis_c(option = "plasma") +  # Color scale for temperature
    theme_tree2()  +# Tree styling
    coord_cartesian(clip = "off")+
    theme(
           # axis.text.x  = element_text(size = 12),
            plot.margin = margin(, 4, , , "cm"),
            legend.position="top"
    )
  
  # Save the plot to a PNG file
  ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
  
}


# Example usage:
plot_tree_with_reconstructed_temps("/home/sheidig/bioenvada/rtrial/CK_0000128b_filtered_NT_checked.afasta.rooted.treefile", "/home/sheidig/bioenvada/rtrial/temp_info.tsv", "/home/sheidig/bioenvada/rtrial/tree_anc.png")

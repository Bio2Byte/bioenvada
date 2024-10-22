# Load required libraries
library(ggtree)
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(ape)
library(RColorBrewer)

#conda install conda-forge::r-ggplot2 conda-forge::r-readr conda-forge::r-dplyr conda-forge::r-stringr bioconda::bioconductor-ggtree  conda-forge::r-ape  conda-forge::r-scales

print("package import done")


# Function to map input data: tree and temperature data
map_input <- function(treefile, temp_file) {
  # Load the tree file
  tree <- treeio::read.tree(treefile)
  
  # Ensure all internal nodes have labels
  if (is.null(tree$node.label) || any(tree$node.label == "")) {
    missing_node_labels <- sum(is.na(tree$node.label) | tree$node.label == "")
    tree$node.label[is.na(tree$node.label) | tree$node.label == ""] <- 
      paste0("Node", LETTERS[1:missing_node_labels])
  }
  
  # Load the temperature data
  temp_data <- read_tsv(temp_file) %>%
    dplyr::rename(strain_label = strain)
  
  # Partial matching of strain names to tree tip labels
  tree_tip_labels <- tree$tip.label
  temp_data$tree_label <- sapply(temp_data$strain_label, function(strain) {
    matching_label <- grep(strain, tree_tip_labels, value = TRUE)
    if (length(matching_label) > 0) {
      return(matching_label[1])  # Return first match
    } else {
      return(NA)
    }
  })
  
  # Filter out unmatched strains
  temp_data <- temp_data %>% filter(!is.na(tree_label))
  
  return(list(tree = tree, temp_data = temp_data))
}

# Function to reconstruct ancestral states
reconstruct_ancestral_states <- function(tree, temp_data) {
  # Create a temperature vector for tree tips
  temp_vector <- rep(NA, length(tree$tip.label))
  # Map temperatures to tree tips
  for (i in seq_along(tree$tip.label)) {
    print(tree$tip.label[i])
    match_idx <- match(tree$tip.label[i], temp_data$tree_label)
    print(match_idx)
    print(temp_data$temp[match_idx])
    if (!is.na(match_idx)) {
      temp_vector[i] <- temp_data$temp[match_idx]
    }
  }
  
  print(temp_vector)
  # Ancestral state reconstruction using ace()
  reconstructed_states <- ace(temp_vector, tree, type='continuous',  method = "pic")
  
  # Combine tip and node temperatures
  full_temp_data <- data.frame(
    label = c(tree$tip.label, tree$node.label),
    temp = c(temp_vector, reconstructed_states$ace)
  )
  
  return(full_temp_data)
}

# Function to create and save the plot
create_plot <- function(tree, full_temp_data, output_file) {
  p <- ggtree(tree) %<+% full_temp_data +
    geom_tippoint(aes(color = temp), size = 3) +
    geom_nodepoint(aes(color = temp), size = 3) +
    geom_tiplab(aes(label = label), size = 3, offset = 0.02) +
    geom_nodelab(aes(label = label), size = 2, hjust = 1.3, vjust = -1) +
    scale_color_viridis_c(option = "plasma") +
    theme_tree2() +
    coord_cartesian(clip = "off") +
    theme(
      plot.margin = margin(,4, , , "cm"),
      legend.position = "top"
    )
  
  # Save the plot
  ggsave(paste(output_file, ".png", sep = ""), plot = p, width = 10, height = 8, dpi = 300)
}


# Main function to orchestrate the workflow
plot_tree_with_reconstructed_temps <- function(treefile, temp_file, output_file) {
  input_data <- map_input(treefile, temp_file)
  
  full_temp_data <- reconstruct_ancestral_states(input_data$tree, input_data$temp_data)
  write.csv(full_temp_data, paste(output_file, "_rec_temps.csv", sep = ""), row.names = FALSE)

  create_plot(input_data$tree, full_temp_data, output_file)
}

# Example usage:
# plot_tree_with_reconstructed_temps("treefile_path", "temp_file_path", "output_file_path")
#plot_tree_with_reconstructed_temps("/home/sheidig/bioenvada/rtrial/CK_0000128b_filtered_NT_checked.afasta.rooted.treefile", "/home/sheidig/bioenvada/rtrial/temp_info.tsv", "/home/sheidig/bioenvada/rtrial/tree_anc")


# Add command line argument handling
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript ancestral_state_rec.R <tree_file> <temp_file> <output_name>")
}

tree_file <- args[1]
temp_file <- args[2]
output_file <- args[3]

# Call the main function with command line arguments
plot_tree_with_reconstructed_temps(tree_file, temp_file, output_file)
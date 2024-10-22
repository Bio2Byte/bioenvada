# Load required libraries
library(ape)
library(ggplot2)
library(dplyr)


# Function to prepare data: tree and trait data
prepare_data <- function(treefile, trait_file) {
  phylo_tree <- read.tree(file = treefile)

  trait_data <- read.table(trait_file, header = TRUE, sep = ",")
  print (trait_data)

  new_trait_data <- trait_data %>% select(species, everything())
  print(new_trait_data)
  return(list(phylo_tree = phylo_tree, new_trait_data = new_trait_data))
}

# Function to perform PIC
perform_pic <- function(trait_data, phylo_tree) {
  pic_temperature <- pic(trait_data$temp, phylo_tree)
  print(pic_temperature)

  results <- list()
  correlation_data <- data.frame(behavior = character(), correlation = numeric(), stringsAsFactors = FALSE)
  
  for (behavior in colnames(trait_data)) {
    if (startsWith(behavior, 'R')){
      
      pic_behavior <- pic(trait_data[[behavior]], phylo_tree)
      correlation <- cor(pic_temperature, pic_behavior)
      
      results[[behavior]] <- list(pic_temperature = pic_temperature, pic_behavior = pic_behavior, correlation = correlation)
      correlation_data <- rbind(correlation_data, data.frame(behavior = behavior, correlation = correlation))
  }}
  
  return(list(results = results, correlation_data = correlation_data))
}

# Function to plot correlations
plot_correlations <- function(correlation_data, output_file) {

  correlation_data$behavior <- factor(correlation_data$behavior, levels = correlation_data$behavior)
  

  every_nth <- function(n) {
    function(x) {
      idx <- seq(1, length(x), by = n)
      labels <- rep("", length(x))
      labels[idx] <- x[idx]
      return(labels)
    }
  }
  
  plot <- ggplot(correlation_data, aes(x = behavior, y = correlation)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Correlation between temperature and biophysical propensity per residue",
         x = "Residue",
         y = "Correlation") +
    scale_x_discrete(labels = every_nth(10)) +
    theme(
      panel.grid.major.x = element_blank(),  # Remove vertical grid lines
      panel.grid.minor = element_blank(),      # Remove minor grid lines
      axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
    scale_y_continuous(limits = c(-1, 1),      # Set y-axis limits from -1 to 1
                       breaks = seq(-1, 1, by = 0.5))   #
  
  ggsave(output_file, plot = plot, width = 10, height = 8, dpi = 300)


}

# Main function to run the workflow
analyze_tree_and_traits <- function(treefile, trait_file, correlation_output_file, plot_output_file) {
  input_data <- prepare_data(treefile, trait_file)
  pic_results <- perform_pic(input_data$new_trait_data, input_data$phylo_tree)
  
  write.csv(pic_results$correlation_data, correlation_output_file, row.names = FALSE)
  plot_correlations(pic_results$correlation_data, plot_output_file)
}

# Add command line argument handling
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script_name.R <treefile> <trait_file>")
}



treefile <- args[1]
trait_file <- args[2]
correlation_output_file <- paste('pic', trait_file, sep = "")
plot_output_file <- paste('pic',trait_file,'.png', sep = "")


# Call the main function with command line arguments
analyze_tree_and_traits(treefile, trait_file, correlation_output_file, plot_output_file)

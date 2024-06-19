Input is a dataframe in which each column is a different fine-mapping method in the format (LD panel.window.method) and each column contains the fine-mapped sig SNPs (pip > 0.5 and parts of a 95% credible set).

All code below is in R.

The order of columns can differ either to assess a given method & window across LD panels, or different fine-mapping methods within a given LD ref panel & window.

### 1. COMBINED HEATMAP FOR 32 FINEMAPPING APPROACHES

```
library(dplyr)
library(readr)
library(ggplot2)

# Read SNP data from CSV file
snp_data <- read_csv("~/Desktop/files_for_figures/all_SNPs_finemap_noduplis.csv")

# Calculate Jaccard index
n_methods <- ncol(snp_data)
jaccard_matrix <- matrix(NA, nrow = n_methods, ncol = n_methods)

for (i in 1:n_methods) {
  for (j in 1:n_methods) {
    snp_set1 <- snp_data[[i]]
    snp_set2 <- snp_data[[j]]
    intersection <- length(intersect(snp_set1, snp_set2))
    union_set <- length(union(snp_set1, snp_set2))
    jaccard_matrix[i, j] <- intersection / union_set
  }
}

# Convert Jaccard matrix to data frame
jaccard_df <- as.data.frame(jaccard_matrix)
colnames(jaccard_df) <- colnames(snp_data)
rownames(jaccard_df) <- colnames(snp_data)

# Create lower triangle matrix
lower_triangle <- function(mat) {
  mat[upper.tri(mat)] <- NA
  return(mat)
}

jaccard_df_new <- lower_triangle(jaccard_df)

# Plot heatmap using ggplot2
heatmap <- ggplot(data = expand.grid(Method1 = colnames(jaccard_df_new), Method2 = colnames(jaccard_df_new)),
                  aes(x = Method1, y = Method2, fill = jaccard_df_new[cbind(match(Method1, colnames(jaccard_df_new)), match(Method2, colnames(jaccard_df_new)))])) +
  geom_tile() +
  # scale_fill_viridis_c(name = "Jaccard Index", na.value = "white") + 
  scale_fill_distiller(palette = "RdBu", direction = - 1, name = "Jaccard Index", na.value = "white", limits = c(0, 1)) +
  labs(x = "", y = "", title = "Jaccard Index of Finemapped SNPs Overlap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), 
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 10),  
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",             
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10)) + 
  coord_fixed() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())  # remove gridlines and background

print(heatmap)
```


### 2. CALCULATE MEAN JACCARD INDEX ACROSS ALL ANALYSES

```
library(dplyr)
library(readr)

# Read SNP data from CSV file
snp_data <- read_csv("~/Desktop/files_for_figures/all_SNPs_finemap_noduplis.csv")
##we can also filter to exclude NoLD approaches and use this df for the rest of analyses
#filtered_snp_data <- snp_data %>% select(-starts_with("NoLD"))

# Calculate Jaccard index
n_methods <- ncol(snp_data)
jaccard_matrix <- matrix(NA, nrow = n_methods, ncol = n_methods)

for (i in 1:n_methods) {
  for (j in 1:n_methods) {
    snp_set1 <- snp_data[[i]]
    snp_set2 <- snp_data[[j]]
    intersection <- length(intersect(snp_set1, snp_set2))
    union_set <- length(union(snp_set1, snp_set2))
    jaccard_matrix[i, j] <- intersection / union_set
  }
}

# Calculate the mean Jaccard index excluding the diagonal
jaccard_values <- jaccard_matrix[lower.tri(jaccard_matrix, diag = FALSE)]
mean_jaccard <- mean(jaccard_values, na.rm = TRUE)

print(paste("Mean Jaccard Index:", mean_jaccard))
#"Mean Jaccard Index: 0.479589396958867"

```

### 3. CALCULATE MEAN JACCARD INDEX PER EACH FINEMAPPING APPROACH

```
library(dplyr)
library(readr)

# Read SNP data from CSV file
snp_data <- read_csv("~/Desktop/files_for_figures/all_SNPs_finemap_noduplis.csv")

# Function to calculate Jaccard index matrix for a subset of methods
calculate_jaccard <- function(data_subset) {
  n_methods <- ncol(data_subset)
  jaccard_matrix <- matrix(NA, nrow = n_methods, ncol = n_methods)
  
  for (i in 1:n_methods) {
    for (j in 1:n_methods) {
      snp_set1 <- data_subset[[i]]
      snp_set2 <- data_subset[[j]]
      intersection <- length(intersect(snp_set1, snp_set2))
      union_set <- length(union(snp_set1, snp_set2))
      jaccard_matrix[i, j] <- intersection / union_set
    }
  }
  
  return(jaccard_matrix)
}

# Function to calculate mean Jaccard index excluding the diagonal
calculate_mean_jaccard <- function(jaccard_matrix) {
  jaccard_values <- jaccard_matrix[lower.tri(jaccard_matrix, diag = FALSE)]
  mean_jaccard <- mean(jaccard_values, na.rm = TRUE)
  return(mean_jaccard)
}

# Loop through SNP data in subsets of 4 columns and calculate mean Jaccard index for each subset
num_columns <- ncol(snp_data)
subset_size <- 4
subset_indices <- seq(1, num_columns, by = subset_size)

for (start_col in subset_indices) {
  end_col <- min(start_col + subset_size - 1, num_columns)
  data_subset <- snp_data[, start_col:end_col]
  method_names <- colnames(data_subset)
  
  jaccard_matrix <- calculate_jaccard(data_subset)
  mean_jaccard <- calculate_mean_jaccard(jaccard_matrix)
  
  print(paste("Mean Jaccard Index for methods", paste(method_names, collapse = ", "), ":", mean_jaccard))
}
```

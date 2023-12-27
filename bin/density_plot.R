# Function to print help and usage instructions
print_help <- function() {
  cat("
Usage: Rscript density_plotter.R <filename> [chrom_string] [merge]

<filename>      : Path to the TE superfamily density table generated using density_table.py.
[chrom_string]  : Optional. Specifies a string pattern to include in the analysis. Chromosomes without this string will be excluded. Useful for excluding scaffolds, contigs, or specific chromosome IDs.
[merge]         : Optional. Boolean flag. If 'merge' is specified, all chromosome plots are merged onto a single page.

This script plots the distribution of EDTA annotations across chromosomes.

Dependencies: ggplot2, dplyr for plotting. R and its packages can be installed with conda.
For installation:
conda config --add channels conda-forge && conda config --add channels CRAN && conda create -n r_env r-base r-ggplot2 r-dplyr

")
  quit(status = 0)
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for no arguments or help request
if (length(args) == 0 || args[1] == '-h' || args[1] == '--help') {
  print_help()
}

filename <- args[1]

# Check if the file has a header
has_header <- readLines(filename, n = 1)
header <- if (grepl("chrom", has_header, fixed = TRUE)) TRUE else FALSE

repeats <- read.table(filename, sep='\t', header = header, 
                      col.names = if (!header) c('chrom', 'coord', 'density', 'type') else NULL)

# Exclude specific types
excluded_types <- c('target_site_duplication', 'repeat_region', 'long_terminal_repeat')
repeats <- subset(repeats, !type %in% excluded_types)

# Replace 'Copia' with 'Ty1' and 'Gypsy' with 'Ty3'
library(dplyr)
repeats <- repeats %>%
    mutate(type = gsub("Copia_LTR_retrotransposon", "Ty1_LTR_retrotransposon", type)) %>%
    mutate(type = gsub("Gypsy_LTR_retrotransposon", "Ty3_LTR_retrotransposon", type))

# Reorder the type factor to move 'Ty1_LTR_retrotransposon' higher up
repeats$type <- factor(repeats$type, levels = c("Ty1_LTR_retrotransposon", setdiff(unique(repeats$type), "Ty1_LTR_retrotransposon")))

# Initialize variables
chrom_string <- NULL
merge_plots <- FALSE

# Check arguments
for (arg in args[-1]) {
  if (tolower(arg) == "merge") {
    merge_plots <- TRUE
  } else {
    chrom_string <- arg
  }
}

# Optional chromosome string filtering
if (!is.null(chrom_string)) {
  repeats <- subset(repeats, grepl(chrom_string, chrom))
}

# Convert chromosome names to a factor with the correct order
chrom_levels <- unique(repeats$chrom)
chrom_levels_sorted <- chrom_levels[order(as.numeric(gsub("[^0-9]", "", chrom_levels)))]
repeats$chrom <- factor(repeats$chrom, levels = chrom_levels_sorted)

# Plotting
library(ggplot2)

# Define a shuffled color palette using hcl.colors with the 'Dark 3' palette
num_categories <- length(unique(repeats$type))
set.seed(46) # Change the seed for a new shuffle
color_palette <- sample(hcl.colors(num_categories, "Dark 3"))

plot_border_theme <- theme(
  panel.border = element_rect(colour = "grey", fill=NA, linewidth=1),
  axis.text.x = element_text(),
  axis.ticks.x = element_line(),
  plot.margin = margin(1, 1, 1, 1, "lines") # Adjust plot margin
)

if (merge_plots) {
  # Merged plot logic
  p <- ggplot(repeats, aes(x = coord / 1e6, y = density * 100, color = type)) + 
       geom_line() +
       facet_wrap(~ chrom, scales = "free_x") + 
       labs(title = "Density of Repeat Types on Chromosomes (%)", 
            x = "Position on chromosome (Mb)", 
            y = "Percent repeat type in window (%)", 
            color = "Type") + 
       theme_minimal() +
       plot_border_theme +
       scale_color_manual(values = color_palette)

  # Save the plot to a PDF file
  pdf("chromosome_density_plots_merged.pdf", width = 8, height = 6)
  print(p)
  dev.off()
} else {
  # Unmerged plot logic
  pdf("chromosome_density_plots.pdf", width = 8, height = 6)
  for (chr in unique(repeats$chrom)) {
    subset_df <- repeats[repeats$chrom == chr,]
    
    if (nrow(subset_df) > 1) {
      p <- ggplot(subset_df, aes(x = coord / 1e6, y = density * 100, group = interaction(chrom, type), color = type)) + 
           geom_line() + 
           labs(title = paste("Density of Repeat Types on", chr, "(%)"), 
                x = "Position on chromosome (Mb)", 
                y = "Percent repeat type in window (%)", 
                color = "Type") +
           theme_minimal() +
           plot_border_theme +
           scale_color_manual(values = color_palette) +
           scale_x_continuous(expand = c(0, 0)) + # Tighten grey plot border for x-axis
           scale_y_continuous(expand = expansion(mult = c(0.005, 0.01))) 
      print(p)
    } else {
      cat(paste("Not enough data points to plot for chromosome:", chr, "\n"))
    }
  }
  dev.off()
}

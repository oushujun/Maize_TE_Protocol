# Function to print help and usage instructions
print_help <- function() {
  cat("
Usage: Rscript density_plotter.R <filename> [chrom_prefix] [merge]

<filename>      : Path to the TE superfamily density table generated using density_table.py.
[chrom_prefix]  : Optional. Specifies the chromosome prefix to include in the analysis. Chromosomes without this prefix will be excluded. Useful for excluding scaffolds & contigs.
[merge]         : Optional. Boolean flag. If 'merge' is specified, all chromosome plots are merged onto a single page.

This script plots the distribution of EDTA annotations across chromosomes.

Dependencies: ggplot2 for plotting. R and its packages can be installed with conda.
For installation:
conda config --add channels conda-forge && conda config --add channels CRAN && conda create -n r_env r-base r-ggplot2

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

# Optional chromosome filtering
if (length(args) > 1 && args[2] != "merge") {
  chrom_prefix <- args[2]
  repeats <- subset(repeats, grepl(paste0("^", chrom_prefix), chrom))
}

# Convert chromosome names to a factor with the correct order
chrom_levels <- unique(repeats$chrom)
chrom_levels_sorted <- chrom_levels[order(as.numeric(gsub("chr", "", chrom_levels)))]
repeats$chrom <- factor(repeats$chrom, levels = chrom_levels_sorted)

# Check for 'merge' flag
merge_plots <- FALSE
if (length(args) >= 3 && tolower(args[3]) == "merge") {
  merge_plots <- TRUE
}

# Plotting
library(ggplot2)

if (merge_plots) {
  p <- ggplot(repeats, aes(x = coord / 1e6, y = density * 100, color = type)) + 
       geom_line() +
       facet_wrap(~ chrom, scales = "free_x") +  # Scales panels by chromosome length
       labs(title = "Density of Repeat Types on Chromosomes (%)", 
     		x = "Position on chromosome (Mb)", 
    		y = "Percent repeat type in window (%)", 
     color = "Type") + 
       theme_minimal() +
       theme(axis.text.x = element_blank(),  # Remove x-axis text
             axis.ticks.x = element_blank())  # Remove x-axis ticks

  # Save the plot to a PDF file
  pdf("chromosome_density_plots_merged.pdf", width = 8, height = 6)
  print(p)
  dev.off()
} else {
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
           theme_minimal()
      print(p)
    } else {
      cat(paste("Not enough data points to plot for chromosome:", chr, "\n"))
    }
  }
  dev.off()
}

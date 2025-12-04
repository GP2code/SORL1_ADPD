library(tidyverse)

# Set working directory
block_folder <- "~/Desktop/${WORK_DIR}"

# List of ancestries 
ancestries <- c("CAH", "FIN", "EAS", "CAS", "SAS", 
                "AMR", "AJ", "MDE", "AFR", "AAC", "EUR")

# Load block files
load_block_file <- function(ancestry) {
  # lowercase ancestry for filename
  ancestry_lower <- tolower(ancestry)
  file_path <- file.path(block_folder, paste0(ancestry_lower, "_case_blocks.blocks.det"))
  read_table(file_path, col_names = TRUE) %>%
    mutate(Ancestry = ancestry)
}

# Combine all ancestry files
all_blocks <- map_dfr(ancestries, load_block_file)
all_blocks$Ancestry <- factor(all_blocks$Ancestry, levels = ancestries)

# Define variant positions and 2-color grouping
variant_lines <- tibble(
  Variant_Pos = c(
    121564878, 121585352,                # Red group
    121496917, 121496918, 121522975,     # Blue group
    121590137, 121605213,
    121614877, 121614890,
    121618854, 121621073,
    121559627
  ),
  Group = c(rep("Red", 2), rep("Blue", 10))
)

# Define labels with correct groups
variant_labels <- tibble(
  Variant_Pos = c(
    121564878, 121585352,                # Red group
    121496917,                           # Blue group
    121522975, 121590137, 121605213,
    121614877,
    121618854, 121621073,
    121559627
  ),
  Label = c(
    "11:121564878", "11:121585352",
    "11:121496917,11:121496918",
    "11:121522975", "11:121590137", "11:121605213",
    "11:121614877,11:121614890",
    "11:121618854", "11:121621073",
    "11:121559627"
  ),
  Group = c(rep("Red", 2), rep("Blue", 8))
)

# Set custom colors
variant_colors <- c("Red" = "red", "Blue" = "blue")

# Plot
plot <- ggplot(all_blocks, aes(x = BP1, xend = BP2, y = Ancestry, yend = Ancestry)) +
  geom_segment(linewidth = 4, color = "black") +
  geom_vline(data = variant_lines, aes(xintercept = Variant_Pos, color = Group), linetype = "dashed") +
  geom_text(
    data = variant_labels,
    mapping = aes(x = Variant_Pos, y = length(ancestries) + 1.5, label = Label, color = Group),
    inherit.aes = FALSE,
    size = 1.5,
    angle = 90,
    vjust = 0,
    hjust = 0
  ) +
  scale_color_manual(values = variant_colors) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 8) +
  labs(
    title = "Haplotype Blocks in SORL1 Region (±100kb) Across Ancestries (Cases)",
    x = "Genomic Position (GRCh38)",
    y = "Ancestry"
  ) +
  theme(
    plot.title = element_text(size = 7.8),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.position = "none",
    plot.margin = margin(t = 80, r = 10, b = 10, l = 10)
  )

# Show plot
print(plot)

# Save plot
ggsave("~/Desktop/haplotype_blocks_colored_cases_variants.png",
       plot, width = 10, height = 6, bg = "white")

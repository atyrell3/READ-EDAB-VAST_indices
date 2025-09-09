# --- 2. Define and Save Predator Lists ---
# This section identifies the key predators (piscivores) and prey species
# based on diet overlap and observation counts.

# Load the diet overlap matrix.
dietoverlap <- read_csv(here("fhdat/tgmat.2022-02-15.csv"))

# Perform hierarchical clustering to identify predator guilds.
# This code block is a good candidate for a dedicated function if it were to be reused.
# The logic uses `hclust` to group predators by diet similarity.
d_dietoverlap <- dist(dietoverlap)
guilds <- hclust(d_dietoverlap, method = "complete")
dend <- as.dendrogram(guilds)
dend <- dendextend::color_branches(dend, k = 6)

labels(dend) <- paste(
  as.character(names(dietoverlap[-1]))[order.dendrogram(dend)],
  "(",
  labels(dend),
  ")",
  sep = ""
)

# Extract the list of piscivorous predators from the clustering results.
# The `partition_leaves` function is a dendextend utility to get nodes from the tree.
# This is a key step, as it determines which fish are considered "piscivores" in the analysis.
pisccomplete <- dendextend::partition_leaves(dend)[[
  dendextend::which_node(
    dend,
    c("Bluefish..S(37)", "Bluefish..M(36)", "Bluefish..L(35)")
  )
]]

# Create a data frame of the selected piscivores for joining.
# `str_remove` and `str_extract` are used to parse the species and size category from the names.
pisccompletedf <- data.frame(
  "COMNAME" = toupper(str_remove(pisccomplete, "\\..*")),
  "SizeCat" = str_remove(str_extract(pisccomplete, "\\..*[:upper:]+"), "\\.."),
  "feedguild" = "pisccomplete"
)

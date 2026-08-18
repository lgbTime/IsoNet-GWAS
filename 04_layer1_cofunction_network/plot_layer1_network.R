#!/usr/bin/env Rscript
#
# plot_network.R — Publication-quality network visualization for GWAS-to-function networks
#
# Reads nodes.tsv and edges.tsv produced by gwas2network.py and renders:
#   1. A full concentric circular network plot
#   2. A focused plot of the top-50 most significant genes + their functions
#
# Usage:
#   Rscript plot_network.R chlorophyll_ab_nodes.tsv chlorophyll_ab_edges.tsv chlorophyll_ab
#

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforce)
  library(ggrepel)
  library(igraph)
  library(RColorBrewer)
  library(viridisLite)
})

# ---- Command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_network.R <nodes.tsv> <edges.tsv> [prefix]")
}
nodes_file <- args[1]
edges_file <- args[2]
prefix     <- ifelse(length(args) >= 3, args[3], tools::file_path_sans_ext(nodes_file))

# ---- Which views to generate (environment variable, comma-separated) ----
#   full          = concentric: phenotype → genes → all annotations
#   focused       = top-50 genes + annotation neighbors (often cluttered)
#   cofunction    = co-function modules (shared terms)
#   orphans       = orphan genes only (no co-function partners)
#   cofunction_orphans = cofunction modules + orphan outer ring
#   Default: full,cofunction,cofunction_orphans (skip focused + orphans-only)
views_env <- Sys.getenv("PLOT_VIEWS")
if (nchar(views_env) == 0) views_env <- "full,cofunction,cofunction_orphans"
views_wanted <- trimws(unlist(strsplit(views_env, ",")))
cat(sprintf("  Views: %s\n", paste(views_wanted, collapse=", ")))

# ---- Read data ----
cat(sprintf("Reading: %s\n", nodes_file))
nodes <- fread(nodes_file)
cat(sprintf("Reading: %s\n", edges_file))
edges <- fread(edges_file)

cat(sprintf("  Nodes: %d  |  Edges: %d\n", nrow(nodes), nrow(edges)))

# ---- Extract phenotype info from data ----
pheno_row <- nodes[type == "phenotype"]
pheno_name     <- pheno_row$id[1]
pheno_display  <- pheno_row$display_name[1]
if (is.na(pheno_display) || pheno_display == "") pheno_display <- gsub("_", " ", pheno_name)
cat(sprintf("  Phenotype: %s\n", pheno_name))

# ---- Prepare data ----
# Classify node categories for coloring
# Use fifelse chain (compatible with data.table >= 1.14)
nodes[, node_type_cat := fifelse(type == "phenotype", "Phenotype",
  fifelse(type %in% c("gene", "transcript"), "Gene / Transcript",
  fifelse(type == "COG_category", "COG Category",
  fifelse(type == "description", "Functional Description",
  fifelse(type == "GO_term", "GO Term",
  fifelse(type == "KEGG_pathway", "KEGG Pathway",
  fifelse(type == "PFAM_domain", "PFAM Domain",
  "Other")))))))]

# Numeric significance for gene nodes
nodes[, sig_num := as.numeric(significance)]
nodes[is.na(sig_num), sig_num := 0]

# ---- Color palettes ----
# Qualitative palette for function categories
func_cats <- setdiff(unique(nodes$node_type_cat), c("Phenotype", "Gene / Transcript"))
n_func_cats <- length(func_cats)
# Use Set3 (up to 12 colors) or Paired
if (n_func_cats <= 6) {
  func_colors <- setNames(brewer.pal(max(3, n_func_cats), "Set2")[1:n_func_cats], func_cats)
} else {
  qual_cols <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Pastel2"), brewer.pal(8, "Accent"))
  func_colors <- setNames(qual_cols[1:n_func_cats], func_cats)
}

# ---- Build igraph object ----
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# ---- Compute custom concentric layout ----
compute_concentric <- function(g, nodes_dt) {
  n_all <- vcount(g)
  # Identify nodes
  pheno_idx <- which(V(g)$type == "phenotype")
  gene_idx  <- which(V(g)$type %in% c("gene", "transcript"))
  func_idx  <- which(!(V(g)$type %in% c("phenotype", "gene", "transcript")))

  n_genes <- length(gene_idx)
  n_funcs <- length(func_idx)

  # --- Adaptive radii based on gene count ---
  #    Base: for 300 genes → gene_r=3, func_r=6
  #    Fewer genes → larger rings to fill the plot better
  #    More genes → slightly smaller rings to avoid overcrowding
  scale_factor <- sqrt(300 / max(n_genes, 30))
  gene_r  <- 3 * scale_factor
  func_r  <- gene_r * 2

  cat(sprintf("      Adaptive radii: gene_r=%.1f  func_r=%.1f  (n_genes=%d)\n",
              gene_r, func_r, n_genes))

  lay <- matrix(0, nrow = n_all, ncol = 2)

  # Phenotype at center
  if (length(pheno_idx) > 0) {
    lay[pheno_idx, ] <- c(0, 0)
  }

  # Genes: circle at radius gene_r, ordered by significance (most sig at top/12 o'clock)
  if (n_genes > 0) {
    gene_order <- gene_idx[order(V(g)$sig_num[gene_idx], decreasing = TRUE)]
    angles <- seq(0, 2 * pi, length.out = n_genes + 1)[1:n_genes]
    angles <- pi/2 - angles
    lay[gene_order, 1] <- gene_r * cos(angles)
    lay[gene_order, 2] <- gene_r * sin(angles)
  }

  # Functions: circle at radius func_r, grouped by category
  if (n_funcs > 0) {
    func_types <- V(g)$node_type_cat[func_idx]
    for (ft in unique(func_types)) {
      ft_idx <- func_idx[func_types == ft]
      n_ft <- length(ft_idx)
      angle_start <- runif(1, 0, 2 * pi)
      ft_angles <- seq(angle_start, angle_start + 2 * pi, length.out = n_ft + 1)[1:n_ft]
      lay[ft_idx, 1] <- func_r * cos(ft_angles)
      lay[ft_idx, 2] <- func_r * sin(ft_angles)
    }
  }

  return(list(layout = lay, gene_r = gene_r, func_r = func_r))
}

cat("Computing concentric layout...\n")
set.seed(42)
layout_result <- compute_concentric(g, nodes)
layout_xy <- layout_result$layout
gene_radius  <- layout_result$gene_r
func_radius  <- layout_result$func_r

# ---- Build plot data frame ----
plot_nodes <- as.data.table(layout_xy)
setnames(plot_nodes, c("x", "y"))
plot_nodes[, name := V(g)$name]
plot_nodes[, type := V(g)$type]
plot_nodes[, node_type_cat := V(g)$node_type_cat]
plot_nodes[, display_name := V(g)$display_name]
plot_nodes[, sig_num := V(g)$sig_num]
plot_nodes[, description := V(g)$description]
plot_nodes[, COG_full_name := V(g)$COG_full_name]

# Label only phenotype + top 25 most significant genes + key functions
gene_nodes <- plot_nodes[type %in% c("gene", "transcript")]
top_genes <- gene_nodes[order(-sig_num)][1:25, name]
# Also label function nodes that connect to top genes
top_gene_ids <- top_genes
func_edges <- edges[source %in% top_gene_ids | target %in% top_gene_ids]
top_funcs <- unique(c(func_edges$source, func_edges$target))
top_funcs <- setdiff(top_funcs, c(top_genes, V(g)$name[V(g)$type == "phenotype"]))
# Limit function labels
if (length(top_funcs) > 40) top_funcs <- sample(top_funcs, 40)

label_nodes <- c(
  V(g)$name[V(g)$type == "phenotype"],
  top_genes,
  top_funcs
)
plot_nodes[, show_label := name %in% label_nodes]

# Clean display names for labels
plot_nodes[, label_text := display_name]
plot_nodes[type == "phenotype", label_text := paste0(pheno_display, "\n(Phenotype)")]
# Truncate long function descriptions
plot_nodes[type == "description" & nchar(label_text) > 50,
           label_text := paste0(substr(label_text, 1, 47), "...")]

# ---- Build edge plot data ----
# We sample edges for the full plot to reduce overplotting
# phenotype->gene edges: keep all
# gene->function edges: keep top genes only for clarity
pheno_name <- pheno_row$id[1]
all_gene_names <- V(g)$name[V(g)$type %in% c("gene", "transcript")]
top_gene_set <- head(all_gene_names[order(V(g)$sig_num[V(g)$type %in% c("gene", "transcript")],
                                          decreasing = TRUE)], 100)

plot_edges <- as.data.table(edges)
plot_edges[, keep := FALSE]
# Always keep phenotype edges
plot_edges[source == pheno_name | target == pheno_name, keep := TRUE]
# Keep annotation edges only for top 100 genes
plot_edges[source %in% top_gene_set | target %in% top_gene_set, keep := TRUE]

# Remove phenotype edges where both ends are pheno (shouldn't happen)
plot_edges <- plot_edges[keep == TRUE]

# Map to coordinates
plot_edges[, from_x := plot_nodes$x[match(source, plot_nodes$name)]]
plot_edges[, from_y := plot_nodes$y[match(source, plot_nodes$name)]]
plot_edges[, to_x   := plot_nodes$x[match(target, plot_nodes$name)]]
plot_edges[, to_y   := plot_nodes$y[match(target, plot_nodes$name)]]
plot_edges <- plot_edges[!is.na(from_x) & !is.na(to_x)]

cat(sprintf("  Plotting %d edges, %d nodes, %d labels\n",
            nrow(plot_edges), nrow(plot_nodes), sum(plot_nodes$show_label)))

# ---- COLORS ----
# Gene color: viridis by significance
gene_sig_range <- range(plot_nodes$sig_num[plot_nodes$type %in% c("gene", "transcript")], na.rm = TRUE)

# Node sizes
gene_scale <- scales::rescale(plot_nodes$sig_num,
  to = c(1.5, 5),
  from = gene_sig_range
)
plot_nodes[, node_size := fifelse(type == "phenotype", 8,
  fifelse(type %in% c("gene", "transcript"), gene_scale,
  fifelse(type == "COG_category", 3,
  fifelse(type == "GO_term", 1.5,
  fifelse(type == "KEGG_pathway", 2,
  fifelse(type == "description", 1.5,
  fifelse(type == "PFAM_domain", 2, 1.5)))))))]

# Edge alpha & size
plot_edges[, edge_alpha := fifelse(type == "association", 0.5, 0.08)]
plot_edges[, edge_size := fifelse(type == "association", 0.4, 0.1)]

# Attach function-category color to annotation edges
# Map each edge's target node to its node_type_cat in the full network
node_cat_map <- setNames(plot_nodes$node_type_cat, plot_nodes$name)
plot_edges[type == "annotation", func_cat := node_cat_map[target]]
# Fallback: if target not found (shouldn't happen), try source
plot_edges[type == "annotation" & (is.na(func_cat) | func_cat == ""),
           func_cat := node_cat_map[source]]
# Set unknown to "Other"
plot_edges[is.na(func_cat), func_cat := "Other"]
# Factor with levels matching the node color palette
plot_edges[, func_cat := factor(func_cat, levels = intersect(names(func_colors), unique(func_cat)))]
# If any levels remain unset, add them
missing_levels <- setdiff(unique(plot_edges$func_cat), levels(plot_edges$func_cat))
if (length(missing_levels) > 0) {
  plot_edges[, func_cat := factor(func_cat, levels = c(levels(func_cat), missing_levels))]
}

# ---- Plot 1: Full concentric network (academic quality) ----
cat("Plotting full network...\n")

# Labels: phenotype + top 35 genes by significance + hub functions (degree >= 2)
all_gene_ids <- plot_nodes[type %in% c("gene", "transcript"), name]
gene_nodes_plot <- plot_nodes[type %in% c("gene", "transcript")]
# Get top genes
gene_order_sig <- gene_nodes_plot[order(-sig_num)]
top_genes_full <- gene_order_sig[1:35, name]
# Also include all genes that have COG annotations (most informative)
genes_with_cog <- gene_nodes_plot[COG_full_name != "" & COG_full_name != "Function unknown", name]
top_genes_full <- unique(c(top_genes_full, head(genes_with_cog, 15)))

# Get function nodes that connect to top-labeled genes
top_gene_set_full <- unique(c(top_genes_full,
  head(gene_order_sig$name, 50)))  # top 50 for edge filtering
func_edges_for_labels <- edges[source %in% top_gene_set_full | target %in% top_gene_set_full]
func_label_candidates <- unique(c(func_edges_for_labels$source, func_edges_for_labels$target))
func_label_candidates <- setdiff(func_label_candidates, c(all_gene_ids, pheno_name))
# Count how many top genes each function connects to
func_conn_count <- sapply(func_label_candidates, function(fn) {
  sum((func_edges_for_labels$source == fn & func_edges_for_labels$target %in% top_gene_set_full) |
      (func_edges_for_labels$target == fn & func_edges_for_labels$source %in% top_gene_set_full))
})
top_func_labels <- names(sort(func_conn_count[func_conn_count >= 2], decreasing = TRUE))
if (length(top_func_labels) > 30) top_func_labels <- head(top_func_labels, 30)

label_nodes_full <- unique(c(
  pheno_name,
  top_genes_full,
  top_func_labels
))
plot_nodes[, show_label_full := name %in% label_nodes_full]
plot_nodes[, label_text_full := display_name]
plot_nodes[type == "phenotype", label_text_full := pheno_display]
# Truncate long function descriptions for readability
plot_nodes[type == "description" & nchar(label_text_full) > 45,
           label_text_full := paste0(substr(label_text_full, 1, 42), "...")]

# Build annotation edge color palette
edge_func_cats <- unique(plot_edges[type == "annotation"]$func_cat)
edge_func_cats <- edge_func_cats[!is.na(edge_func_cats)]
edge_colors <- func_colors[names(func_colors) %in% edge_func_cats]
if (length(edge_colors) == 0) edge_colors <- c("Other" = "grey50")

# Precompute all adaptive sizes
ring_label_size  <- scales::rescale(gene_radius, to = c(1.8, 3.2), from = c(2, 5))
gene_node_size   <- scales::rescale(gene_radius, to = c(2.2, 3.2), from = c(2, 5))
pheno_node_size  <- scales::rescale(gene_radius, to = c(10, 16), from = c(2, 5))
pheno_label_vjust <- -(scales::rescale(gene_radius, to = c(2, 3.2), from = c(2, 5)))
pheno_label_size  <- scales::rescale(gene_radius, to = c(3, 4.5), from = c(2, 5))
plot_limit        <- func_radius + 1.5

# Adaptive ring guides
ring_data <- rbind(
  data.frame(r = gene_radius,   label = "Gene / Transcript layer", stringsAsFactors = FALSE),
  data.frame(r = func_radius,   label = "Function annotation layer", stringsAsFactors = FALSE)
)
theta <- seq(0, 2*pi, length.out = 200)
ring_paths <- rbindlist(lapply(seq_len(nrow(ring_data)), function(i) {
  data.table(r = ring_data$r[i], theta = theta,
             x = ring_data$r[i] * cos(theta),
             y = ring_data$r[i] * sin(theta),
             ring_label = ring_data$label[i])
}))

ring_label_df <- data.frame(
  x = c(0, 0),
  y = c(gene_radius + 0.18, func_radius + 0.35),
  label = ring_data$label,
  stringsAsFactors = FALSE
)

p_full <- ggplot() +
  # 0. Concentric ring guides
  geom_path(
    data = ring_paths,
    aes(x = x, y = y, group = r),
    color = "grey80", linewidth = 0.3, linetype = "dotted", alpha = 0.7
  ) +
  # 0b. Ring labels (positioned just outside each ring)
  geom_text(
    data = ring_label_df,
    aes(x = x, y = y, label = label),
    size = ring_label_size,
    color = "grey50", fontface = "italic", hjust = 0.5
  ) +
  # 1. Annotation edges — subtle, function-category colored
  geom_segment(
    data = plot_edges[type == "annotation"],
    aes(x = from_x, y = from_y, xend = to_x, yend = to_y, color = func_cat),
    alpha = 0.16, linewidth = 0.14
  ) +
  # 2. Association edges — p-value encoded by linewidth (thickness) + alpha (opacity)
  #    Both scaled continuously: strongest hit = thickest + most opaque
  geom_segment(
    data = plot_edges[type == "association"],
    aes(x = from_x, y = from_y, xend = to_x, yend = to_y,
        linewidth = weight, alpha = weight),
    color = "#E65C00"    # deep orange-red, stands out from annotation edges
  ) +
  # 3. Function nodes (outer ring)
  geom_point(
    data = plot_nodes[!type %in% c("phenotype", "gene", "transcript")],
    aes(x = x, y = y, color = node_type_cat),
    alpha = 0.6, stroke = 0.05, size = 1.5
  ) +
  # 4. Function labels
  geom_text_repel(
    data = plot_nodes[show_label_full == TRUE & !type %in% c("phenotype", "gene", "transcript")],
    aes(x = x, y = y, label = label_text_full, color = node_type_cat),
    size = 1.5, max.overlaps = 60, box.padding = 0.2, point.padding = 0.1,
    segment.color = "grey80", segment.size = 0.15, fontface = "italic",
    show.legend = FALSE
  ) +
  # 5. Gene nodes (inner ring) — fill = -log10(p), size adapts to n_genes
  geom_point(
    data = plot_nodes[type %in% c("gene", "transcript")],
    aes(x = x, y = y, fill = sig_num),
    size = gene_node_size, shape = 21, color = "grey30", stroke = 0.15, alpha = 0.88
  ) +
  # --- Scales ---
  scale_color_manual(
    name = "Annotation type",
    values = func_colors,
    guide = guide_legend(order = 1, override.aes = list(size = 3.5, alpha = 0.85))
  ) +
  scale_fill_viridis_c(
    name = expression(-log[10](italic(p))~"(gene)"),
    option = "inferno", direction = -1,
    limits = gene_sig_range,
    guide = guide_colorbar(order = 2, barwidth = 0.7, barheight = 4,
                           title.position = "top")
  ) +
  scale_linewidth_continuous(
    name = expression(-log[10](italic(p))~"(link)"),
    range = c(0.15, 1.5),
    guide = guide_legend(order = 3)
  ) +
  scale_alpha_continuous(
    range = c(0.08, 0.65),
    guide = "none"
  ) +
  # 6. Phenotype node (center) — circle with gold fill + double ring
  geom_point(
    data = plot_nodes[type == "phenotype"],
    aes(x = x, y = y),
    size = pheno_node_size, shape = 21, fill = "#FFD700", color = "#B8860B", stroke = 2.2
  ) +
  # 6b. Inner accent ring on phenotype
  geom_point(
    data = plot_nodes[type == "phenotype"],
    aes(x = x, y = y),
    size = pheno_node_size * 0.58, shape = 21, fill = NA, color = "#8B6914", stroke = 0.6
  ) +
  # 7. Gene labels
  geom_text_repel(
    data = plot_nodes[show_label_full == TRUE & type %in% c("gene", "transcript")],
    aes(x = x, y = y, label = name),
    size = 2.0, color = "grey15", fontface = "bold",
    max.overlaps = 60, box.padding = 0.3, point.padding = 0.15,
    segment.color = "grey70", segment.size = 0.2,
    show.legend = FALSE
  ) +
  # 8. Phenotype label — positioned above the center node
  geom_label(
    data = plot_nodes[type == "phenotype"],
    aes(x = x, y = y, label = label_text_full),
    size = pheno_label_size,
    fontface = "bold", fill = "#FFF8DC", color = "#8B6914",
    label.padding = unit(0.4, "lines"),
    vjust = pheno_label_vjust,
    linewidth = 0.3
  ) +
  coord_fixed(
    xlim = c(-plot_limit, plot_limit),
    ylim = c(-plot_limit, plot_limit),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title       = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 9, hjust = 0.5, color = "grey35",
                                    margin = margin(t = 3, b = 5)),
    plot.margin      = margin(15, 15, 15, 15),
    legend.position  = "right",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.box       = "vertical",
    legend.margin    = margin(0, 0, 0, 0)
  ) +
  labs(
    title = paste0("GWAS Network Architecture: ", pheno_display),
    subtitle = sprintf(
      "%d significant genes/transcripts (inner ring)  |  %d functional annotation terms (outer ring)  |  %d annotation edges",
      length(unique(plot_nodes$name[plot_nodes$type %in% c("gene","transcript")])),
      sum(!plot_nodes$type %in% c("phenotype","gene","transcript")),
      nrow(plot_edges[type == "annotation"])
    )
  )

# Save full network
pdf_file <- paste0(prefix, "_network_full.pdf")
png_file <- paste0(prefix, "_network_full.png")

if ("full" %in% views_wanted) {
  ggsave(pdf_file, p_full, width = 18, height = 16, dpi = 150, device = "pdf")
  cat(sprintf("  -> %s\n", pdf_file))
  ggsave(png_file, p_full, width = 18, height = 16, dpi = 200, device = "png")
  cat(sprintf("  -> %s\n", png_file))
} else { cat("  skipped (full not in PLOT_VIEWS)\n") }


# ---- Plot 2: Focused top-50 network (more detail) ----
cat("Plotting focused top-50 network...\n")

top50_genes <- head(all_gene_names[order(V(g)$sig_num[V(g)$type %in% c("gene","transcript")],
                                          decreasing = TRUE)], 50)

# Subgraph: phenotype + top50 genes + their direct function neighbors
# neighborhood() returns vertex indices — convert to names
func_neighbor_idx <- unique(unlist(neighborhood(g, nodes = top50_genes, order = 1)))
func_neighbor_names <- V(g)$name[func_neighbor_idx]
sub_nodes <- unique(c(pheno_name, top50_genes, func_neighbor_names))
# Keep it manageable
if (length(sub_nodes) > 300) {
  # Prioritize functions that connect to multiple top genes
  func_in_sub <- setdiff(sub_nodes, c(pheno_name, top50_genes))
  edge_counts <- sapply(func_in_sub, function(fn) {
    sum((edges$source %in% top50_genes & edges$target == fn) |
          (edges$target %in% top50_genes & edges$source == fn))
  })
  func_in_sub <- names(sort(edge_counts, decreasing = TRUE))[1:min(200, length(func_in_sub))]
  sub_nodes <- c(pheno_name, top50_genes, func_in_sub)
}

sg <- induced_subgraph(g, sub_nodes)
cat(sprintf("  Subgraph: %d nodes, %d edges\n", vcount(sg), ecount(sg)))

# Compute layout: FR with phenotype pinned at center
set.seed(42)
sg_layout <- layout_with_fr(sg)
# Fix phenotype at (0,0)
pheno_sg_idx <- which(V(sg)$name == pheno_name)
if (length(pheno_sg_idx) > 0) {
  sg_layout[pheno_sg_idx, ] <- c(0, 0)
}

sg_nodes <- as.data.table(sg_layout)
setnames(sg_nodes, c("x", "y"))
sg_nodes[, name := V(sg)$name]
sg_nodes[, type := V(sg)$type]
sg_nodes[, node_type_cat := V(sg)$node_type_cat]
sg_nodes[, sig_num := V(sg)$sig_num]
sg_nodes[, description := V(sg)$description]
sg_nodes[, display_name := V(sg)$display_name]

# Scale layout
max_r <- max(sqrt(sg_nodes$x^2 + sg_nodes$y^2))
sg_nodes[, x := x / max_r * 5]
sg_nodes[, y := y / max_r * 5]
if (length(pheno_sg_idx) > 0) {
  sg_nodes[pheno_sg_idx, c("x", "y") := .(0, 0)]
}

# Build edges
sg_edges <- as.data.table(as_data_frame(sg, what = "edges"))
sg_edges[, from_x := sg_nodes$x[match(from, sg_nodes$name)]]
sg_edges[, from_y := sg_nodes$y[match(from, sg_nodes$name)]]
sg_edges[, to_x   := sg_nodes$x[match(to, sg_nodes$name)]]
sg_edges[, to_y   := sg_nodes$y[match(to, sg_nodes$name)]]
sg_edges <- sg_edges[!is.na(from_x) & !is.na(to_x)]

# Labels: phenotype + all top50 genes + functions connected to >=2 genes
func_deg <- table(c(sg_edges$from, sg_edges$to))
label50_names <- c(
  pheno_name,
  top50_genes,
  names(func_deg[names(func_deg) %in% setdiff(sub_nodes, c(pheno_name, top50_genes)) &
                  func_deg >= 2])
)
sg_nodes[, show_label := name %in% label50_names]
sg_nodes[, label_text := display_name]
sg_nodes[type == "phenotype", label_text := pheno_display]
sg_nodes[type == "description" & nchar(label_text) > 40,
         label_text := paste0(substr(label_text, 1, 37), "...")]

# Node sizes for subgraph
sg_gene_sig_range <- range(sg_nodes$sig_num[sg_nodes$type %in% c("gene", "transcript")], na.rm = TRUE)
if (diff(sg_gene_sig_range) == 0) sg_gene_sig_range <- c(0, 1)

sg_gene_scale <- scales::rescale(sg_nodes$sig_num,
  to = c(2.5, 7), from = sg_gene_sig_range)
sg_func_scale <- scales::rescale(
  as.numeric(func_deg[sg_nodes$name]),
  to = c(1.5, 4),
  from = c(1, max(func_deg, na.rm = TRUE))
)
sg_nodes[, node_size := fifelse(type == "phenotype", 10,
  fifelse(type %in% c("gene", "transcript"), sg_gene_scale,
  sg_func_scale))]
sg_nodes[is.na(node_size), node_size := 1.5]

# Colors
sg_edges[, edge_alpha := fifelse(type == "association", 0.6, 0.15)]
sg_edges[, edge_size := fifelse(type == "association", 0.5, 0.15)]

p_focus <- ggplot() +
  # Annotation edges — colored by function type
  geom_segment(
    data = sg_edges[type == "annotation"],
    aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
    alpha = 0.25, color = "grey50", linewidth = 0.2
  ) +
  # Association edges
  geom_segment(
    data = sg_edges[type == "association"],
    aes(x = from_x, y = from_y, xend = to_x, yend = to_y,
        alpha = weight),
    color = "#E69F00", linewidth = 0.7
  ) +
  # Function nodes
  geom_point(
    data = sg_nodes[!type %in% c("phenotype", "gene", "transcript")],
    aes(x = x, y = y, color = node_type_cat, size = node_size),
    alpha = 0.75, stroke = 0.15
  ) +
  # Gene nodes
  geom_point(
    data = sg_nodes[type %in% c("gene", "transcript")],
    aes(x = x, y = y, fill = sig_num, size = node_size),
    shape = 21, color = "grey30", stroke = 0.3, alpha = 0.9
  ) +
  # Phenotype
  geom_point(
    data = sg_nodes[type == "phenotype"],
    aes(x = x, y = y),
    size = 15, shape = 23, fill = "#FFD700", color = "#B8860B", stroke = 2
  ) +
  # Gene labels
  geom_text(
    data = sg_nodes[show_label == TRUE & type %in% c("gene", "transcript")],
    aes(x = x, y = y, label = name),
    size = 2.5, hjust = 0.5, vjust = -1.5, color = "grey20", fontface = "italic"
  ) +
  # Function labels
  geom_text(
    data = sg_nodes[show_label == TRUE & !type %in% c("phenotype", "gene", "transcript")],
    aes(x = x, y = y, label = label_text),
    size = 1.9, hjust = 0.5, vjust = -0.8, color = "grey40"
  ) +
  # Phenotype label
  geom_label(
    data = sg_nodes[type == "phenotype"],
    aes(x = x, y = y, label = label_text),
    size = 4.5, fontface = "bold", fill = "#FFF8DC", color = "#8B6914",
    label.padding = unit(0.5, "lines"), vjust = -2.5
  ) +
  scale_fill_viridis_c(
    name = "-log10(p)",
    option = "inferno", direction = -1,
    guide = guide_colorbar(order = 1, barwidth = 0.8, barheight = 5)
  ) +
  scale_color_manual(
    name = "Function Category",
    values = func_colors,
    guide = guide_legend(order = 2, override.aes = list(size = 3, alpha = 0.9))
  ) +
  scale_size_continuous(
    name = "Node importance",
    guide = guide_legend(order = 3)
  ) +
  scale_alpha_continuous(
    name = "GWAS\nassociation",
    guide = guide_legend(order = 4)
  ) +
  coord_fixed(clip = "off") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    plot.margin   = margin(20, 20, 20, 20),
    legend.position = "right",
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 8)
  ) +
  labs(
    title = "GWAS Network (Focused): Top 50 Significant Genes",
    subtitle = sprintf(
      "%s GWAS -- top 50 associations + functional annotations  |  %d nodes, %d edges",
      pheno_display,
      nrow(sg_nodes), nrow(sg_edges)
    )
  )

if ("focused" %in% views_wanted) {
  pdf_file2 <- paste0(prefix, "_network_focused.pdf")
  png_file2 <- paste0(prefix, "_network_focused.png")
  ggsave(pdf_file2, p_focus, width = 18, height = 16, dpi = 150, device = "pdf")
  cat(sprintf("  -> %s\n", pdf_file2))
  ggsave(png_file2, p_focus, width = 18, height = 16, dpi = 200, device = "png")
  cat(sprintf("  -> %s\n", png_file2))
} else { cat("  skipped (focused not in PLOT_VIEWS)\n") }


# ---- Plot 3: Publication-grade Co-function gene network ----
# Genes linked by shared GO / KEGG / PFAM / COG annotations
# Filters out overly generic root-level GO terms
cat("Building co-function gene network...\n")

# Generic root-level terms — filter only the truly universal GO level-1 terms
# that connect everything and provide no functional discrimination.
# NOTE: GO terms in the network have format GO:GO:XXXXXXX (Python f'GO:{go}')
# where the go value from Emapper already contains the GO: prefix.
ROOT_GO_TERMS <- c(
  "GO:GO:0003674",  # molecular_function
  "GO:GO:0008150",  # biological_process  
  "GO:GO:0005575",  # cellular_component
  "GO:GO:0003824",  # catalytic activity
  "GO:GO:0005488",  # binding
  "GO:GO:0005515",  # protein binding
  "GO:GO:0005622",  # intracellular
  "GO:GO:0005623",  # cell
  "GO:GO:0009987"   # cellular process
)
# COG:S (Function unknown) also adds noise — filter it from edge building
FILTER_COG_S <- TRUE

# Standard COG A-Z category full names
COG_DEFAULTS <- c(
  A="RNA processing and modification", B="Chromatin structure and dynamics",
  C="Energy production and conversion", D="Cell cycle control, cell division, chromosome partitioning",
  E="Amino acid transport and metabolism", F="Nucleotide transport and metabolism",
  G="Carbohydrate transport and metabolism", H="Coenzyme transport and metabolism",
  I="Lipid transport and metabolism", J="Translation, ribosomal structure and biogenesis",
  K="Transcription", L="Replication, recombination and repair",
  M="Cell wall/membrane/envelope biogenesis", N="Cell motility",
  O="Posttranslational modification, protein turnover, chaperones",
  P="Inorganic ion transport and metabolism", Q="Secondary metabolites biosynthesis, transport and catabolism",
  R="General function prediction only", S="Function unknown",
  T="Signal transduction mechanisms", U="Intracellular trafficking, secretion, and vesicular transport",
  V="Defense mechanisms", W="Extracellular structures", X="Mobilome: prophages, transposons",
  Y="Nuclear structure", Z="Cytoskeleton")

gene_names_all <- V(g)$name[V(g)$type %in% c("gene","transcript")]
gene_sig_all  <- V(g)$sig_num[V(g)$type %in% c("gene","transcript")]
names(gene_sig_all) <- gene_names_all

# Build descriptive labels for genes
node_desc_map <- setNames(
  ifelse(nodes$description != "" & nodes$description != "-" & !is.na(nodes$description),
         nodes$description, nodes$display_name),
  nodes$id
)
# Truncate long descriptions
node_desc_map <- sapply(node_desc_map, function(d) {
  if (nchar(d) > 55) paste0(substr(d, 1, 52), "...") else d
})

# Build gene↔function adjacency
func_edges_dt <- as.data.table(edges[type == "annotation"])
gene_co_pairs <- list()

for (i in seq_len(nrow(func_edges_dt))) {
  s <- func_edges_dt$source[i]
  t <- func_edges_dt$target[i]
  if (s %in% gene_names_all && !t %in% gene_names_all) {
    gene_co_pairs[[length(gene_co_pairs) + 1]] <- data.table(gene = s, func = t)
  } else if (t %in% gene_names_all && !s %in% gene_names_all) {
    gene_co_pairs[[length(gene_co_pairs) + 1]] <- data.table(gene = t, func = s)
  }
}

if (length(gene_co_pairs) > 0) {
  gf <- rbindlist(gene_co_pairs)

  # Filter out root GO terms that connect too broadly (noise)
  gf <- gf[!func %in% ROOT_GO_TERMS]
  # Also filter COG:S (Function unknown) — no biological signal for co-function
  if (FILTER_COG_S) {
    gf <- gf[!grepl("^COG:\\[S\\] Function unknown", func)]
  }
  cat(sprintf("  Gene-function pairs after filtering: %d\n", nrow(gf)))

  gf_list <- split(gf$gene, gf$func)

  # Build gene-gene pairs + track which functions are shared
  co_edges_list <- list()
  for (fn in names(gf_list)) {
    genes_in_fn <- unique(gf_list[[fn]])
    if (length(genes_in_fn) >= 2) {
      comb <- combn(sort(genes_in_fn), 2, simplify = FALSE)
      for (p in comb) {
        co_edges_list[[length(co_edges_list) + 1]] <- data.table(
          gene1 = p[1], gene2 = p[2], shared_func = fn
        )
      }
    }
  }

  if (length(co_edges_list) > 0) {
    co_long <- rbindlist(co_edges_list)

    # Aggregate: count shared functions per gene pair + collect function names
    co_agg <- co_long[, .(
      shared_count = .N,
      shared_funcs = paste(sort(unique(shared_func)), collapse = ";;")
    ), by = .(gene1, gene2)]

    # --- Quality filter: keep edges with biological signal ---
    # COG/PFAM/KEGG term → confident at 1 shared (specific curated annotation)
    # GO/Description term → needs ≥2 shared (needs corroboration)
    has_strong_term <- function(sfs) {
      terms <- unlist(strsplit(sfs, ";;"))
      any(grepl("^COG:(?!\\[S\\])", terms, perl = TRUE) |  # COG non-S
          grepl("^PFAM:", terms) |
          grepl("^KEGG:", terms))
    }
    co_agg[, strong := sapply(shared_funcs, has_strong_term)]
    # Keep if: (strong term and shared_count>=1) OR (non-strong but shared_count>=2)
    co_agg_sub <- co_agg[(strong == TRUE & shared_count >= 1) |
                         (strong == FALSE & shared_count >= 2)]
    threshold_used <- 2

    # Fallback
    if (nrow(co_agg_sub) < 5) {
      co_agg_sub <- co_agg[shared_count >= 1]
      threshold_used <- 1
    }
    cat(sprintf("  Co-function pairs: %d total, %d after strong-term filter, %d plotted\n",
                nrow(co_agg), nrow(co_agg[(strong == TRUE & shared_count >= 1) |
                                          (strong == FALSE & shared_count >= 2)]),
                nrow(co_agg_sub)))

    # Build a per-gene lookup table (plain data.table, no igraph dependency)
    genes_with_edges <- unique(c(co_agg_sub$gene1, co_agg_sub$gene2))
    # Also include orphan genes (no co-function partners) from top-N
    all_top_genes <- gene_names_all  # from GWAS top-N
    orphan_genes  <- setdiff(all_top_genes, genes_with_edges)
    cat(sprintf("  Genes with co-function edges: %d  |  Orphan genes (no edges): %d\n",
                length(genes_with_edges), length(orphan_genes)))

    # All genes: edged + orphans (orphans shown for completeness)
    genes_in_co <- unique(c(genes_with_edges, orphan_genes))
    gene_lookup <- data.table(
      name = genes_in_co,
      sig_num = as.numeric(gene_sig_all[genes_in_co]),
      description = as.character(node_desc_map[genes_in_co])
    )
    gene_lookup[is.na(sig_num), sig_num := min(gene_sig_all, na.rm = TRUE)]

    # Build igraph — orphans become isolated vertices (degree=0)
    g_co <- graph_from_data_frame(co_agg_sub, vertices = genes_in_co, directed = FALSE)

    # Attach attributes from lookup
    V(g_co)$sig_num      <- gene_lookup$sig_num[match(V(g_co)$name, gene_lookup$name)]
    V(g_co)$description  <- gene_lookup$description[match(V(g_co)$name, gene_lookup$name)]

    # Edge shared_count from co_agg_sub
    edge_key <- paste(co_agg_sub$gene1, co_agg_sub$gene2)
    edge_key_g <- paste(ends(g_co, E(g_co))[,1], ends(g_co, E(g_co))[,2])
    edge_key_g2 <- paste(ends(g_co, E(g_co))[,2], ends(g_co, E(g_co))[,1])
    match_idx <- match(edge_key_g, edge_key)
    match_idx[is.na(match_idx)] <- match(edge_key_g2[is.na(match_idx)], edge_key)
    E(g_co)$shared_count <- co_agg_sub$shared_count[match_idx]
    E(g_co)$shared_funcs  <- co_agg_sub$shared_funcs[match_idx]

    # Community detection — only on connected vertices
    set.seed(42)
    # Track orphans (degree 0) before clustering
    is_orphan <- degree(g_co) == 0
    n_orphans  <- sum(is_orphan)

    # Build a clean community vector (no NAs): orphans=0, connected=1..n_comm
    community_vec <- integer(vcount(g_co))
    if (n_orphans > 0 && any(!is_orphan)) {
      # Cluster only the connected component
      g_connected <- induced_subgraph(g_co, V(g_co)[!is_orphan])
      comm_conn <- cluster_louvain(g_connected)
      community_vec[!is_orphan] <- comm_conn$membership
      community_vec[is_orphan]  <- 0L
      n_comm <- length(unique(comm_conn$membership))
    } else if (n_orphans == vcount(g_co)) {
      n_comm <- 0
    } else {
      comm <- cluster_louvain(g_co)
      community_vec <- comm$membership
      n_comm <- length(unique(comm$membership))
    }
    # Store in igraph (no NAs)
    V(g_co)$community <- community_vec
    n_connected <- vcount(g_co) - n_orphans
    cat(sprintf("  Co-function network: %d genes (%d connected + %d orphan), %d edges, %d Louvain modules\n",
                vcount(g_co), n_connected, n_orphans, ecount(g_co), max(n_comm, 0)))

    # --- Characterize each module: pick the most specific shared function ---
    # Also determine the most significant gene per module (for dedup / labeling)
    all_shared_funcs <- strsplit(E(g_co)$shared_funcs, ";;")
    comm_by_name <- setNames(V(g_co)$community, V(g_co)$name)

    # Per-module top gene (most significant) — uses gene_lookup, robust
    module_top_genes <- sapply(seq_len(n_comm), function(m) {
      m_mask <- V(g_co)$community == m
      m_names <- V(g_co)$name[m_mask]
      if (length(m_names) == 0) return("?")
      # Lookup sig_num from gene_lookup (more reliable than V(g_co)$sig_num)
      m_rows <- gene_lookup[name %in% m_names]
      if (nrow(m_rows) == 0) return(m_names[1])
      tg <- m_rows[order(-sig_num)][1, name]
      if (length(tg) == 0 || is.na(tg)) tg <- m_names[1]
      tg
    })

    module_chars <- lapply(seq_len(n_comm), function(m) {
      m_genes <- V(g_co)$name[V(g_co)$community == m]
      if (length(m_genes) == 0) return(paste0("Module_", m))

      # Build COG letter→full name mapping from the nodes data
      # nodes data.table has columns: id, type, COG_category, COG_full_name, ...
      cog_full_names <- nodes[type %in% c("gene","transcript") & COG_full_name != "" & COG_full_name != "-",
                              unique(COG_full_name)]
      # Format from gwas2network.py: "[K] Transcription" or "[ATY] [A] ...; [T] ...; [Y] ..."
      # For single-letter: extract letter and name
      cog_map <- list()
      for (cfn in cog_full_names) {
        if (grepl("^\\[([A-Z])\\] (.+)$", cfn)) {
          letter <- gsub("^\\[([A-Z])\\].*", "\\1", cfn)
          if (!letter %in% names(cog_map)) {
            cog_map[[letter]] <- gsub("^\\[[A-Z]\\] ", "", cfn)
          }
        }
      }
      for (l in names(COG_DEFAULTS)) {
        if (!l %in% names(cog_map)) cog_map[[l]] <- COG_DEFAULTS[l]
      }

      # --- Frequency-based naming ---
      gene_cogs <- gene_lookup[name %in% m_genes]
      # nodes data.table uses 'id' not 'name' for gene identifier
      m_cog_raw <- nodes[type %in% c("gene","transcript") & id %in% m_genes, COG_category]
      cog_letters <- unlist(strsplit(paste(m_cog_raw, collapse = ""), ""))
      cog_letters <- cog_letters[cog_letters != "-" & cog_letters != ""]

      if (length(cog_letters) > 0) {
        cog_freq <- sort(table(cog_letters), decreasing = TRUE)
        top_letters <- names(cog_freq)
        best_letter <- head(top_letters[top_letters != "S"], 1)
        if (length(best_letter) == 0) best_letter <- top_letters[1]
        cog_name <- cog_map[[best_letter]]
        if (is.null(cog_name) || is.na(cog_name)) cog_name <- paste0("COG:", best_letter)
        best <- paste0("[", best_letter, "] ", cog_name)
      } else {
        # Fallback: use edge-based function terms
        edge_ends <- ends(g_co, E(g_co))
        m_edge_mask <- (comm_by_name[edge_ends[,1]] == m |
                        comm_by_name[edge_ends[,2]] == m)
        m_funcs <- unique(unlist(all_shared_funcs[m_edge_mask]))
        pfam_terms <- grep("^PFAM:", m_funcs, value = TRUE)
        kegg_terms <- grep("^KEGG:", m_funcs, value = TRUE)
        desc_terms <- grep("^DESC:", m_funcs, value = TRUE)
        candidates <- c(
          gsub("^PFAM:", "", pfam_terms),
          gsub("^KEGG:", "", kegg_terms),
          gsub("^DESC:", "", desc_terms)
        )
        best <- head(candidates[candidates != "Function unknown"], 1)
        if (length(best) == 0) best <- head(candidates, 1)
        if (length(best) == 0) best <- paste0("Module_", m)
      }
      best[1]
    })
    module_names <- unlist(module_chars)
    # Truncate
    module_names <- sapply(module_names, function(x) {
      if (nchar(x) > 60) paste0(substr(x, 1, 57), "...") else x
    })
    names(module_names) <- as.character(seq_len(n_comm))

    # --- Make names unique — append top gene of each module ---
    # Always show top gene for richer module annotation
    seen <- character()
    for (i in seq_along(module_names)) {
      nm <- module_names[i]
      tg <- module_top_genes[i]
      if (nm %in% seen) {
        module_names[i] <- paste0(nm, " (", tg, ")")
      } else if (grepl("^Module_\\d+$", nm)) {
        module_names[i] <- paste0("Co-function cluster (", tg, ")")
      }
      seen <- c(seen, nm)
    }

    # --- Filter out "Function unknown" modules -------------------------------
    # Modules dominated by COG:S ("Function unknown") provide no biological
    # signal worth visualizing — they just add visual noise.  Remove those
    # genes from the co-function graph and report what was dropped.
    unknown_mask <- grepl("^\\[S\\]", module_names)
    unknown_genes_dropped <- character(0)
    if (any(unknown_mask)) {
      unknown_comms <- which(unknown_mask)
      unknown_genes_dropped <- V(g_co)$name[V(g_co)$community %in% unknown_comms]
      cat(sprintf(
        "  Dropping %d 'Function unknown' module%s (%d genes): %s\n    Genes: %s\n",
        length(unknown_comms), if (length(unknown_comms) > 1) "s" else "",
        length(unknown_genes_dropped),
        paste(sapply(unknown_comms, function(m) sprintf("M%d", m)), collapse = ", "),
        paste(head(unknown_genes_dropped, 10), collapse = ", ")))
      # Remove these vertices from the co-function graph entirely
      g_co <- delete_vertices(g_co, unknown_genes_dropped)
      # Update orphans / communities for the remaining graph
      is_orphan           <- degree(g_co) == 0
      n_orphans           <- sum(is_orphan)
      n_connected         <- vcount(g_co) - n_orphans
      # Renumber non-orphan communities to stay contiguous 1..n_comm
      if (n_connected > 1) {
        g_conn_only   <- induced_subgraph(g_co, V(g_co)[!is_orphan])
        comm_new      <- cluster_louvain(g_conn_only)
        n_comm        <- length(unique(comm_new$membership))
        community_vec <- integer(vcount(g_co))
        community_vec[!is_orphan] <- comm_new$membership
        community_vec[is_orphan]  <- 0L
      } else {
        n_comm        <- if (n_connected > 0) 1L else 0L
        community_vec <- ifelse(is_orphan, 0L, 1L)
      }
      V(g_co)$community <- community_vec
      # Refresh shared-funcs and community lookup after vertex deletion
      all_shared_funcs <- strsplit(E(g_co)$shared_funcs, ";;")
      comm_by_name      <- setNames(V(g_co)$community, V(g_co)$name)
      # Rebuild COG map for re-characterization (needed inside lapply closure)
      cog_full_names2 <- nodes[type %in% c("gene","transcript") & COG_full_name != "" & COG_full_name != "-",
                               unique(COG_full_name)]
      cog_map2 <- list()
      for (cfn in cog_full_names2) {
        if (grepl("^\\[([A-Z])\\] (.+)$", cfn)) {
          letter <- gsub("^\\[([A-Z])\\].*", "\\1", cfn)
          if (!letter %in% names(cog_map2)) {
            cog_map2[[letter]] <- gsub("^\\[[A-Z]\\] ", "", cfn)
          }
        }
      }
      for (l in names(COG_DEFAULTS)) {
        if (!l %in% names(cog_map2)) cog_map2[[l]] <- COG_DEFAULTS[l]
      }
      # Re-characterize the remaining modules
      module_chars <- lapply(seq_len(max(n_comm, 0)), function(m) {
        m_genes <- V(g_co)$name[V(g_co)$community == m]
        if (length(m_genes) == 0) return(paste0("Module_", m))
        m_cog_raw <- nodes[type %in% c("gene","transcript") & id %in% m_genes, COG_category]
        cog_letters <- unlist(strsplit(paste(m_cog_raw, collapse = ""), ""))
        cog_letters <- cog_letters[cog_letters != "-" & cog_letters != ""]
        if (length(cog_letters) > 0) {
          cog_freq <- sort(table(cog_letters), decreasing = TRUE)
          top_letters <- names(cog_freq)
          best_letter <- head(top_letters[top_letters != "S"], 1)
          if (length(best_letter) == 0) best_letter <- top_letters[1]
          cog_name <- cog_map2[[best_letter]]
          if (is.null(cog_name) || is.na(cog_name)) cog_name <- paste0("COG:", best_letter)
          paste0("[", best_letter, "] ", cog_name)
        } else {
          edge_ends <- ends(g_co, E(g_co))
          m_edge_mask <- (comm_by_name[edge_ends[,1]] == m | comm_by_name[edge_ends[,2]] == m)
          m_funcs <- unique(unlist(all_shared_funcs[m_edge_mask]))
          pfam_terms <- grep("^PFAM:", m_funcs, value = TRUE)
          kegg_terms <- grep("^KEGG:", m_funcs, value = TRUE)
          desc_terms <- grep("^DESC:", m_funcs, value = TRUE)
          candidates <- c(gsub("^PFAM:", "", pfam_terms),
                          gsub("^KEGG:", "", kegg_terms),
                          gsub("^DESC:", "", desc_terms))
          best <- head(candidates[candidates != "Function unknown"], 1)
          if (length(best) == 0) best <- head(candidates, 1)
          if (length(best) == 0) best <- paste0("Module_", m)
          best[1]
        }
      })
      module_names <- unlist(module_chars)
      module_names <- sapply(module_names, function(x) {
        if (nchar(x) > 60) paste0(substr(x, 1, 57), "...") else x
      })
      names(module_names) <- as.character(seq_len(max(n_comm, 0)))
      # Recompute module top genes
      module_top_genes <- sapply(seq_len(max(n_comm, 0)), function(m) {
        m_names <- V(g_co)$name[V(g_co)$community == m]
        if (length(m_names) == 0) return("?")
        m_rows <- gene_lookup[name %in% m_names]
        if (nrow(m_rows) == 0) return(m_names[1])
        tg <- m_rows[order(-sig_num)][1, name]
        if (length(tg) == 0 || is.na(tg)) tg <- m_names[1]
        tg
      })
      # Make module names unique
      seen <- character()
      for (i in seq_along(module_names)) {
        nm <- module_names[i]
        tg <- module_top_genes[i]
        if (nm %in% seen) {
          module_names[i] <- paste0(nm, " (", tg, ")")
        } else if (grepl("^Module_\\d+$", nm)) {
          module_names[i] <- paste0("Co-function cluster (", tg, ")")
        }
        seen <- c(seen, nm)
      }
    }

    cat("  Module characterization:\n")
    if (n_comm > 0) {
      for (m in seq_len(n_comm)) {
        m_genes <- V(g_co)$name[V(g_co)$community == m]
        cat(sprintf("    Module %d (%d genes): %s\n", m, length(m_genes), module_names[m]))
      }
    }
    if (n_orphans > 0) {
      orphan_list <- V(g_co)$name[is_orphan]
      cat(sprintf("    Orphans (%d genes, no co-function partners): %s\n",
                  n_orphans, paste(head(orphan_list, 8), collapse = ", ")))
    }
    if (length(unknown_genes_dropped) > 0) {
      cat(sprintf("    Dropped (unknown function): %d genes\n", length(unknown_genes_dropped)))
    }

    # --- Layout: Kamada-Kawai for connected component, outer ring for orphans ---
    set.seed(42)
    if (n_connected > 2) {
      g_conn <- induced_subgraph(g_co, V(g_co)[!is_orphan])
      co_layout_conn <- layout_with_kk(g_conn)
    } else if (n_connected > 0) {
      g_conn <- induced_subgraph(g_co, V(g_co)[!is_orphan])
      co_layout_conn <- layout_with_fr(g_conn)
    } else {
      co_layout_conn <- matrix(nrow = 0, ncol = 2)
    }

    # Build full layout: connected genes at center, orphans in outer ring
    co_layout <- matrix(0, nrow = vcount(g_co), ncol = 2)
    if (n_connected > 0) {
      # Scale connected layout to radius ~3
      max_r_conn <- max(sqrt(co_layout_conn[,1]^2 + co_layout_conn[,2]^2))
      if (max_r_conn > 0) {
        co_layout_conn <- co_layout_conn / max_r_conn * 3.5
      }
      co_layout[!is_orphan, ] <- co_layout_conn
    }
    # Orphans: place in an outer ring, ordered by significance
    if (n_orphans > 0) {
      orphan_idx <- which(is_orphan)
      orphan_sigs <- V(g_co)$sig_num[is_orphan]
      orphan_order <- orphan_idx[order(-orphan_sigs)]
      orphan_angles <- seq(0, 2 * pi, length.out = n_orphans + 1)[1:n_orphans]
      orphan_r <- if (n_connected > 0) 5.5 else 4
      co_layout[orphan_order, 1] <- orphan_r * cos(orphan_angles)
      co_layout[orphan_order, 2] <- orphan_r * sin(orphan_angles)
    }

    co_nodes <- as.data.table(co_layout)
    setnames(co_nodes, c("x", "y"))
    co_nodes[, name := V(g_co)$name]
    co_nodes[, sig_num := V(g_co)$sig_num]
    community_raw <- V(g_co)$community
    co_nodes[, community := factor(community_raw)]
    co_nodes[, description := V(g_co)$description]
    co_nodes[, is_orphan := community_raw == 0]
    # Module label: orphans = generic; connected = named via module_names
    co_nodes[, module_label := ifelse(community_raw == 0,
                                       "No co-function partners",
                                       module_names[as.character(community_raw)])]

    # Normalize layout
    max_r_co <- max(sqrt(co_nodes$x^2 + co_nodes$y^2))
    if (max_r_co > 0) {
      co_nodes[, x := x / max_r_co * 5]
      co_nodes[, y := y / max_r_co * 5]
    }

    # Build edge data
    co_edges_plot <- as.data.table(as_data_frame(g_co, what = "edges"))
    co_edges_plot[, from_x := co_nodes$x[match(from, co_nodes$name)]]
    co_edges_plot[, from_y := co_nodes$y[match(from, co_nodes$name)]]
    co_edges_plot[, to_x   := co_nodes$x[match(to, co_nodes$name)]]
    co_edges_plot[, to_y   := co_nodes$y[match(to, co_nodes$name)]]
    co_edges_plot <- co_edges_plot[!is.na(from_x) & !is.na(to_x)]

    # --- Smart labeling: top degree per module + top significance overall ---
    co_deg <- degree(g_co)
    co_nodes[, degree := co_deg[name]]
    # Top 3 per module by degree
    top_per_module <- co_nodes[, .SD[order(-degree)][1:min(3, .N)], by = community]
    top_sig_global <- co_nodes[order(-sig_num)][1:min(15, nrow(co_nodes))]
    label_ids <- unique(c(top_per_module$name, top_sig_global$name))
    co_nodes[, show_label := name %in% label_ids]
    # Dual label: description (if available) + gene ID in smaller font below
    co_nodes[, has_desc := description != "" & !is.na(description) & description != name]
    co_nodes[, label_desc := ifelse(has_desc, description, "")]
    co_nodes[, label_id   := name]  # always show gene ID for labeled nodes

    # --- Extract dominant COG letter from each module name ---
    # module_names look like: "[K] Transcription" or "[O] Posttranslational... (TCONS_xxx)"
    module_cog_letter <- sapply(module_names, function(nm) {
      m <- regmatches(nm, regexpr("(?<=\\[)[A-Z](?=\\])", nm, perl = TRUE))
      if (length(m) == 0) "?" else m[1]
    })
    # For modules with multi-letter names like "[E] Amino acid...; [G] Carbohydrate..."
    # the first letter is already the dominant one from our naming logic

    # --- Color per COG letter (same letter = same color), orphans get grey ---
    if (n_comm > 0) {
      unique_letters <- unique(module_cog_letter)
    } else {
      unique_letters <- character(0)
    }
    n_unique <- length(unique_letters)
    academic_palette <- c(
      "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
      "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
      "#7E6148", "#B09C85", "#008B8B", "#CD853F"
    )
    if (n_unique <= length(academic_palette)) {
      cog_colors <- setNames(academic_palette[1:max(n_unique, 1)], unique_letters)
    } else {
      cog_colors <- setNames(viridis(n_unique, option = "turbo"), unique_letters)
    }
    # Map community → COG letter → color
    if (n_comm > 0) {
      comm_colors <- cog_colors[module_cog_letter]
      names(comm_colors) <- as.character(seq_len(n_comm))
    } else {
      comm_colors <- character(0)
    }

    # Legend: show unique COG categories + orphans
    if (n_comm > 0) {
      legend_labels <- sapply(unique_letters, function(letter) {
        mods_with_letter <- which(module_cog_letter == letter)
        short_name <- gsub(" \\[.+", "", module_names[mods_with_letter[1]])
        short_name <- gsub("^\\[[A-Z]\\] ", "", short_name)
        if (nchar(short_name) > 55) short_name <- paste0(substr(short_name, 1, 52), "...")
        sprintf("[%s] %s  (%d module%s)", letter, short_name,
                length(mods_with_letter), ifelse(length(mods_with_letter) > 1, "s", ""))
      })
    } else {
      legend_labels <- character(0)
    }
    if (n_orphans > 0) {
      # Orphan color is stored separately for the orphan plot
      orphan_color <- "#D0D0D0"
    }

    # Attach cog_letter to co_nodes
    co_nodes[, cog_letter := ifelse(is_orphan, "orphan",
                                     module_cog_letter[as.character(community)])]

    co_sig_range <- range(co_nodes$sig_num, na.rm = TRUE)

    # --- Build publication plot (connected genes only) ---
    # Community hull data — skip orphans
    hull_data <- co_nodes[!is_orphan, .SD[chull(x, y)], by = .(community, cog_letter)]
    # Connected-only node/edge sets
    co_nodes_conn <- co_nodes[!is_orphan]
    co_edges_conn <- co_edges_plot  # edges never involve orphans by construction

    p_cofunc <- ggplot() +
      # 1. Community hulls (light fill, colored by COG letter)
      geom_polygon(
        data = hull_data,
        aes(x = x, y = y, fill = cog_letter, group = community),
        alpha = 0.10, color = NA, show.legend = FALSE
      ) +
      # 2. Edges: darker/thicker for more shared functions
      geom_segment(
        data = co_edges_conn,
        aes(x = from_x, y = from_y, xend = to_x, yend = to_y,
            linewidth = shared_count, alpha = shared_count),
        color = "#404040"
      ) +
      # 3. Gene nodes: filled by COG letter, sized by GWAS significance
      geom_point(
        data = co_nodes_conn,
        aes(x = x, y = y, fill = cog_letter, size = sig_num),
        shape = 21, color = "white", stroke = 0.4, alpha = 0.92
      ) +
      # 4. High-degree hub nodes get a gold highlight ring
      geom_point(
        data = co_nodes_conn[degree >= 3],
        aes(x = x, y = y, size = sig_num),
        shape = 21, fill = NA, color = "#B8860B", stroke = 1.5, alpha = 0.7
      ) +
      # 5a. Labels: gene ID (bold, above node)
      geom_text_repel(
        data = co_nodes_conn[show_label == TRUE],
        aes(x = x, y = y, label = label_id),
        size = 2.8, fontface = "bold", max.overlaps = 50,
        box.padding = 0.5, point.padding = 0.2,
        segment.color = "grey60", segment.size = 0.25,
        nudge_y = 0.12, show.legend = FALSE,
        color = "grey15"
      ) +
      # 5b. Labels: functional description (italic, below ID, only if available)
      geom_text_repel(
        data = co_nodes_conn[show_label == TRUE & has_desc == TRUE],
        aes(x = x, y = y, label = label_desc),
        color = "grey35", size = 2.0, fontface = "italic",
        max.overlaps = 50, box.padding = 0.5, point.padding = 0.2,
        segment.color = NA, nudge_y = -0.18, show.legend = FALSE
      ) +
      # 6. Module name labels (at hull centroid)
      geom_label(
        data = co_nodes_conn[, .(
          x = mean(x), y = mean(y),
          module_label = module_label[1],
          cog_letter = cog_letter[1]
        ), by = community],
        aes(x = x, y = y, label = module_label, fill = cog_letter),
        color = "white", size = 2.6, fontface = "bold",
        alpha = 0.88, label.padding = unit(0.3, "lines"),
        linewidth = 0.2, show.legend = FALSE
      ) +
      # --- Scales ---
      scale_fill_manual(
        name = "COG category",
        values = cog_colors,
        labels = legend_labels,
        guide = guide_legend(order = 1, override.aes = list(size = 5, alpha = 0.9))
      ) +
      scale_size_continuous(
        name = expression(-log[10](italic(p))),
        range = c(2.5, 10),
        guide = guide_legend(order = 2)
      ) +
      scale_linewidth_continuous(
        name = "Shared\nfunctions",
        range = c(0.35, 3.0),
        guide = guide_legend(order = 3, override.aes = list(alpha = 0.6))
      ) +
      scale_alpha_continuous(
        name = "Shared\nfunctions",
        range = c(0.18, 0.65),
        guide = guide_legend(order = 3)
      ) +
      coord_fixed(clip = "off") +
      theme_void() +
      theme(
        plot.background  = element_rect(fill = "white", color = NA),
        plot.title       = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(size = 9.5, hjust = 0.5, color = "grey35",
                                        margin = margin(t = 4, b = 10)),
        plot.margin      = margin(25, 25, 25, 25),
        legend.position  = "right",
        legend.title     = element_text(size = 9, face = "bold"),
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.45, "cm"),
        legend.spacing.y = unit(0.15, "cm"),
        legend.box       = "vertical",
        legend.margin    = margin(0, 0, 0, 0)
      ) +
      labs(
        title = "Co-Function Network of GWAS-Significant Transcripts",
        subtitle = sprintf(
          paste0(pheno_display, " -- %d connected transcripts linked by shared GO, KEGG, PFAM and COG terms  |  %d edges  |  %d functional modules"),
          nrow(co_nodes_conn), nrow(co_edges_conn), n_comm
        )
      )

    # --- Exports ---
    if ("cofunction" %in% views_wanted) {
      pdf_file3 <- paste0(prefix, "_network_cofunction.pdf")
      png_file3 <- paste0(prefix, "_network_cofunction.png")
      tiff_file3 <- paste0(prefix, "_network_cofunction.tiff")
      ggsave(pdf_file3, p_cofunc, width = 16, height = 14, dpi = 150, device = "pdf")
      cat(sprintf("  -> %s\n", pdf_file3))
      ggsave(png_file3, p_cofunc, width = 16, height = 14, dpi = 300, device = "png")
      cat(sprintf("  -> %s\n", png_file3))
      ggsave(tiff_file3, p_cofunc, width = 16, height = 14, dpi = 300, device = "tiff",
             compression = "lzw")
      cat(sprintf("  -> %s\n", tiff_file3))
    } else { cat("  skipped (cofunction not in PLOT_VIEWS)\n") }

    # --- Export co-function edge table with shared function names ---
    co_edge_export <- as.data.table(as_data_frame(g_co, what = "edges"))
    co_edge_export[, shared_function_count := shared_count]
    co_edge_out <- paste0(prefix, "_cofunction_edges.tsv")
    fwrite(co_edge_export, co_edge_out, sep = "\t")
    cat(sprintf("  -> %s  (%d edges)\n", co_edge_out, nrow(co_edge_export)))

    # --- Module summary (connected genes only, skip orphans) ---
    module_summary <- co_nodes[!is_orphan, .(
      n_genes = .N,
      mean_sig = round(mean(sig_num, na.rm = TRUE), 2),
      max_sig  = round(max(sig_num, na.rm = TRUE), 2),
      characterizing_function = module_label[1],
      top_genes = paste(name[order(-sig_num)][1:min(5, .N)], collapse = ", ")
    ), by = community][order(-n_genes)]
    module_out <- paste0(prefix, "_cofunction_modules.tsv")
    fwrite(module_summary, module_out, sep = "\t")
    cat(sprintf("  -> %s  (%d modules)\n", module_out, nrow(module_summary)))

    # Print module summary to console
    cat("\n  Co-Function Module Summary:\n")
    for (i in seq_len(nrow(module_summary))) {
      row <- module_summary[i]
      cat(sprintf("    M%d: %s  |  %d genes  |  mean -log10(p)=%.2f  |  max=%.2f\n",
                  i, row$characterizing_function,
                  row$n_genes, row$mean_sig, row$max_sig))
    }

    # --- Orphan-only plot (genes with no co-function partners) ---
    if (n_orphans > 0) {
      cat("\n  Building orphan gene plot...\n")
      co_nodes_orph <- co_nodes[is_orphan == TRUE]

      # Label top 20 most significant orphans
      co_nodes_orph[, show_label := sig_num >= sort(sig_num, decreasing = TRUE)[min(20, .N)]]
      co_nodes_orph[, label_text := name]

      # Color orphans by their COG category
      # Build orphan COG mapping from nodes data
      orphan_cog_raw <- nodes[type %in% c("gene","transcript") & id %in% co_nodes_orph$name, .(id, COG_category)]
      orphan_cog_map <- setNames(orphan_cog_raw$COG_category, orphan_cog_raw$id)
      # Extract first letter for coloring (or S for unknown)
      co_nodes_orph[, cog_letter := sapply(name, function(g) {
        raw <- orphan_cog_map[g]
        if (is.na(raw) || raw == "" || raw == "-") return("?")
        substr(raw, 1, 1)
      })]

      # Build orphan color mapping using the same palette as connected genes
      orphan_unique_letters <- unique(co_nodes_orph$cog_letter)
      orphan_palette <- c(
        "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
        "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
        "#7E6148", "#B09C85", "#008B8B", "#CD853F", "#A0A0A0"
      )
      if (length(orphan_unique_letters) <= length(orphan_palette)) {
        orphan_colors <- setNames(orphan_palette[1:length(orphan_unique_letters)], orphan_unique_letters)
      } else {
        orphan_colors <- setNames(viridis(length(orphan_unique_letters), option = "turbo"), orphan_unique_letters)
      }

      # Legend: COG letter → full name
      orphan_legend <- sapply(names(orphan_colors), function(letter) {
        cn <- COG_DEFAULTS[letter]
        if (is.na(cn)) cn <- paste0("COG:", letter)
        n <- sum(co_nodes_orph$cog_letter == letter)
        sprintf("[%s] %s  (%d gene%s)", letter, cn, n, ifelse(n > 1, "s", ""))
      })

      p_orphan <- ggplot(co_nodes_orph) +
        geom_point(
          aes(x = x, y = y, fill = cog_letter, size = sig_num),
          shape = 21, color = "grey40", stroke = 0.3, alpha = 0.9
        ) +
        geom_text_repel(
          data = co_nodes_orph[show_label == TRUE],
          aes(x = x, y = y, label = label_text),
          size = 2.5, color = "grey20", fontface = "italic",
          max.overlaps = 30, box.padding = 0.4, point.padding = 0.2,
          segment.color = "grey60", segment.size = 0.2,
          show.legend = FALSE
        ) +
        scale_fill_manual(
          name = "COG category",
          values = orphan_colors,
          labels = orphan_legend,
          guide = guide_legend(order = 1, override.aes = list(size = 4, alpha = 0.9))
        ) +
        scale_size_continuous(
          name = expression(-log[10](italic(p))),
          range = c(2, 7),
          guide = guide_legend(order = 2)
        ) +
        coord_fixed(clip = "off") +
        theme_void() +
        theme(
          plot.background  = element_rect(fill = "white", color = NA),
          plot.title       = element_text(size = 15, face = "bold", hjust = 0.5),
          plot.subtitle    = element_text(size = 9.5, hjust = 0.5, color = "grey35",
                                          margin = margin(t = 4, b = 10)),
          plot.margin      = margin(25, 25, 25, 25),
          legend.position  = "right",
          legend.title     = element_text(size = 9, face = "bold"),
          legend.text      = element_text(size = 8)
        ) +
        labs(
          title = "Orphan GWAS-Significant Transcripts",
          subtitle = sprintf(
            paste0(pheno_display, " -- %d genes with no co-function partners  |  significant but unlinked in the annotation network"),
            n_orphans
          )
        )

      if ("orphans" %in% views_wanted) {
        pdf_file4 <- paste0(prefix, "_network_orphans.pdf")
        png_file4 <- paste0(prefix, "_network_orphans.png")
        ggsave(pdf_file4, p_orphan, width = 14, height = 12, dpi = 150, device = "pdf")
        cat(sprintf("  -> %s\n", pdf_file4))
        ggsave(png_file4, p_orphan, width = 14, height = 12, dpi = 300, device = "png")
        cat(sprintf("  -> %s\n", png_file4))
      }

      # Export orphan gene list
      orphan_export <- co_nodes_orph[order(-sig_num), .(name, sig_num, description)]
      orphan_out <- paste0(prefix, "_orphan_genes.tsv")
      fwrite(orphan_export, orphan_out, sep = "\t")
      cat(sprintf("  -> %s  (%d genes)\n", orphan_out, nrow(orphan_export)))

      # --- Combined plot: connected modules (center) + orphans (outer ring) ---
      cat("  Building combined cofunction + orphan plot...\n")

      # Scale connected positions — 75% area (radius ~7.5, orphan rim at 8.7)
      max_r_conn2 <- max(sqrt(co_nodes_conn$x^2 + co_nodes_conn$y^2))
      if (max_r_conn2 > 0) {
        co_nodes_conn_combined <- copy(co_nodes_conn)
        co_nodes_conn_combined[, x := x / max_r_conn2 * 7.5]
        co_nodes_conn_combined[, y := y / max_r_conn2 * 7.5]
        # Rebuild edge coordinates
        co_edges_combined <- copy(co_edges_conn)
        co_edges_combined[, from_x := co_nodes_conn_combined$x[match(from, co_nodes_conn_combined$name)]]
        co_edges_combined[, from_y := co_nodes_conn_combined$y[match(from, co_nodes_conn_combined$name)]]
        co_edges_combined[, to_x   := co_nodes_conn_combined$x[match(to, co_nodes_conn_combined$name)]]
        co_edges_combined[, to_y   := co_nodes_conn_combined$y[match(to, co_nodes_conn_combined$name)]]
        co_edges_combined <- co_edges_combined[!is.na(from_x) & !is.na(to_x)]
      }

      # Orphans in thin outer ring
      co_nodes_orph_combined <- copy(co_nodes_orph)
      n_o <- nrow(co_nodes_orph_combined)
      orphan_angles2 <- pi/2 - seq(0, 2 * pi, length.out = n_o + 1)[1:n_o]
      co_nodes_orph_combined[, x := 8.7 * cos(orphan_angles2)]
      co_nodes_orph_combined[, y := 8.7 * sin(orphan_angles2)]

      # Show ALL orphan gene IDs
      co_nodes_orph_combined[, show_label_orph := sig_num >= sort(sig_num, decreasing = TRUE)[min(25, .N)]]

      # Combined hull data (scaled)
      hull_data_combined <- co_nodes_conn_combined[, .SD[chull(x, y)], by = .(community, cog_letter)]

      # Thin separator ring
      ring_theta <- seq(0, 2*pi, length.out = 200)
      ring_path_combined <- data.table(x = 8.1 * cos(ring_theta), y = 8.1 * sin(ring_theta))

      # Combine orphan and connected colors for unified legend
      combined_legend_labels <- legend_labels
      if (length(orphan_legend) > 0) {
        orphan_legend_entry <- c("orphan" = sprintf("No partners  (%d gene%s)",
                                    n_orphans, ifelse(n_orphans > 1, "s", "")))
        combined_legend_labels <- c(combined_legend_labels, orphan_legend_entry)
      }

      p_combined <- ggplot() +
        # Separator ring
        geom_path(
          data = ring_path_combined,
          aes(x = x, y = y),
          color = "grey75", linewidth = 0.4, linetype = "dashed", alpha = 0.8
        ) +
        # Ring label
        annotate("text", x = 0, y = 8.85, label = "Orphan genes (no co-function partners)",
                 size = 2.6, color = "grey55", fontface = "italic") +
        annotate("text", x = 0, y = -1.2, label = "Co-function modules",
                 size = 3.0, color = "grey50", fontface = "italic") +
        # Connected edges
        geom_segment(
          data = co_edges_combined,
          aes(x = from_x, y = from_y, xend = to_x, yend = to_y,
              linewidth = shared_count, alpha = shared_count),
          color = "#404040"
        ) +
        # Connected hulls
        geom_polygon(
          data = hull_data_combined,
          aes(x = x, y = y, fill = cog_letter, group = community),
          alpha = 0.10, color = NA, show.legend = FALSE
        ) +
        # Connected gene nodes
        geom_point(
          data = co_nodes_conn_combined,
          aes(x = x, y = y, fill = cog_letter, size = sig_num),
          shape = 21, color = "white", stroke = 0.35, alpha = 0.92
        ) +
        # Orphan gene nodes — same size scale as connected genes (based on -log10(p))
        geom_point(
          data = co_nodes_orph_combined,
          aes(x = x, y = y, size = sig_num),
          fill = "#D0D0D0", color = "grey50", shape = 21, stroke = 0.3, alpha = 0.72
        ) +
        # Connected labels (top hubs only)
        geom_text_repel(
          data = co_nodes_conn_combined[degree >= 3 & show_label == TRUE],
          aes(x = x, y = y, label = label_id),
          size = 2.5, fontface = "bold", max.overlaps = 40,
          box.padding = 0.3, segment.color = "grey60", segment.size = 0.2,
          show.legend = FALSE, color = "grey15"
        ) +
        # Orphan labels — top significant ones only (others shown as dots)
        geom_text_repel(
          data = co_nodes_orph_combined[show_label_orph == TRUE],
          aes(x = x, y = y, label = name),
          size = 1.8, color = "grey35", fontface = "italic",
          max.overlaps = 40, box.padding = 0.3, point.padding = 0.1,
          segment.color = "grey70", segment.size = 0.15,
          force = 1, show.legend = FALSE
        ) +
        # Module name labels
        geom_label(
          data = co_nodes_conn_combined[, .(
            x = mean(x), y = mean(y),
            module_label = module_label[1],
            cog_letter = cog_letter[1]
          ), by = community],
          aes(x = x, y = y, label = module_label, fill = cog_letter),
          color = "white", size = 2.8, fontface = "bold",
          alpha = 0.88, label.padding = unit(0.30, "lines"),
          linewidth = 0.15, show.legend = FALSE
        ) +
        scale_fill_manual(
          name = "COG category",
          values = c(cog_colors, c("orphan" = "#D0D0D0")),
          labels = combined_legend_labels,
          guide = guide_legend(order = 1, override.aes = list(size = 4, alpha = 0.9))
        ) +
        scale_size_continuous(
          name = expression(-log[10](italic(p))),
          range = c(2.5, 8),
          guide = guide_legend(order = 2, override.aes = list(fill = "grey40", alpha = 0.8))
        ) +
        scale_linewidth_continuous(
          name = "Shared\nfunctions",
          range = c(0.2, 2.0),
          guide = guide_legend(order = 3)
        ) +
        scale_alpha_continuous(
          name = "Shared\nfunctions",
          range = c(0.15, 0.55),
          guide = guide_legend(order = 3)
        ) +
        coord_fixed(xlim = c(-10, 10), ylim = c(-10, 10), clip = "off") +
        theme_void() +
        theme(
          plot.background  = element_rect(fill = "white", color = NA),
          plot.title       = element_text(size = 15, face = "bold", hjust = 0.5),
          plot.subtitle    = element_text(size = 9.5, hjust = 0.5, color = "grey35",
                                          margin = margin(t = 4, b = 10)),
          plot.margin      = margin(25, 25, 25, 25),
          legend.position  = "right",
          legend.title     = element_text(size = 9, face = "bold"),
          legend.text      = element_text(size = 8)
        ) +
        labs(
          title = "Co-Function Network with Orphan Genes",
          subtitle = sprintf(
            paste0(pheno_display, " -- %d connected genes (%d modules) + %d orphans (no co-function partners)  |  %d edges"),
            nrow(co_nodes_conn_combined), n_comm, n_orphans, nrow(co_edges_combined)
          )
        )

      if ("cofunction_orphans" %in% views_wanted) {
        pdf_file5 <- paste0(prefix, "_network_cofunction_with_orphans.pdf")
        png_file5 <- paste0(prefix, "_network_cofunction_with_orphans.png")
        ggsave(pdf_file5, p_combined, width = 18, height = 16, dpi = 150, device = "pdf")
        cat(sprintf("  -> %s\n", pdf_file5))
        ggsave(png_file5, p_combined, width = 18, height = 16, dpi = 300, device = "png")
        cat(sprintf("  -> %s\n", png_file5))
      }
    }

  } else {
    cat("  WARNING: No co-function gene pairs found\n")
  }
} else {
  cat("  WARNING: No annotation edges — skipping co-function network\n")
}

cat("\nDone! Generated:\n")
cat(sprintf("  %s\n", pdf_file))
cat(sprintf("  %s\n", png_file))
if (exists("pdf_file2")) {
  cat(sprintf("  %s\n", pdf_file2))
  cat(sprintf("  %s\n", png_file2))
}
if (exists("pdf_file3")) {
  cat(sprintf("  %s\n", pdf_file3))
  cat(sprintf("  %s\n", png_file3))
  cat(sprintf("  %s\n", tiff_file3))
  cat(sprintf("  %s\n", co_edge_out))
  cat(sprintf("  %s\n", module_out))
}
if (exists("pdf_file4")) {
  cat(sprintf("  %s\n", pdf_file4))
  cat(sprintf("  %s\n", png_file4))
  cat(sprintf("  %s\n", orphan_out))
}
if (exists("pdf_file5")) {
  cat(sprintf("  %s\n", pdf_file5))
  cat(sprintf("  %s\n", png_file5))
}

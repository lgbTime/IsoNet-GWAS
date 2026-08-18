#!/usr/bin/env Rscript
#
# plot_cofunction_network.R — Publication-quality co-function network plot
#
# Reads nodes.tsv / edges.tsv / modules.tsv produced by gwas2network.py
# --build-network and renders a co-function network in the style of
# plot_network.R:
#
#   - Full network: all modules + singletons, concentric/modular layout
#   - Focused view: top modules + their connected isoforms (FR layout)
#   - Module map: each module as a colored hull with labeled hubs
#
# Usage:
#   Rscript plot_cofunction_network.R \
#       metabolites_cofunction_nodes.tsv \
#       metabolites_cofunction_edges.tsv \
#       metabolites_cofunction_modules.tsv \
#       [metabolites_network_ai_analysis.json | NONE] \
#       [output_prefix]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(igraph)
  library(RColorBrewer)
  library(viridisLite)
  library(ggrepel)
})

# ---- Command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_cofunction_network.R <nodes.tsv> <edges.tsv> <modules.tsv> [ai.json] [prefix]")
}
nodes_file   <- args[1]
edges_file   <- args[2]
modules_file <- args[3]
ai_file      <- if (length(args) >= 4 && args[4] != "NONE") args[4] else NULL
prefix       <- if (length(args) >= 5) args[5] else "cofunction_network"

# ---- Read data ----
cat(sprintf("Reading: %s\n", nodes_file))
nodes_raw <- fread(nodes_file)
cat(sprintf("Reading: %s\n", edges_file))
edges_raw <- fread(edges_file)
cat(sprintf("Reading: %s\n", modules_file))
modules <- fread(modules_file)

cat(sprintf("  Nodes: %d  |  Edges: %d  |  Modules: %d\n",
            nrow(nodes_raw), nrow(edges_raw), nrow(modules)))

# ---- Load AI analysis if available ----
ai_data <- NULL
if (!is.null(ai_file) && file.exists(ai_file)) {
  cat(sprintf("Reading AI analysis: %s\n", ai_file))
  ai_data <- jsonlite::fromJSON(ai_file, simplifyVector = FALSE)
}

# ---- Build igraph ----
all_node_ids <- unique(c(edges_raw$source, edges_raw$target))
# Include nodes that appear in nodes.tsv but not in edges (singletons)
extra_ids <- setdiff(nodes_raw$Isoform, all_node_ids)
all_node_ids <- c(all_node_ids, extra_ids)

g <- graph_from_data_frame(edges_raw, vertices = all_node_ids, directed = FALSE)

# Attach node attributes from nodes_raw
idx <- match(V(g)$name, nodes_raw$Isoform)
V(g)$Isoform      <- V(g)$name
V(g)$ref_Gene     <- nodes_raw$ref_Gene[idx]
V(g)$Pvalue       <- nodes_raw$Pvalue[idx]
V(g)$Annotation   <- nodes_raw$Annotation[idx]
V(g)$PFAM         <- nodes_raw$PFAM[idx]
V(g)$GO_desc      <- nodes_raw$GO_description[idx]
V(g)$KEGG_path    <- nodes_raw$KEGG_Pathway[idx]
V(g)$COG_category <- nodes_raw$COG_category[idx]
V(g)$Metabolites  <- nodes_raw$Metabolites[idx]
V(g)$Chr          <- nodes_raw$Chr[idx]

# -log10(p) for sizing / coloring — careful with p=0
V(g)$sig_num <- -log10(pmax(V(g)$Pvalue, 1e-300, na.rm = TRUE))
V(g)$sig_num[is.na(V(g)$sig_num)] <- 0

# ---- Module membership ----
iso2mod <- list()
mod_meta <- list()  # module_id → list(size, metabolites)
for (i in seq_len(nrow(modules))) {
  iso_list <- unlist(strsplit(modules$isoforms[i], ",\\s*"))
  for (iso in iso_list) {
    iso2mod[[iso]] <- modules$module_id[i]
  }
  mod_meta[[modules$module_id[i]]] <- list(
    size = modules$size[i],
    metabolites = modules$metabolites[i]
  )
}
V(g)$module <- sapply(V(g)$name, function(x) {
  m <- iso2mod[[x]]
  if (is.null(m)) "singleton" else m
})

# Classify: big modules (≥3), small modules (2), singletons (1)
mod_sizes <- table(V(g)$module)
big_mods   <- names(mod_sizes[mod_sizes >= 3])
small_mods <- names(mod_sizes[mod_sizes == 2])
singletons <- names(mod_sizes[mod_sizes == 1])

n_big    <- length(big_mods)
n_small  <- length(small_mods)
n_single <- length(singletons)
cat(sprintf("  Modules: %d big (≥3)  |  %d small (=2)  |  %d singletons\n",
            n_big, n_small, n_single))

# ---- Node categories for coloring ----
n_big_mod_nodes <- sum(V(g)$module %in% big_mods)
cat(sprintf("  Nodes in big modules: %d  |  singletons: %d\n",
            n_big_mod_nodes, n_single))

# ---- Academic color palette for modules ----
academic_palette <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85", "#008B8B", "#CD853F",
  "#8B5A2B", "#556B2F", "#8B008B", "#2F4F4F",
  "#FF6F00", "#1B9E77", "#D95F02", "#7570B3"
)
if (n_big > 0) {
  big_colors <- setNames(
    academic_palette[1:min(n_big, length(academic_palette))],
    big_mods
  )
} else {
  big_colors <- character(0)
}
singleton_color <- "#D0D0D0"
small_color     <- "#B0B0B0"

# ---- Layout: FR with moderate inter-component spacing ----
# Standard FR works fine but spreads disconnected components arbitrarily.
# We use layout_with_fr with custom weights to keep things reasonable.
set.seed(42)
lay <- layout_with_fr(g, niter = 3000)
# Normalize to [-10, 10]
max_r <- max(sqrt(lay[,1]^2 + lay[,2]^2))
if (max_r > 0) {
  lay <- lay / max_r * 10
}

# ---- Build plot data.table ----
plot_nodes <- data.table(
  x       = lay[, 1],
  y       = lay[, 2],
  name    = V(g)$name,
  module  = V(g)$module,
  sig_num = V(g)$sig_num,
  annot   = V(g)$Annotation,
  pfam    = V(g)$PFAM,
  go_desc = V(g)$GO_desc,
  kegg    = V(g)$KEGG_path,
  cog     = V(g)$COG_category,
  degree  = degree(g),
  mets    = V(g)$Metabolites,
  chr     = V(g)$Chr
)
plot_nodes[, node_cat := fifelse(module %in% big_mods, "module_member",
                          fifelse(module %in% small_mods, "small_module",
                          "singleton"))]

# ---- Edge coordinates ----
plot_edges <- as.data.table(edges_raw)
plot_edges[, from_x := plot_nodes$x[match(source, plot_nodes$name)]]
plot_edges[, from_y := plot_nodes$y[match(source, plot_nodes$name)]]
plot_edges[, to_x   := plot_nodes$x[match(target, plot_nodes$name)]]
plot_edges[, to_y   := plot_nodes$y[match(target, plot_nodes$name)]]
plot_edges <- plot_edges[!is.na(from_x) & !is.na(to_x)]

# Edge weight for thickness/alpha
max_wt <- max(plot_edges$weight, 3)
plot_edges[, edge_lw := scales::rescale(weight, to = c(0.15, 2.0), from = c(1, max_wt))]
plot_edges[, edge_alpha := scales::rescale(weight, to = c(0.10, 0.55), from = c(1, max_wt))]

# ---- Module hulls ----
hull_data <- plot_nodes[module %in% big_mods, .SD[chull(x, y)], by = module]

# ---- Labels: top hubs per big module + top significance overall ----
top_hubs <- plot_nodes[module %in% big_mods & degree >= 2,
                       .SD[order(-degree)][1:min(3, .N)], by = module]
top_sig  <- plot_nodes[order(-sig_num)][1:min(20, .N)]
label_ids <- unique(c(top_hubs$name, top_sig$name))
plot_nodes[, show_label := name %in% label_ids]

# Clean up label text — use gene name if available, else isoform
plot_nodes[, label_text := name]

# ---- AI module names for centroids ----
module_names <- data.table(module = character(), ai_name = character(),
                            x = numeric(), y = numeric())
if (!is.null(ai_data) && !is.null(ai_data$modules)) {
  for (am in ai_data$modules) {
    mid <- am$module_id
    nm  <- am$name
    if (!is.null(mid) && !is.null(nm)) {
      mn <- plot_nodes[module == mid]
      if (nrow(mn) > 0) {
        module_names <- rbind(module_names, data.table(
          module = mid, ai_name = nm,
          x = mean(mn$x), y = mean(mn$y)
        ))
      }
    }
  }
  cat(sprintf("  AI names loaded for %d modules\n", nrow(module_names)))
}

# ---- Build module display labels — use best annotation of the hub gene ----
# Instead of raw term counts, pick the best-annotated hub gene (highest degree
# or most connected to annotation terms) and use its eggNOG Description + PFAM
# as the module label. This gives human-readable labels like:
#   "M2 PPR: Pentatricopeptide repeat-containing protein | PPR,PPR_2,PPR_3"
build_module_label <- function(m) {
  mm  <- mod_meta[[m]]
  mod_isos <- unlist(strsplit(modules[module_id == m]$isoforms, ",\\s*"))
  if (length(mod_isos) == 0 || is.na(mod_isos[1])) {
    return(sprintf("%s  (%d isoforms)", m, mm$size))
  }
  # Get annotations for module members from nodes
  mod_nodes <- plot_nodes[name %in% mod_isos]
  if (nrow(mod_nodes) == 0) return(sprintf("%s  (%d isoforms)", m, mm$size))

  # Best hub = highest degree within module
  mod_nodes <- mod_nodes[order(-degree)]
  best <- mod_nodes[1]

  # Build label: AI name (if available) + best gene annotation
  ai_name <- module_names[module == m]$ai_name
  if (length(ai_name) == 0 || is.na(ai_name) || ai_name == "") ai_name <- NULL

  parts <- c()
  # Part 1: AI functional name
  if (!is.null(ai_name)) parts <- c(parts, ai_name)

  # Part 2: best PFAM domains
  pfam <- best$pfam
  if (!is.na(pfam) && nchar(pfam) > 0 && pfam != "-") {
    pfam_clean <- gsub(",", ", ", pfam)
    if (nchar(pfam_clean) > 50) pfam_clean <- paste0(substr(pfam_clean, 1, 47), "...")
    parts <- c(parts, pfam_clean)
  }

  # Part 3: COG category from GFF annotation
  cog <- best$cog
  if (!is.na(cog) && nchar(cog) > 0 && cog != "-" && cog != "S") {
    # COG letter codes are expanded in the nodes file already
    parts <- c(parts, paste0("COG:", cog))
  }

  # Part 4: best eggNOG description (truncated)
  annot <- best$annot
  if (!is.na(annot) && nchar(annot) > 3 && annot != "-") {
    annot_clean <- gsub("\\s+", " ", annot)
    if (nchar(annot_clean) > 60) annot_clean <- paste0(substr(annot_clean, 1, 57), "...")
    parts <- c(parts, annot_clean)
  }

  # Part 5: top GO term
  go <- best$go_desc
  if (!is.na(go) && nchar(go) > 0) {
    go_first <- strsplit(go, ";\\s*")[[1]][1]
    if (!is.na(go_first) && nchar(go_first) < 50)
      parts <- c(parts, go_first)
  }

  if (length(parts) == 0) return(sprintf("%s  (%d isoforms)", m, mm$size))
  label <- paste(parts, collapse = "  |  ")
  if (nchar(label) > 120) label <- paste0(substr(label, 1, 117), "...")
  paste0(m, " ", label)
}

# Build legend labels and module display labels
if (n_big > 0) {
  mod_display_labels <- sapply(big_mods, build_module_label)
  names(mod_display_labels) <- big_mods
  legend_labels <- mod_display_labels
  names(legend_labels) <- big_mods
} else {
  legend_labels <- character(0)
  mod_display_labels <- character(0)
}

# ---- Significance range ----
sig_range <- range(plot_nodes$sig_num, na.rm = TRUE)
if (diff(sig_range) <= 0) sig_range <- c(0, 1)

# ===================================================================
# PLOT 1 — Full network
# ===================================================================
cat("Plotting full co-function network...\n")

# Pre-compute module centroids for label placement
mod_centroids <- plot_nodes[module %in% big_mods,
  .(x = mean(x), y = mean(y), size = .N, mod = module[1]), by = module]
# Attach enriched display labels
mod_centroids[, display_label := ifelse(
  module %in% names(mod_display_labels),
  mod_display_labels[module], module
)]

p_full <- ggplot() +
  # 1. Edges — continuous weight by shared annotation count
  geom_segment(
    data = plot_edges,
    aes(x = from_x, y = from_y, xend = to_x, yend = to_y,
        linewidth = weight, alpha = weight),
    color = "#505050"
  ) +
  # 2. Singleton nodes (grey, small, background)
  geom_point(
    data = plot_nodes[node_cat == "singleton"],
    aes(x = x, y = y, size = sig_num),
    fill = singleton_color, color = "grey70", shape = 21, stroke = 0.15, alpha = 0.45
  ) +
  # 3. Small module nodes
  geom_point(
    data = plot_nodes[node_cat == "small_module"],
    aes(x = x, y = y, size = sig_num),
    fill = small_color, color = "grey60", shape = 21, stroke = 0.2, alpha = 0.55
  ) +
  # 4. Big module hulls — light fill
  geom_polygon(
    data = hull_data,
    aes(x = x, y = y, fill = module, group = module),
    alpha = 0.09, color = NA, show.legend = FALSE
  ) +
  # 5. Big module nodes — colored fill, white border
  geom_point(
    data = plot_nodes[node_cat == "module_member"],
    aes(x = x, y = y, fill = module, size = sig_num),
    shape = 21, color = "white", stroke = 0.35, alpha = 0.90
  ) +
  # 6. Hub highlight ring (degree ≥ 3 in big modules)
  geom_point(
    data = plot_nodes[node_cat == "module_member" & degree >= 3],
    aes(x = x, y = y, size = sig_num),
    shape = 21, fill = NA, color = "#B8860B", stroke = 1.3, alpha = 0.60
  ) +
  # 7. Module name labels at centroid (with top characterizing terms)
  geom_label(
    data = mod_centroids,
    aes(x = x, y = y, label = display_label, fill = mod),
    color = "white", size = 2.6, fontface = "bold",
    alpha = 0.85, label.padding = unit(0.28, "lines"),
    linewidth = 0.2, show.legend = FALSE
  ) +
  # 8. AI module names (italic, below centroid labels)
  geom_text(
    data = module_names,
    aes(x = x, y = y, label = ai_name),
    size = 2.3, color = "grey30", fontface = "italic",
    nudge_y = -0.45, show.legend = FALSE
  ) +
  # 9. Isoform labels (ggrepel) — top hubs only
  geom_text_repel(
    data = plot_nodes[show_label == TRUE],
    aes(x = x, y = y, label = label_text),
    size = 2.0, color = "grey15", fontface = "bold",
    max.overlaps = 80, box.padding = 0.3, point.padding = 0.15,
    segment.color = "grey60", segment.size = 0.2,
    show.legend = FALSE
  ) +
  # ---- Scales ----
  scale_fill_manual(
    name = "Functional module",
    values = big_colors,
    labels = legend_labels,
    guide = guide_legend(order = 1, override.aes = list(size = 4, alpha = 0.85))
  ) +
  scale_size_continuous(
    name = expression(-log[10](italic(p))),
    range = c(1.5, 8),
    guide = guide_legend(order = 2, override.aes = list(fill = "grey40"))
  ) +
  scale_linewidth_continuous(
    name = "Shared\nannotations",
    range = c(0.15, 2.5),
    guide = guide_legend(order = 3, override.aes = list(alpha = 0.5))
  ) +
  scale_alpha_continuous(
    name = "Shared\nannotations",
    range = c(0.10, 0.55),
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
    legend.key.size  = unit(0.4, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.box       = "vertical",
    legend.margin    = margin(0, 0, 0, 0)
  ) +
  labs(
    title = "Co-Function Network of GWAS-Significant Isoforms",
    subtitle = sprintf(
      paste0("%d isoforms  |  %d edges (shared GO / PFAM)  |  %d functional modules  |  ",
             "node size = -log10(p)  |  edge thickness = shared terms"),
      nrow(plot_nodes), nrow(plot_edges), n_big + n_small
    )
  )

# ---- Export full ----
pdf_full <- paste0(prefix, "_network_cofunction.pdf")
png_full <- paste0(prefix, "_network_cofunction.png")
ggsave(pdf_full, p_full, width = 18, height = 16, dpi = 150, device = "pdf")
cat(sprintf("  → %s\n", pdf_full))
ggsave(png_full, p_full, width = 18, height = 16, dpi = 300, device = "png")
cat(sprintf("  → %s\n", png_full))

# ===================================================================
# PLOT 2 — Focused: top modules only, FR layout, more detail
# ===================================================================
if (n_big > 0) {
  cat("Plotting focused view (top modules)...\n")

  # Subgraph: keep nodes in big modules + their edges
  focused_nodes <- plot_nodes[node_cat == "module_member"]
  if (nrow(focused_nodes) > 4) {
    # Build subgraph from these nodes
    f_ids <- focused_nodes$name
    f_edges <- plot_edges[source %in% f_ids & target %in% f_ids]
    if (nrow(f_edges) > 0) {
      f_g <- graph_from_data_frame(f_edges, vertices = f_ids, directed = FALSE)

      # Attributes
      fidx <- match(V(f_g)$name, plot_nodes$name)
      V(f_g)$sig_num   <- plot_nodes$sig_num[fidx]
      V(f_g)$module    <- plot_nodes$module[fidx]
      V(f_g)$annot     <- plot_nodes$annot[fidx]
      V(f_g)$degree    <- degree(f_g)

      set.seed(42)
      f_lay <- layout_with_kk(f_g)
      max_fr <- max(sqrt(f_lay[,1]^2 + f_lay[,2]^2))
      if (max_fr > 0) f_lay <- f_lay / max_fr * 6

      f_nodes <- data.table(
        x = f_lay[, 1], y = f_lay[, 2],
        name = V(f_g)$name,
        module = V(f_g)$module,
        sig_num = V(f_g)$sig_num,
        degree = V(f_g)$degree
      )
      f_edges_dt <- as.data.table(f_edges)
      f_edges_dt[, from_x := f_nodes$x[match(source, f_nodes$name)]]
      f_edges_dt[, from_y := f_nodes$y[match(source, f_nodes$name)]]
      f_edges_dt[, to_x   := f_nodes$x[match(target, f_nodes$name)]]
      f_edges_dt[, to_y   := f_nodes$y[match(target, f_nodes$name)]]
      f_edges_dt <- f_edges_dt[!is.na(from_x) & !is.na(to_x)]
      if (nrow(f_edges_dt) > 0) {
        f_edges_dt[, edge_lw := scales::rescale(weight, to = c(0.3, 3.0),
                                                from = c(1, max(weight, 3)))]
        f_edges_dt[, edge_alpha := scales::rescale(weight, to = c(0.15, 0.65),
                                                   from = c(1, max(weight, 3)))]
      }

      f_hulls <- f_nodes[, .SD[chull(x, y)], by = module]
      f_hub_ids <- f_nodes[degree >= 2, name]
      f_nodes[, show_label := name %in% f_hub_ids]
      f_centroids <- f_nodes[, .(x = mean(x), y = mean(y), mod = module[1]),
                             by = module]
      f_centroids[, display_label := ifelse(
        module %in% names(mod_display_labels),
        mod_display_labels[module], module
      )]

      sig_f_range <- range(f_nodes$sig_num, na.rm = TRUE)
      if (diff(sig_f_range) <= 0) sig_f_range <- c(0, 1)

      p_focus <- ggplot() +
        geom_segment(
          data = f_edges_dt,
          aes(x = from_x, y = from_y, xend = to_x, yend = to_y,
              linewidth = weight, alpha = weight),
          color = "#505050"
        ) +
        geom_polygon(
          data = f_hulls,
          aes(x = x, y = y, fill = module, group = module),
          alpha = 0.08, color = NA, show.legend = FALSE
        ) +
        geom_point(
          data = f_nodes,
          aes(x = x, y = y, fill = module, size = sig_num),
          shape = 21, color = "white", stroke = 0.4, alpha = 0.90
        ) +
        geom_point(
          data = f_nodes[degree >= 3],
          aes(x = x, y = y, size = sig_num),
          shape = 21, fill = NA, color = "#B8860B", stroke = 1.5, alpha = 0.55
        ) +
        geom_label(
          data = f_centroids,
          aes(x = x, y = y, label = display_label, fill = mod),
          color = "white", size = 2.8, fontface = "bold",
          alpha = 0.85, label.padding = unit(0.3, "lines"),
          linewidth = 0.2, show.legend = FALSE
        ) +
        geom_text_repel(
          data = f_nodes[show_label == TRUE],
          aes(x = x, y = y, label = name),
          size = 2.4, color = "grey15", fontface = "bold",
          max.overlaps = 60, box.padding = 0.3, point.padding = 0.15,
          segment.color = "grey60", segment.size = 0.2,
          show.legend = FALSE
        ) +
        scale_fill_manual(
          name = "Module",
          values = big_colors,
          labels = legend_labels,
          guide = guide_legend(order = 1, override.aes = list(size = 4, alpha = 0.85))
        ) +
        scale_size_continuous(
          name = expression(-log[10](italic(p))),
          range = c(2.5, 10),
          guide = guide_legend(order = 2)
        ) +
        scale_linewidth_continuous(
          name = "Shared annotations",
          range = c(0.3, 3.0),
          guide = guide_legend(order = 3)
        ) +
        scale_alpha_continuous(
          name = "Shared annotations",
          range = c(0.15, 0.65),
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
          legend.key.size  = unit(0.4, "cm"),
          legend.spacing.y = unit(0.1, "cm"),
          legend.box       = "vertical",
          legend.margin    = margin(0, 0, 0, 0)
        ) +
        labs(
          title = "Co-Function Network (Focused): Top Functional Modules",
          subtitle = sprintf(
            paste0("%d isoforms in %d modules  |  %d co-function edges  |  ",
                   "node size = -log10(p)"),
            nrow(f_nodes), n_big, nrow(f_edges_dt)
          )
        )

      pdf_focus <- paste0(prefix, "_network_cofunction_focused.pdf")
      png_focus <- paste0(prefix, "_network_cofunction_focused.png")
      ggsave(pdf_focus, p_focus, width = 18, height = 16, dpi = 150, device = "pdf")
      cat(sprintf("  → %s\n", pdf_focus))
      ggsave(png_focus, p_focus, width = 18, height = 16, dpi = 300, device = "png")
      cat(sprintf("  → %s\n", png_focus))
    }
  }
}

cat(sprintf("\n[DONE] Network plots saved to %s\n", prefix))

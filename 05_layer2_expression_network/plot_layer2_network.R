#!/usr/bin/env Rscript
#
# plot_layer2_network.R — Publication-quality co-expression network plot
#   Style: matches V1 plot_network.R (white bg, scientific legends, TIFF)
#
# Usage:
#   Rscript plot_layer2_network.R <edges.tsv> <ranking.tsv> [cofunc_nodes.tsv] [prefix]
#   HIGHLIGHT_GENES="LOC111886547" LABEL_DENSITY="dense" Rscript ...

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforce)
  library(ggrepel)
  library(igraph)
  library(RColorBrewer)
  library(viridisLite)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript plot_layer2_network.R <edges.tsv> <ranking.tsv> [cofunc_nodes.tsv] [prefix] [module_id] [ai_analysis.json]")
edges_file       <- args[1]
ranking_file     <- args[2]
cofunc_nodes_file <- if (length(args) >= 3 && file.exists(args[3])) args[3] else NULL
prefix           <- if (length(args) >= 4) args[4] else "layer2_network"
module_id_arg    <- if (length(args) >= 5) args[5] else "M?"
ai_json          <- if (length(args) >= 6 && file.exists(args[6])) args[6] else NULL

# Environment
highlight_genes <- character(0)
hl_env <- Sys.getenv("HIGHLIGHT_GENES")
if (nchar(hl_env) > 0) highlight_genes <- trimws(unlist(strsplit(hl_env, ",")))
label_density <- Sys.getenv("LABEL_DENSITY")
if (nchar(label_density) == 0) label_density <- "dense"

# ---- Resolve module display name from AI analysis -------------------------
module_display_name <- module_id_arg
if (!is.null(ai_json) && requireNamespace("jsonlite", quietly = TRUE)) {
  ai_data <- jsonlite::fromJSON(ai_json)
  mods <- ai_data$modules  # data.frame: rows = modules
  if (is.data.frame(mods) && nrow(mods) > 0) {
    for (i in seq_len(nrow(mods))) {
      if (mods[i, "module_id"] == module_id_arg) {
        ai_name <- mods[i, "name"]
        if (!is.null(ai_name) && nchar(ai_name) > 0) {
          module_display_name <- sprintf("%s: %s", module_id_arg, ai_name)
        }
        break
      }
    }
  }
}
# Fallback: try reading module_id from env if arg not given
if (module_id_arg == "M?") {
  env_mid <- Sys.getenv("MODULE_ID")
  if (nchar(env_mid) > 0) module_display_name <- env_mid
}
cat(sprintf("  Module label: %s\n", module_display_name))

# ---- Read data -----------------------------------------------------------
cat(sprintf("Reading: %s\n", edges_file))
edges <- fread(edges_file)
cat(sprintf("Reading: %s\n", ranking_file))
ranking <- fread(ranking_file)

# ---- GWAS seeds from co-function nodes ------------------------------------
gwas_seeds <- character(0)
if (!is.null(cofunc_nodes_file)) {
  cofunc <- fread(cofunc_nodes_file)
  seed_locs  <- cofunc$Isoform[grep("^LOC", cofunc$Isoform)]
  seed_tcons <- cofunc$Isoform[grep("^TCONS", cofunc$Isoform)]
  matched <- character(0)
  for (i in seq_len(nrow(ranking))) {
    gene <- ranking$gene_id[i]
    if (gene %in% seed_locs) { matched <- c(matched, gene); next }
    qids <- trimws(unlist(strsplit(ranking$query_ids[i], ",")))
    isos <- if ("isoforms" %in% names(ranking)) ranking$isoforms[i] else ""
    all_ids <- c(qids, if (nchar(isos) > 0) trimws(unlist(strsplit(isos, ","))) else character(0))
    if (any(all_ids %in% seed_tcons)) matched <- c(matched, gene)
  }
  gwas_seeds <- unique(c(matched, intersect(seed_locs, ranking$gene_id)))
}

# ---- Rising hubs ----------------------------------------------------------
rising_genes <- ranking[rank_change >= 15 & gwas_pvalue > 0.01, gene_id]

# ---- Build igraph ---------------------------------------------------------
setnames(edges, "source_gene", "source"); setnames(edges, "target_gene", "target")
edges$weight <- abs(edges$rho)
all_nodes <- unique(c(edges$source, edges$target))
ranked_nodes <- intersect(all_nodes, ranking$gene_id)
edges_sub <- edges[source %in% ranked_nodes & target %in% ranked_nodes]
g <- graph_from_data_frame(edges_sub[, .(source, target, weight, rho)],
                           vertices = ranked_nodes, directed = FALSE)

idx <- match(V(g)$name, ranking$gene_id)
V(g)$kME            <- ranking$kME[idx]
V(g)$combined_score <- ranking$combined_score[idx]
V(g)$rank_change    <- ranking$rank_change[idx]
V(g)$gwas_pvalue    <- ranking$gwas_pvalue[idx]
V(g)$evidence       <- ranking$evidence[idx]
V(g)$neg_log10_p    <- ranking$neg_log10_p[idx]
V(g)$query_ids      <- ranking$query_ids[idx]

V(g)$category <- "expanded"
V(g)$category[V(g)$name %in% gwas_seeds]     <- "GWAS seed"
V(g)$category[V(g)$name %in% highlight_genes] <- "Highlight"
V(g)$category[V(g)$name %in% rising_genes &
              V(g)$category == "expanded"]    <- "Rising hub"

cat_colors <- c(
  "GWAS seed"   = "#E41A1C",
  "Highlight"   = "#FF7F00",
  "Rising hub"  = "#4DAF4A",
  "expanded"    = "#B0B0B0"
)

# ---- Layout ----------------------------------------------------------------
set.seed(42)
lo <- layout_with_fr(g, weights = E(g)$weight^2, niter = 5000)

# ---- Node dataframe --------------------------------------------------------
node_df <- data.table(
  name     = V(g)$name, x = lo[, 1], y = lo[, 2],
  kME      = V(g)$kME, category = V(g)$category,
  combined = V(g)$combined_score, gwas_p = V(g)$gwas_pvalue,
  neg_log10_p = V(g)$neg_log10_p, rank_change = V(g)$rank_change,
  evidence = V(g)$evidence, query = V(g)$query_ids, degree = degree(g)
)
node_df$sig_scaled <- rescale(pmax(node_df$combined, -3), to = c(1.5, 7),
                               from = c(-3, max(node_df$combined, 1)))
node_df$size <- node_df$sig_scaled

# ---- Edge dataframe --------------------------------------------------------
el <- get.edgelist(g, names = FALSE)
edge_df <- data.table(from_idx = el[,1], to_idx = el[,2], rho = abs(E(g)$rho))
edge_df$x_from <- lo[edge_df$from_idx, 1]; edge_df$y_from <- lo[edge_df$from_idx, 2]
edge_df$x_to   <- lo[edge_df$to_idx, 1];   edge_df$y_to   <- lo[edge_df$to_idx, 2]

# ---- Labeling --------------------------------------------------------------
node_df$label <- ""
node_df$label[node_df$category %in% c("GWAS seed", "Highlight")] <-
  node_df$name[node_df$category %in% c("GWAS seed", "Highlight")]
if (label_density %in% c("normal", "dense", "all")) {
  if (length(highlight_genes) > 0) {
    hl_idx <- which(V(g)$name %in% highlight_genes)
    nbrs <- unique(unlist(adjacent_vertices(g, hl_idx)))
    nbr_names <- V(g)$name[nbrs]
    node_df$label[node_df$name %in% nbr_names & node_df$label == ""] <-
      node_df$name[node_df$name %in% nbr_names & node_df$label == ""]
  }
  n_rising <- if(label_density == "dense") 20 else 10
  rising_kme <- node_df[category == "Rising hub"][order(-kME)][1:n_rising, name]
  node_df$label[node_df$name %in% rising_kme & node_df$label == ""] <-
    node_df$name[node_df$name %in% rising_kme & node_df$label == ""]
  n_deg <- if(label_density == "dense") 10 else 5
  top_deg <- node_df[order(-degree)][1:n_deg, name]
  node_df$label[node_df$name %in% top_deg & node_df$label == ""] <-
    node_df$name[node_df$name %in% top_deg & node_df$label == ""]
}
if (label_density %in% c("dense", "all"))
  node_df$label[node_df$kME > 0.7 & node_df$label == ""] <- node_df$name[node_df$kME > 0.7 & node_df$label == ""]
if (label_density == "all")
  node_df$label[node_df$label == ""] <- node_df$name[node_df$label == ""]

n_labeled <- sum(node_df$label != "")
cat(sprintf("  Nodes: %d  Edges: %d  Labels: %d  GWAS: %d  Highlight: %d  Rising: %d  Expanded: %d\n",
            nrow(node_df), nrow(edge_df), n_labeled,
            sum(node_df$category == "GWAS seed"), sum(node_df$category == "Highlight"),
            sum(node_df$category == "Rising hub"), sum(node_df$category == "expanded")))

# ---- Theme (V1 style) -----------------------------------------------------
theme_v1 <- theme_void() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  plot.title       = element_text(size = 15, face = "bold", hjust = 0.5),
  plot.subtitle    = element_text(size = 10, hjust = 0.5, color = "grey35", margin = margin(t = 3, b = 5)),
  plot.margin      = margin(15, 15, 15, 15),
  legend.position  = "right",
  legend.title     = element_text(size = 9, face = "bold"),
  legend.text      = element_text(size = 8),
  legend.key.size  = unit(0.5, "cm"),
  legend.spacing.y = unit(0.1, "cm"),
  legend.box       = "vertical",
  legend.margin    = margin(0, 0, 0, 0)
)

# ---- Plot 1: Full co-expression network -----------------------------------
cat("Plot 1: full network...\n")

# Hulls around core module
core_nodes <- node_df[kME > 0.5]
if (nrow(core_nodes) >= 3) {
  core_nodes$tight <- core_nodes$kME > 0.7
  hull_core <- core_nodes[, .SD[chull(x, y)], by = tight]
  hull_colors <- c("FALSE" = "#4DAF4A33", "TRUE" = "#4DAF4A66")
}

p1 <- ggplot() +
  geom_segment(data = edge_df,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to, alpha = rho),
               color = "grey55", linewidth = 0.45)

if (nrow(core_nodes) >= 3 && exists("hull_core")) {
  p1 <- p1 + geom_polygon(data = hull_core,
                           aes(x = x, y = y, group = tight, color = tight),
                           fill = NA, linewidth = 0.55, linetype = "dashed") +
    scale_color_manual(values = hull_colors, guide = "none")
}

p1 <- p1 +
  geom_point(data = node_df,
             aes(x = x, y = y, size = size, fill = category),
             shape = 21, color = "grey30", stroke = 0.3) +
  geom_text_repel(data = node_df[label != ""],
                  aes(x = x, y = y, label = label),
                  size = 2.6, max.overlaps = 50, box.padding = 0.3,
                  segment.size = 0.25, segment.color = "grey50",
                  force = 1.5, min.segment.length = 0.1) +
  scale_fill_manual(name = "Category", values = cat_colors,
                    guide = guide_legend(override.aes = list(size = 3.5))) +
  scale_alpha_continuous(range = c(0.15, 0.7), guide = "none") +
  scale_size_continuous(range = c(2.2, 8), guide = "none") +
  coord_fixed() +
  theme_v1 +
  labs(title = paste(module_display_name, "Co-expression Network"),
       subtitle = sprintf("%d genes  |  %d edges  |  |r| >= 0.65, FDR < 0.05",
                          nrow(node_df), nrow(edge_df)))

pdf_full <- paste0(prefix, "_coexpression_network.pdf")
ggsave(pdf_full, p1, width = 12, height = 11, dpi = 150)
ggsave(sub("\\.pdf$", ".png", pdf_full), p1, width = 12, height = 11, dpi = 200)
ggsave(sub("\\.pdf$", ".tiff", pdf_full), p1, width = 12, height = 11, dpi = 300, compression = "lzw")
cat(sprintf("  -> %s\n", pdf_full))

# ---- Plot 2: Focused (largest component) -----------------------------------
cat("Plot 2: focused view...\n")
comps <- components(g)
lg <- which(comps$membership == which.max(comps$csize))
g_focus <- induced_subgraph(g, lg)
deg_f <- degree(g_focus); g_focus <- induced_subgraph(g_focus, names(deg_f[deg_f >= 1]))

if (vcount(g_focus) >= 3) {
  set.seed(42); lo2 <- layout_with_fr(g_focus, weights = E(g_focus)$weight^2, niter = 8000)
  el2 <- get.edgelist(g_focus, names = FALSE)
  ef2 <- data.table(from_idx = el2[,1], to_idx = el2[,2], rho = abs(E(g_focus)$rho))
  ef2$x_from <- lo2[ef2$from_idx,1]; ef2$y_from <- lo2[ef2$from_idx,2]
  ef2$x_to   <- lo2[ef2$to_idx,1];   ef2$y_to   <- lo2[ef2$to_idx,2]

  nf2 <- data.table(name = V(g_focus)$name, x = lo2[,1], y = lo2[,2],
                    kME = V(g_focus)$kME, category = V(g_focus)$category,
                    combined = V(g_focus)$combined_score,
                    degree = degree(g_focus), gwas_p = V(g_focus)$gwas_pvalue)
  nf2$sig_scaled <- rescale(pmax(nf2$combined, -3), to = c(2, 9), from = c(-3, max(nf2$combined, 1)))
  nf2$size <- nf2$sig_scaled

  nf2$label <- ""
  nf2$label[nf2$category %in% c("GWAS seed", "Highlight")] <-
    nf2$name[nf2$category %in% c("GWAS seed", "Highlight")]
  n_rising <- if(label_density == "dense") 15 else 8
  tr <- nf2[category == "Rising hub"][order(-kME)][1:n_rising, name]
  nf2$label[nf2$name %in% tr & nf2$label == ""] <- nf2$name[nf2$name %in% tr & nf2$label == ""]
  td <- nf2[order(-degree)][1:8, name]
  nf2$label[nf2$name %in% td & nf2$label == ""] <- nf2$name[nf2$name %in% td & nf2$label == ""]

  nf2$annot <- nf2$label
  key_set <- c(gwas_seeds, highlight_genes, tr[1:5])
  for (i in seq_len(nrow(nf2)))
    if (nf2$label[i] != "" && nf2$name[i] %in% key_set)
      nf2$annot[i] <- sprintf("%s  (kME=%.2f, deg=%d)", nf2$name[i], nf2$kME[i], nf2$degree[i])

  p2 <- ggplot() +
    geom_segment(data = ef2,
                 aes(x = x_from, y = y_from, xend = x_to, yend = y_to, alpha = rho),
                 color = "grey55", linewidth = 0.55) +
    geom_point(data = nf2,
               aes(x = x, y = y, size = size, fill = category),
               shape = 21, color = "grey30", stroke = 0.3) +
    geom_text_repel(data = nf2[label != ""],
                    aes(x = x, y = y, label = annot),
                    size = 2.8, max.overlaps = 60, box.padding = 0.5,
                    segment.size = 0.3, segment.color = "grey40", force = 3,
                    min.segment.length = 0.1) +
    scale_fill_manual(name = "Category", values = cat_colors,
                      guide = guide_legend(override.aes = list(size = 3.5))) +
    scale_alpha_continuous(range = c(0.15, 0.75), guide = "none") +
    scale_size_continuous(range = c(2.5, 11), guide = "none") +
    coord_fixed() +
    theme_v1 +
    labs(title = paste(module_display_name, "— Largest Connected Component"),
         subtitle = sprintf("%d genes  |  %d edges  |  ● GWAS seed  ★ highlight  ▲ rising hub",
                            nrow(nf2), nrow(ef2)))

  pdf_focus <- paste0(prefix, "_coexpression_focused.pdf")
  ggsave(pdf_focus, p2, width = 14, height = 13, dpi = 150)
  ggsave(sub("\\.pdf$", ".png", pdf_focus), p2, width = 14, height = 13, dpi = 200)
  ggsave(sub("\\.pdf$", ".tiff", pdf_focus), p2, width = 14, height = 13, dpi = 300, compression = "lzw")
  cat(sprintf("  -> %s\n", pdf_focus))
}

# ---- Plot 3: kME vs GWAS scatter -------------------------------------------
cat("Plot 3: evidence scatter...\n")
node_df$sig_label <- ""
node_df$sig_label[node_df$category != "expanded"] <- node_df$name[node_df$category != "expanded"]
node_df$sig_label[node_df$category == "Rising hub" & node_df$kME < 0.7] <- ""

max_x <- max(node_df$neg_log10_p, na.rm = TRUE)
if (max_x < 0.1) max_x <- 5

p3 <- ggplot(node_df, aes(x = neg_log10_p, y = kME)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "grey70", alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey70", alpha = 0.5) +
  annotate("text", x = -log10(0.05) + 0.05, y = 0.98, label = "p = 0.05", size = 3.5, color = "grey50", hjust = 0) +
  annotate("text", x = max_x * 0.92, y = 0.52, label = "kME = 0.5", size = 3.5, color = "grey50", hjust = 1) +
  geom_point(aes(size = combined, fill = category),
             shape = 21, color = "grey30", stroke = 0.3, alpha = 0.85) +
  geom_text_repel(data = node_df[sig_label != ""],
                  aes(label = name), size = 2.6, max.overlaps = 25,
                  box.padding = 0.3, segment.size = 0.25) +
  scale_fill_manual(name = "Category", values = cat_colors,
                    guide = guide_legend(override.aes = list(size = 3))) +
  scale_size_continuous(range = c(1.5, 7), guide = "none") +
  labs(x = expression(-log[10](p[GWAS])), y = "kME (module eigengene connectivity)",
       title = paste(module_display_name, "— GWAS vs Expression Evidence"),
       subtitle = sprintf("%d genes", nrow(node_df))) +
  theme_bw(base_size = 12) + theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

pdf_scatter <- paste0(prefix, "_evidence_scatter.pdf")
ggsave(pdf_scatter, p3, width = 10, height = 8, dpi = 150)
ggsave(sub("\\.pdf$", ".png", pdf_scatter), p3, width = 10, height = 8, dpi = 200)
cat(sprintf("  -> %s\n", pdf_scatter))

# ---- Plot 4: kME bar chart (V1 style) --------------------------------------
cat("Plot 4: kME bar chart...\n")
top20 <- node_df[order(-kME)][1:min(20, nrow(node_df))]
top20$name <- factor(top20$name, levels = rev(top20$name))

p4 <- ggplot(top20, aes(x = kME, y = name, fill = category)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(name = "Category", values = cat_colors,
                    guide = guide_legend(override.aes = list(size = 3))) +
  labs(x = "kME (module connectivity)", y = "",
       title = paste(module_display_name, "— Top 20 Genes by Module Connectivity"),
       subtitle = sprintf("%d genes ranked by kME", nrow(node_df))) +
  theme_bw(base_size = 12) + theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

pdf_bar <- paste0(prefix, "_kME_barchart.pdf")
ggsave(pdf_bar, p4, width = 10, height = 7, dpi = 150)
ggsave(sub("\\.pdf$", ".png", pdf_bar), p4, width = 10, height = 7, dpi = 200)
cat(sprintf("  -> %s\n", pdf_bar))

# ---- Summary --------------------------------------------------------------
cat(sprintf("\n  ● GWAS seed=%d  |  ★ highlight=%d  |  ▲ rising=%d  |  · expanded=%d\n",
            sum(node_df$category == "GWAS seed"), sum(node_df$category == "Highlight"),
            sum(node_df$category == "Rising hub"), sum(node_df$category == "expanded")))
cat("Done.\n")

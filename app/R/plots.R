# Plots for the TOPOSCORE app. Colour-blind-safe palette (Okabe-Ito) for classes.
library(ggplot2)

TOPO_COLOURS <- c("SIG1" = "#D55E00", "Gray" = "#999999", "SIG2" = "#0072B2",
                  "SIG1+" = "#D55E00", "SIG2+" = "#0072B2")
BRAND <- c(coral = "#F24E57", grey = "#554A48")

# ---- reference distribution (discovery cohort, Derosa et al. 2024) ------------------------------
plot_reference_density <- function(scores, clin, results = NULL, highlight = character(0),
                                   lims = TOPOSCORE_CUTOFFS) {
  ref <- merge(scores, clin, by = "Sample_id")
  ref$OS_category <- factor(ifelse(ref$OS < 12, "OS < 12 months", "OS >= 12 months"),
                            levels = c("OS < 12 months", "OS >= 12 months"))
  ymax <- max(vapply(split(ref$TOPOB01, ref$OS_category), function(x) max(density(x, from = 0, to = 1)$y), 0))

  p <- ggplot(ref, aes(x = TOPOB01, fill = OS_category)) +
    geom_density(alpha = 0.55, colour = NA) +
    geom_vline(xintercept = unname(lims), colour = "grey20", linewidth = 0.5) +
    scale_fill_manual(values = c("OS < 12 months" = "#4D4D4D", "OS >= 12 months" = "#009E73"),
                      name = sprintf("Reference cohort (n = %d)", nrow(ref))) +
    annotate("text", x = mean(c(0, lims[[1]])), y = ymax * 0.9, label = "SIG1\n68% OS < 12 months",
             colour = TOPO_COLOURS[["SIG1"]], size = 4, fontface = "bold") +
    annotate("text", x = mean(lims), y = ymax * 1.08, label = "Grey zone", colour = "grey35", size = 4) +
    annotate("text", x = mean(c(lims[[2]], 1)), y = ymax * 0.9, label = "SIG2\n21% OS < 12 months",
             colour = TOPO_COLOURS[["SIG2"]], size = 4, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), name = "S-score") +
    scale_y_continuous(name = "Density", expand = expansion(mult = c(0.02, 0.25))) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")

  if (!is.null(results) && nrow(results)) {
    p <- p + geom_rug(data = results, aes(x = S_score), inherit.aes = FALSE, sides = "b",
                      colour = BRAND[["coral"]], alpha = 0.6, length = unit(0.04, "npc"))
    hl <- results[results$Sample_id %in% highlight, , drop = FALSE]
    if (nrow(hl)) {
      p <- p + geom_vline(data = hl, aes(xintercept = S_score),
                          colour = BRAND[["coral"]], linetype = "dashed", linewidth = 0.8)
      if (nrow(hl) <= 30) {
        hl$y <- ymax * 1.2
        p <- p + ggrepel::geom_text_repel(data = hl, aes(x = S_score, y = y, label = Sample_id),
                                          inherit.aes = FALSE, colour = BRAND[["coral"]], size = 3.6,
                                          direction = "x", angle = 90, hjust = 0, segment.colour = NA,
                                          max.overlaps = Inf, seed = 1)
      }
    }
  }
  p
}

# ---- SIG species heatmap -----------------------------------------------------------------------------
# presence: samples x species logical matrix (attr(res, "presence")); res: results table
plot_sig_heatmap <- function(res, presence, crosswalk, not_in_table = character(0), samples = NULL) {
  if (!is.null(samples) && length(samples)) {
    keep <- rownames(presence) %in% samples
    presence <- presence[keep, , drop = FALSE]; res <- res[res$Sample_id %in% samples, , drop = FALSE]
  }
  ord_s <- res$Sample_id[order(res$S_score, res$Sample_id)]
  sig <- setNames(crosswalk$SIG, crosswalk$species_Jan21)
  prev <- colMeans(presence)
  sp_order <- unlist(lapply(c("SIG1", "SIG2"), function(g) {
    s <- names(sig)[sig == g]; s[order(prev[s], s)] }))

  long <- data.frame(Sample_id = rep(rownames(presence), times = ncol(presence)),
                     species = rep(colnames(presence), each = nrow(presence)),
                     detected = as.vector(presence), stringsAsFactors = FALSE)
  long$state <- ifelse(long$species %in% not_in_table, "Not in table",
                       ifelse(long$detected, "Detected", "Not detected"))
  long$group <- sig[long$species]
  cls <- data.frame(Sample_id = res$Sample_id, species = "TOPOSCORE", state = res$Toposcore,
                    group = "Call", stringsAsFactors = FALSE)
  d <- rbind(long[, names(cls)], cls)
  d$Sample_id <- factor(d$Sample_id, levels = ord_s)
  d$species <- factor(gsub("_", " ", d$species), levels = gsub("_", " ", c(sp_order, "TOPOSCORE")))
  d$group <- factor(d$group, levels = c("Call", "SIG1", "SIG2"))
  d$state <- factor(d$state, levels = c("Detected", "Not detected", "Not in table", "SIG1+", "SIG2+"))

  n <- length(ord_s)
  ggplot(d, aes(x = Sample_id, y = species, fill = state)) +
    geom_tile(colour = if (n <= 80) "white" else NA, linewidth = 0.3) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c("Detected" = BRAND[["grey"]], "Not detected" = "#EEEEEE",
                                 "Not in table" = "#F8D3D5", "SIG1+" = TOPO_COLOURS[["SIG1+"]],
                                 "SIG2+" = TOPO_COLOURS[["SIG2+"]]), drop = FALSE, name = NULL) +
    labs(x = sprintf("Samples (n = %d, ordered by S-score)", n), y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = if (n <= 80) element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8) else element_blank(),
          axis.text.y = element_text(face = "italic", size = 8.5),
          panel.grid = element_blank(), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold"), legend.position = "top")
}

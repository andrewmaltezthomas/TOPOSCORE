# TOPOSCORE calculation (Derosa et al., Cell 2024) for profiles from several profilers/databases.
#
#   res <- calculate_toposcore("profile.tsv", "GTDB_r220")
#   res                         # one row per sample (same columns as before, plus Akk_abundance)
#   attr(res, "presence")       # samples x 82 logical matrix: which SIG species were detected
#   attr(res, "diagnostics")    # database, SIG species absent from the table, Akkermansia source, ...
#
# Differences from the previous implementation (results are otherwise identical):
#   * one code path for every database: format-specific parsing, then a single scoring step
#   * the crosswalk is passed in (default: data/toposcore_crosswalk.tsv) rather than read from fixed paths
#   * non-abundance columns (NCBI_tax_id, bracken_taxid, genome_size, ...) are dropped, never read as samples
#   * works for a single sample; MetaPhlAn species-level input no longer fails when SGB9226 is present
#   * Meteor2 works when msp_0025 is absent, and uses the Meteor2 Akkermansia cut-off (4.52%)
#   * SIG species with no alias anywhere in the table are reported (diagnostics), not only warned about

TOPOSCORE_CUTOFFS <- c(SIG1 = 0.5351, SIG2 = 0.7911)   # S-score class boundaries
TOPOSCORE_N_SIG   <- c(SIG1 = 37, SIG2 = 45)           # published signature sizes (S-score denominators)

# One row per supported input. akk_id: exact taxon id used for Akkermansia; akk_fallback: used if absent.
toposcore_databases <- data.frame(
  id           = c("species_Jan21", "SGB_Jan21", "SGB_Jun23", "SGB_Jan25", "GTDB_r207", "GTDB_r220", "UHGG_v2_0_2", "MSP_Meteor2"),
  label        = c("Species - MetaPhlAn 4 vJan21", "SGB - MetaPhlAn 4 vJan21", "SGB - MetaPhlAn 4 vJun23",
                   "SGB - MetaPhlAn 4 vJan25", "Species - GTDB r207", "Species - GTDB r220",
                   "Species - UHGG v2.0.2 (Kraken2/Bracken)", "MSP - Meteor2"),
  format       = c("metaphlan_species", rep("metaphlan_sgb", 3), rep("gtdb", 3), "meteor"),
  akk_id       = c("SGB9226", "SGB9226", "SGB9226", "SGB9226", rep("s__Akkermansia muciniphila", 3), "msp_0025"),
  akk_fallback = c("Akkermansia_muciniphila", rep("s__Akkermansia_muciniphila", 3), rep(NA, 3), "muciniphila$"),
  akk_cutoff   = c(rep(4.799, 7), 4.52),
  stringsAsFactors = FALSE
)

load_toposcore_crosswalk <- function(path = file.path("data", "toposcore_crosswalk.tsv")) {
  xw <- read.delim(path, colClasses = "character", check.names = FALSE)
  need <- c("species_Jan21", "SIG")
  if (!all(need %in% names(xw))) stop("Crosswalk must contain columns: ", paste(need, collapse = ", "))
  n <- table(factor(xw$SIG, levels = names(TOPOSCORE_N_SIG)))
  if (!all(n == TOPOSCORE_N_SIG)) {
    warning(sprintf("Crosswalk has %d SIG1 / %d SIG2 species; the S-score uses the published %d / %d.",
                    n[["SIG1"]], n[["SIG2"]], TOPOSCORE_N_SIG[["SIG1"]], TOPOSCORE_N_SIG[["SIG2"]]))
  }
  xw
}

# ---- reading ------------------------------------------------------------------------------------
# Returns list(m = samples x taxa numeric matrix, lineage = taxon names as given, extra = extra columns)
read_profile_table <- function(input, id_col) {
  if (is.data.frame(input)) {
    df <- as.data.frame(input, check.names = FALSE)
  } else {
    lines <- readLines(input, warn = FALSE)
    # keep the header even when written as "#clade_name" (single MetaPhlAn profiles); drop other comments
    hdr <- grep(paste0("^#?", id_col, "\t"), lines)[1]
    if (is.na(hdr)) stop("No '", id_col, "' column found. Is this the right input type?")
    lines <- c(sub("^#", "", lines[hdr]), grep("^#", lines[-seq_len(hdr)], value = TRUE, invert = TRUE))
    df <- read.delim(text = lines, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  }
  if (!id_col %in% names(df)) stop("No '", id_col, "' column found. Is this the right input type?")
  annot <- c(id_col, "NCBI_tax_id", "clade_taxid", "additional_species", "species", "bracken_taxid",
             "genome_size", "taxonomy_id", "msp_taxonomy")
  sample_cols <- setdiff(names(df), annot)
  not_num <- sample_cols[!vapply(df[sample_cols], is.numeric, logical(1))]
  if (length(not_num)) stop("Non-numeric sample column(s): ", paste(head(not_num, 5), collapse = ", "))
  if (!length(sample_cols)) stop("No sample columns found.")
  if (anyDuplicated(sample_cols)) stop("Duplicated sample names: ", paste(unique(sample_cols[duplicated(sample_cols)]), collapse = ", "))
  m <- t(as.matrix(df[sample_cols]))
  m[is.na(m)] <- 0
  dimnames(m) <- list(sample_cols, as.character(df[[id_col]]))
  list(m = m, extra = df[intersect(names(df), annot)])
}

sum_duplicate_taxa <- function(m) {
  if (!anyDuplicated(colnames(m))) return(m)
  t(rowsum(t(m), group = colnames(m), reorder = FALSE))
}

# Converts the profile to samples x taxon-id matrix (ids as used in the crosswalk) plus Akkermansia abundance
parse_profile <- function(input, db, abundance_scale = c("auto", "percent", "fraction")) {
  abundance_scale <- match.arg(abundance_scale)
  fmt <- db$format
  p <- read_profile_table(input, if (fmt == "meteor") "msp" else "clade_name")
  m <- p$m; lineage <- colnames(m)
  level_of <- function(sep) lengths(strsplit(lineage, sep, fixed = TRUE))
  last_of  <- function(sep) vapply(strsplit(lineage, sep, fixed = TRUE), function(x) x[length(x)], "")
  akk_col  <- function(mat, id)   # keep sample names, also for a single sample
    if (!is.na(id) && id %in% colnames(mat)) setNames(as.numeric(mat[, id]), rownames(mat)) else NULL

  if (fmt %in% c("metaphlan_sgb", "metaphlan_species")) {
    if (!any(grepl("|", lineage, fixed = TRUE))) stop("Taxonomy is not in MetaPhlAn format (k__|p__|...).")
    lev <- level_of("|"); last <- last_of("|")
    sgb <- m[, lev == 8, drop = FALSE]; colnames(sgb) <- sub("^t__", "", last[lev == 8])
    spp <- m[, lev == 7, drop = FALSE]; colnames(spp) <- last[lev == 7]
    akk <- akk_col(sgb, db$akk_id)
    akk_source <- db$akk_id
    if (is.null(akk)) { akk <- akk_col(spp, sub("^(s__)?", "s__", db$akk_fallback)); akk_source <- db$akk_fallback }
    if (fmt == "metaphlan_sgb") {
      if (!ncol(sgb)) stop("No SGB-level rows (t__SGB...) in this table. Choose the species-level MetaPhlAn option, ",
                           "or merge MetaPhlAn profiles that include SGBs.")
      dup <- unique(colnames(sgb)[duplicated(colnames(sgb))])
      if (length(dup)) stop("Duplicated SGB ids (check the profile merging): ", paste(dup, collapse = ", "))
      ids <- sgb
    } else {
      colnames(spp) <- sub("^s__", "", colnames(spp))
      ids <- sum_duplicate_taxa(spp)
    }
  } else if (fmt == "gtdb") {
    if (!any(grepl(";", lineage, fixed = TRUE))) stop("Taxonomy is not in GTDB format (d__;p__;...).")
    lev <- level_of(";")
    if (!any(lev == 7)) stop("No species-level rows (d__;...;s__) in this table.")
    ids <- m[, lev == 7, drop = FALSE]; colnames(ids) <- last_of(";")[lev == 7]
    ids <- sum_duplicate_taxa(ids)
    akk <- akk_col(ids, db$akk_id); akk_source <- db$akk_id
  } else if (fmt == "meteor") {
    ids <- sum_duplicate_taxa(m)
    akk <- akk_col(ids, db$akk_id); akk_source <- db$akk_id
    if (is.null(akk) && "msp_taxonomy" %in% names(p$extra)) {
      hit <- as.character(p$extra$msp[grepl(db$akk_fallback, p$extra$msp_taxonomy)])
      if (length(hit)) { akk <- rowSums(ids[, hit, drop = FALSE]); akk_source <- paste(hit, collapse = "+") }
    }
  } else stop("Unknown format: ", fmt)

  if (is.null(akk)) { akk <- setNames(rep(0, nrow(m)), rownames(m)); akk_source <- "not in table (counted as 0)" }

  # Abundance units: MetaPhlAn is always %, GTDB/Meteor tables may hold fractions
  to_percent <- switch(abundance_scale, percent = FALSE, fraction = TRUE,
                       auto = fmt %in% c("gtdb", "meteor") && max(ids) <= 1.05 && all(rowSums(ids) <= 1.05))
  if (to_percent) { ids <- ids * 100; akk <- akk * 100 }
  list(m = ids, akk = akk, akk_source = akk_source, scaled_to_percent = to_percent)
}

# Guess the input type from the first lines, to catch a wrong database choice before scoring
detect_profile_format <- function(path, n = 400) {
  l <- readLines(path, n = n, warn = FALSE)
  l <- l[!grepl("^#", l) | grepl("^#clade_name", l)]
  if (length(l) && grepl("^msp\t", l[1])) return("meteor")
  ids <- sub("\t.*", "", l[-1])
  if (any(grepl("|", ids, fixed = TRUE)) && any(grepl("^k__", ids))) return("metaphlan")
  if (any(grepl(";", ids, fixed = TRUE)) && any(grepl("^d__", ids))) return("gtdb")
  NA_character_
}

# ---- scoring ------------------------------------------------------------------------------------
calculate_toposcore <- function(input, taxonomy_database, crosswalk = load_toposcore_crosswalk(),
                                abundance_scale = c("auto", "percent", "fraction"),
                                cutoffs = TOPOSCORE_CUTOFFS) {
  db <- toposcore_databases[toposcore_databases$id == taxonomy_database, ]
  if (nrow(db) != 1) stop("taxonomy_database must be one of: ", paste(toposcore_databases$id, collapse = ", "))
  if (!taxonomy_database %in% names(crosswalk)) stop("Crosswalk has no column '", taxonomy_database, "'")

  prof <- parse_profile(input, db, match.arg(abundance_scale))
  m <- prof$m

  aliases <- lapply(strsplit(crosswalk[[taxonomy_database]], ","), function(a) {
    a <- trimws(a); a[!is.na(a) & a != "" & a != "-"]
  })
  in_table <- vapply(aliases, function(a) any(a %in% colnames(m)), logical(1))
  presence <- matrix(vapply(aliases, function(a) {
    a <- intersect(a, colnames(m))
    if (length(a)) rowSums(m[, a, drop = FALSE] > 0) > 0 else rep(FALSE, nrow(m))
  }, logical(nrow(m))), nrow = nrow(m), dimnames = list(rownames(m), crosswalk$species_Jan21))

  is1 <- crosswalk$SIG == "SIG1"; is2 <- crosswalk$SIG == "SIG2"
  taxa_str <- function(cols) vapply(seq_len(nrow(presence)), function(i)
    paste(colnames(presence)[cols][presence[i, cols]], collapse = ";"), "")

  sig1 <- rowSums(presence[, is1, drop = FALSE]); sig2 <- rowSums(presence[, is2, drop = FALSE])
  S <- round((sig2 / TOPOSCORE_N_SIG[["SIG2"]] - sig1 / TOPOSCORE_N_SIG[["SIG1"]] + 1) / 2, 3)
  sig_class <- ifelse(S <= cutoffs[["SIG1"]], "SIG1", ifelse(S >= cutoffs[["SIG2"]], "SIG2", "Gray"))
  akk <- unname(prof$akk[rownames(m)])
  akk_status <- ifelse(akk >= db$akk_cutoff, "High", ifelse(akk == 0, "Zero", "Normal"))
  topo <- ifelse(sig_class == "SIG1", "SIG1+", ifelse(sig_class == "SIG2", "SIG2+",
                 ifelse(akk_status %in% c("High", "Zero"), "SIG1+", "SIG2+")))

  res <- data.frame(Sample_id = rownames(m), Akk_status = akk_status, SIG1_count = unname(sig1),
                    SIG2_count = unname(sig2), S_score = unname(S), SIG_class = sig_class, Toposcore = topo,
                    SIG1_taxa = taxa_str(is1), SIG2_taxa = taxa_str(is2), Akk_abundance = akk,
                    stringsAsFactors = FALSE)
  res <- tibble::as_tibble(res)
  attr(res, "presence") <- presence
  attr(res, "diagnostics") <- list(
    database = db$label, n_samples = nrow(m),
    sig1_not_in_table = crosswalk$species_Jan21[is1 & !in_table],
    sig2_not_in_table = crosswalk$species_Jan21[is2 & !in_table],
    akk_source = prof$akk_source, akk_cutoff = db$akk_cutoff, scaled_to_percent = prof$scaled_to_percent)
  res
}

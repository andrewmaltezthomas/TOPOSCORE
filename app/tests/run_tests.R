# Self-contained checks for the TOPOSCORE function and app.
# Run from the app folder:  Rscript tests/run_tests.R
app_dir <- if (file.exists("app.R")) "." else ".."
old_wd <- setwd(app_dir); on.exit(setwd(old_wd), add = TRUE)
suppressMessages({ library(shiny); source("R/toposcore_calc.R"); source("R/plots.R") })

n_fail <- 0
ok <- function(label, cond) {
  cat(if (isTRUE(cond)) "PASS" else { n_fail <<- n_fail + 1; "FAIL" }, "-", label, "\n")
}
quiet <- function(x) suppressWarnings(suppressMessages(x))
td <- tempfile("toposcore_tests_"); dir.create(td)
xw <- load_toposcore_crosswalk()
mpa <- "data/mpa4_vJan21_test_data.tsv"; syl <- "data/sylph_GTDBr220_test_data.tsv"

# ---- crosswalk and formats ---------------------------------------------------------------------
ok("crosswalk has 37 SIG1 and 45 SIG2 species", all(table(xw$SIG) == c(37, 45)))
ok("every database option has a crosswalk column", all(toposcore_databases$id %in% names(xw)))
ok("format detection: MetaPhlAn", identical(detect_profile_format(mpa), "metaphlan"))
ok("format detection: GTDB", identical(detect_profile_format(syl), "gtdb"))

# ---- scores on the bundled examples ------------------------------------------------------------
check_scores <- function(r) {
  s <- round(0.5 * (r$SIG2_count / 45 - r$SIG1_count / 37 + 1), 3)   # classified after rounding, as in the original app
  cls <- ifelse(s <= TOPOSCORE_CUTOFFS[["SIG1"]], "SIG1", ifelse(s >= TOPOSCORE_CUTOFFS[["SIG2"]], "SIG2", "Gray"))
  topo <- ifelse(cls == "SIG1", "SIG1+", ifelse(cls == "SIG2", "SIG2+",
                 ifelse(r$Akk_status %in% c("Zero", "High"), "SIG1+", "SIG2+")))
  isTRUE(all.equal(r$S_score, s, tolerance = 1e-6)) && identical(r$SIG_class, cls) && identical(r$Toposcore, topo)
}
for (db in c("species_Jan21", "SGB_Jan21")) {
  r <- quiet(calculate_toposcore(mpa, db, xw))
  ok(sprintf("MetaPhlAn %s: 5 samples, consistent S-score and classes", db), nrow(r) == 5 && check_scores(r))
}
r_g <- quiet(calculate_toposcore(syl, "GTDB_r220", xw))
ok("GTDB r220: 5 samples, consistent S-score and classes", nrow(r_g) == 5 && check_scores(r_g))
ok("presence matrix matches the counts",
   all(rowSums(attr(r_g, "presence")[, xw$species_Jan21[xw$SIG == "SIG1"]]) == r_g$SIG1_count))

# ---- fractions are scaled to percent -----------------------------------------------------------
x <- read.delim(syl, check.names = FALSE); x[-1] <- x[-1] / 100
write.table(x, file.path(td, "frac.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
r_f <- quiet(calculate_toposcore(file.path(td, "frac.tsv"), "GTDB_r220", xw))
ok("fractions give the same calls as percentages",
   identical(r_f$Toposcore, r_g$Toposcore) && isTRUE(attr(r_f, "diagnostics")$scaled_to_percent))

# ---- single sample (crashed in the previous version) -------------------------------------------
x1 <- read.delim(syl, check.names = FALSE)[, 1:2]
write.table(x1, file.path(td, "one.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
r1 <- quiet(calculate_toposcore(file.path(td, "one.tsv"), "GTDB_r220", xw))
ok("single-sample table works and keeps its name", nrow(r1) == 1 && r1$Sample_id == names(x1)[2])

# ---- Meteor2 without the reference Akkermansia MSP ---------------------------------------------
m_ids <- strsplit(xw$MSP_Meteor2[xw$SIG == "SIG2"][1:3], ",")
mt <- data.frame(msp = c(sapply(m_ids, `[`, 1), "msp_0999"),
                 msp_taxonomy = c("a", "b", "c", "k__Bacteria;s__Akkermansia muciniphila"),
                 S1 = c(.5, .4, .05, .05), S2 = c(.98, 0, .02, 0))
write.table(mt, file.path(td, "meteor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
r_m <- quiet(calculate_toposcore(file.path(td, "meteor.tsv"), "MSP_Meteor2", xw))
ok("Meteor2: Akkermansia found by taxonomy when msp_0025 is absent",
   attr(r_m, "diagnostics")$akk_source == "msp_0999" && identical(r_m$Akk_status, c("High", "Zero")))

# ---- wrong input is refused with a clear error -------------------------------------------------
ok("GTDB table with a MetaPhlAn SGB option gives an error",
   inherits(try(quiet(calculate_toposcore(syl, "SGB_Jan21", xw)), silent = TRUE), "try-error"))

# ---- plots build -------------------------------------------------------------------------------
ref_s <- read.csv("data/scores_disc.csv", check.names = FALSE); ref_c <- read.csv("data/clin_disc.csv", check.names = FALSE)
pdf(NULL)
ok("score distribution plot builds",
   inherits(ggplot2::ggplot_build(plot_reference_density(ref_s, ref_c, r_g, r_g$Sample_id[1:2])), "ggplot_built"))
d <- attr(r_g, "diagnostics")
ok("signature heatmap builds", inherits(ggplot2::ggplot_build(
  plot_sig_heatmap(r_g, attr(r_g, "presence"), xw, c(d$sig1_not_in_table, d$sig2_not_in_table))), "ggplot_built"))

# ---- app server ---------------------------------------------------------------------------------
app_env <- new.env(); quiet(sys.source("app.R", envir = app_env, keep.source = FALSE))
testServer(app_env$server, {
  session$setInputs(file = list(datapath = syl, name = "x.tsv"), database = "SGB_Jan21"); session$setInputs(calculate = 1)
  ok("app refuses a database that does not match the file", is.null(results()))
  session$setInputs(database = "GTDB_r220"); session$setInputs(calculate = 2)
  ok("app calculates on the example file", identical(results()$Toposcore, r_g$Toposcore))
  ok("app CSV download has all samples", nrow(read.csv(output$dl_results)) == 5)
})

cat(sprintf("\n%d failure(s)\n", n_fail)); if (n_fail > 0) quit(status = 1)

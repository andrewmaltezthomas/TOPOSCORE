# Lists the packages the TOPOSCORE app needs and reports any that are missing (does not install them).
pkgs <- c("shiny", "bslib", "DT", "ggplot2", "ggrepel", "tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Missing packages: ", paste(missing, collapse = ", "),
          "\nInstall with: install.packages(c(", paste0('"', missing, '"', collapse = ", "), "))")
} else {
  message("All TOPOSCORE app packages are installed.")
}
if (packageVersion("bslib") < "0.6.0") message("bslib >= 0.6.0 is recommended.")

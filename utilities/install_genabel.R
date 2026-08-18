#!/usr/bin/env Rscript
# =============================================================================
# install_genabel.R — install GenABEL + GenABEL.data from the bundled tarballs
# -----------------------------------------------------------------------------
# GenABEL was removed from CRAN, so install.packages("GenABEL") fails. The
# source tarballs are bundled in resources/r_packages/ so the installation is
# guaranteed to work offline / without hunting the CRAN archive.
#
# Requirements:
#   - a working R toolchain (gcc/g++/gfortran; on Windows: Rtools)
#   - GenABEL.data is installed FIRST (GenABEL depends on it)
#
# Usage:
#   Rscript utilities/install_genabel.R
#   # or point at a different directory containing the two tarballs:
#   Rscript utilities/install_genabel.R /path/to/tarballs
#
# Tip: the conda route (environment.yml, r-genabel from conda-forge) avoids
# compilation entirely and is the recommended option when conda is available.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
    pkg_dir <- args[1]
} else {
    # default: <repo>/resources/r_packages  (script lives in <repo>/utilities)
    script_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- script_args[grep("^--file=", script_args)]
    script_dir <- if (length(file_arg) >= 1) {
        dirname(sub("^--file=", "", file_arg[1]))
    } else {
        getwd()
    }
    pkg_dir <- file.path(script_dir, "..", "resources", "r_packages")
}

tarballs <- c(
    file.path(pkg_dir, "GenABEL.data_1.0.0.tar.gz"),
    file.path(pkg_dir, "GenABEL_1.8-0.tar.gz")
)

for (tarball in tarballs) {
    if (!file.exists(tarball)) {
        stop("tarball not found: ", tarball)
    }
    cat("Installing:", tarball, "\n")
    install.packages(tarball, repos = NULL, type = "source")
}

if (requireNamespace("GenABEL", quietly = TRUE)) {
    cat("\nGenABEL installed successfully.\n")
} else {
    cat("\nWARNING: GenABEL did not load after installation — check the build log.\n")
}

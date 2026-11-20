# Fit each profile with the same catalog, seed, and sampler settings as the accuracy panel.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L)
suppressPackageStartupMessages(library(signeR))
stopifnot(as.character(packageVersion("signeR")) == "2.14.0")
read_matrix <- function(path) {
    as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
}
counts <- read_matrix(args[1])
reference <- read_matrix(args[2])
background <- read_matrix(args[3])
seeds <- read_matrix(args[4])
method <- args[5]
stopifnot(identical(rownames(counts), rownames(reference)),
          identical(rownames(reference), rownames(background)),
          identical(colnames(counts), rownames(seeds)), method %in% colnames(seeds))
design <- if (method == "signeR+background") cbind(reference, background) else reference
weights <- matrix(0, ncol(counts), ncol(reference) + 1L,
                  dimnames = list(colnames(counts), c(colnames(reference), "Background")))
for (sample in colnames(counts)) {
    set.seed(seeds[sample, method])
    fit <- signeR(M = t(counts[, sample, drop = FALSE]), P = design, fixedP = TRUE,
                  main_burn = 10000, main_eval = 2000, estimate_hyper = FALSE,
                  parallelization = "sequential")
    stopifnot(identical(fit$SignExposures@signames, colnames(design)))
    exposures <- as.numeric(Median_exp(fit$SignExposures)[, 1])
    stopifnot(length(exposures) == ncol(design), all(is.finite(exposures)),
              all(exposures >= 0), sum(exposures) > 0)
    known <- seq_len(ncol(reference))
    fraction <- exposures / sum(exposures)
    weights[sample, ] <- c(fraction[known], sum(fraction[-known]))
    message("Completed ", sample, " at ", format(Sys.time(), tz = "UTC", usetz = TRUE))
}
write.table(weights, args[6], sep = "\t", quote = FALSE, col.names = NA)

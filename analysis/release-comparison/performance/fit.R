# Source the upstream fixed-catalog fitter and its numerical helpers.
args <- commandArgs(trailingOnly = TRUE)
for (name in c("internal.scripts.R", "golden_section_search.R", "whichSignatures.R")) {
    source(file.path(args[1], "R", name))
}
read_matrix <- function(path) {
    as.data.frame(t(read.delim(path, row.names = 1, check.names = FALSE)))
}
counts <- read_matrix(args[2])
reference <- read_matrix(args[3])
stopifnot(identical(names(counts), names(reference)))
profiles <- counts / rowSums(counts)
fits <- lapply(rownames(profiles), function(sample) {
    fit <- whichSignatures(profiles, sample, reference)
    c(unlist(fit$weights, use.names = FALSE), fit$unknown)
})
weights <- do.call(rbind, fits)
colnames(weights) <- c(rownames(reference), "Background")
rownames(weights) <- rownames(counts)
write.table(weights, args[4], sep = "\t", quote = FALSE, col.names = NA)

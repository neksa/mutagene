# Exchange labelled SBS96 count matrices with the MutaGene CLI.
write_mutagene_matrix <- function(counts, path) {
    counts <- as.matrix(counts)
    if (!is.numeric(counts) || nrow(counts) != 96L || is.null(rownames(counts)) ||
        is.null(colnames(counts)) || anyDuplicated(rownames(counts)) ||
        anyDuplicated(colnames(counts)) || any(!is.finite(counts)) ||
        any(counts < 0) || any(counts != floor(counts))) {
        stop("Provide a labelled, finite, nonnegative integer SBS96 matrix")
    }
    frame <- data.frame(Type=rownames(counts), counts, check.names=FALSE)
    write.table(frame, path, sep="\t", quote=FALSE, row.names=FALSE)
    invisible(path)
}

fit_mutagene <- function(counts, output, signatures="COSMICv2", executable="mutagene") {
    input <- tempfile(fileext=".tsv")
    on.exit(unlink(input), add=TRUE)
    write_mutagene_matrix(counts, input)
    status <- system2(executable, c("signature", "-f", "MATRIX", "-i", shQuote(input),
                                   "-s", shQuote(signatures), "-o", shQuote(output)))
    if (status != 0L) stop("MutaGene signature fitting fails")
    read.delim(output, check.names=FALSE, stringsAsFactors=FALSE)
}

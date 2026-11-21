# Match pentanucleotide counts to the original SBS96 cohort by specimen ID.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
load(args[1])
tri <- read.delim(args[2], row.names = 1, check.names = FALSE)
ids <- sub(".*::", "", names(tri))
source_ids <- sub(".*__", "", names(lego1536.SNV))
stopifnot(!anyDuplicated(ids), all(vapply(ids, function(id) sum(source_ids == id) == 1, logical(1))))
counts <- lego1536.SNV[, match(ids, source_ids)]
names(counts) <- names(tri)
k <- rownames(counts)
row.names(counts) <- paste0(substr(k, 3, 4), "[", substr(k, 1, 1), ">", substr(k, 2, 2), "]", substr(k, 5, 6))
stopifnot(nrow(counts) == 1536, ncol(counts) == 146)
collapsed <- rowsum(counts, substr(row.names(counts), 2, 8))
stopifnot(all(collapsed[row.names(tri), ] == as.matrix(tri)))
write.table(data.frame(Type = row.names(counts), counts, check.names = FALSE),
            args[3], sep = "\t", row.names = FALSE, quote = FALSE)

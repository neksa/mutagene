# Fit the pinned signeR package without changing its inference implementation.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)
root <- normalizePath(args[1])
.libPaths(c(normalizePath(args[2]), .libPaths()))
workers <- as.integer(args[3])
suppressPackageStartupMessages(library(signeR))
suppressPackageStartupMessages(library(posterior))
stopifnot(as.character(packageVersion('signeR')) == '2.14.0')
tasks <- read.delim(file.path(root, 'tasks.tsv'), check.names=FALSE)
reference <- as.matrix(read.delim(file.path(root, 'reference.tsv'), row.names=1, check.names=FALSE))
background <- as.matrix(read.delim(file.path(root, 'background.tsv'), row.names=1, check.names=FALSE))
catalogs <- as.matrix(read.delim(file.path(root, 'catalogs.tsv'), row.names=1, check.names=FALSE))
stopifnot(identical(rownames(reference), rownames(catalogs)))
stopifnot(identical(rownames(reference), rownames(background)))
stopifnot(max(abs(colSums(reference)-1)) < 1e-12)
writeLines(capture.output(sessionInfo()), file.path(root, 'session-info.txt'))
write.csv(installed.packages()[,c('Package','Version','LibPath')], file.path(root, 'r-packages.csv'), row.names=FALSE)
writeLines(capture.output(dput(formals(signeR))), file.path(root, 'signeR-defaults.txt'))
writeLines(deparse(body(signeR)), file.path(root, 'signeR-body.R'))
fit_one <- function(index) {
    task <- tasks[index,]
    prefix <- file.path(root, 'fits', sprintf('%04d', task$task_id))
    if(file.exists(paste0(prefix, '.json'))) stop('Refusing to overwrite a completed fit')
    begin <- proc.time()[['elapsed']]
    warnings <- character()
    log <- file(paste0(prefix, '.log'), open='wt')
    sink(log)
    result <- tryCatch(withCallingHandlers({
        design <- if(task$background) cbind(reference, background) else reference
        counts <- catalogs[,task$sample,drop=FALSE]
        set.seed(task$seed)
        fit <- signeR(M=t(counts), P=design, fixedP=TRUE, main_burn=10000,
                      main_eval=2000, estimate_hyper=FALSE, parallelization='sequential')
        exposures <- as.numeric(Median_exp(fit$SignExposures)[,1])
        stopifnot(length(exposures)==ncol(design), all(is.finite(exposures)),
                  all(exposures>=0), sum(exposures)>0)
        stopifnot(identical(fit$SignExposures@signames, colnames(design)))
        weights <- exposures/sum(exposures)
        draws <- t(fit$SignExposures@Exp[,1,])
        stopifnot(all(is.finite(draws)), all(rowSums(draws)>0))
        proportions <- draws/rowSums(draws)
        colnames(proportions) <- colnames(design)
        diagnostic <- summarise_draws(as_draws_matrix(proportions), 'rhat', 'ess_bulk', 'ess_tail', 'mcse_mean')
        active <- weights*sum(counts)>=10
        active_diagnostic <- diagnostic[active,,drop=FALSE]
        if(task$diagnostic) saveRDS(proportions, paste0(prefix,'.draws.rds'), compress='xz')
        write.table(data.frame(signature=colnames(design), median_count=exposures, weight=weights),
                    paste0(prefix,'.estimates.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
        write.table(diagnostic, paste0(prefix,'.diagnostics.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
        list(error=NULL, weights=weights, raw_total=sum(exposures), raw_total_ratio=sum(exposures)/sum(counts),
             single_chain_split_rhat_max=max(diagnostic$rhat,na.rm=TRUE),
             active_split_rhat_max=if(any(active)) max(active_diagnostic$rhat,na.rm=TRUE) else NA_real_,
             active_ess_bulk_min=if(any(active)) min(active_diagnostic$ess_bulk,na.rm=TRUE) else NA_real_,
             mcse_mean_max=max(diagnostic$mcse_mean,na.rm=TRUE),
             known_count=ncol(reference), posterior_draws=2000L)
    }, warning=function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
    }), error=function(e) list(error=conditionMessage(e)))
    sink(); close(log)
    result$seconds <- proc.time()[['elapsed']] - begin
    result$warnings <- warnings
    result$task_id <- task$task_id
    result$sample <- task$sample
    result$background <- task$background
    result$chain <- task$chain
    result$seed <- task$seed
    jsonlite::write_json(result, paste0(prefix,'.json'), auto_unbox=TRUE, pretty=TRUE, null='null', digits=17)
    cat(sprintf('Completed %d/%d: %s, background=%s, chain=%d, %.1fs, error=%s\n',
                index,nrow(tasks),task$sample,task$background,task$chain,result$seconds,
                if(is.null(result$error)) 'none' else result$error), file=stdout())
    flush.console()
    result$task_id
}
completed <- parallel::mclapply(seq_len(nrow(tasks)), fit_one, mc.cores=workers,
                                mc.preschedule=FALSE, mc.set.seed=FALSE)
stopifnot(!any(vapply(completed, inherits, logical(1), 'try-error')))
cat('Finished', length(completed), 'signeR fits\n')

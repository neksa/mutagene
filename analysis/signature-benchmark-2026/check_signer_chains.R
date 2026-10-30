# Compare independent chains without filtering benchmark estimates by convergence.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
root <- normalizePath(args[1])
.libPaths(c(normalizePath(args[2]), .libPaths()))
suppressPackageStartupMessages(library(posterior))
tasks <- read.delim(file.path(root,'tasks.tsv'),check.names=FALSE)
catalogs <- as.matrix(read.delim(file.path(root,'catalogs.tsv'),row.names=1,check.names=FALSE))
repeats <- tasks[tasks$chain==2,]
rows <- list()
for(i in seq_len(nrow(repeats))) {
    repeat_task <- repeats[i,]
    first <- tasks[tasks$sample==repeat_task$sample & tasks$background==repeat_task$background & tasks$chain==1,]
    stopifnot(nrow(first)==1)
    prefix1 <- file.path(root,'fits',sprintf('%04d',first$task_id))
    prefix2 <- file.path(root,'fits',sprintf('%04d',repeat_task$task_id))
    a <- readRDS(paste0(prefix1,'.draws.rds'))
    b <- readRDS(paste0(prefix2,'.draws.rds'))
    stopifnot(identical(dim(a),dim(b)),identical(colnames(a),colnames(b)))
    chains <- array(0,dim=c(nrow(a),2,ncol(a)),dimnames=list(NULL,NULL,colnames(a)))
    chains[,1,] <- a
    chains[,2,] <- b
    diagnostic <- summarise_draws(as_draws_array(chains),'rhat','ess_bulk','ess_tail','mcse_mean')
    path <- file.path(root,sprintf('chain-diagnostics-%04d.tsv',first$task_id))
    write.table(diagnostic,path,sep='\t',row.names=FALSE,quote=FALSE)
    m1 <- jsonlite::read_json(paste0(prefix1,'.json'),simplifyVector=TRUE)
    m2 <- jsonlite::read_json(paste0(prefix2,'.json'),simplifyVector=TRUE)
    active <- (m1$weights+m2$weights)/2*sum(catalogs[,repeat_task$sample])>=10
    rows[[i]] <- data.frame(sample=repeat_task$sample,background=as.logical(repeat_task$background),
        first_task=first$task_id,repeat_task=repeat_task$task_id,
        exposure_l1=sum(abs(m1$weights-m2$weights)),
        rhat_max=max(diagnostic$rhat,na.rm=TRUE),
        active_rhat_max=if(any(active)) max(diagnostic$rhat[active],na.rm=TRUE) else NA_real_,
        active_ess_bulk_min=if(any(active)) min(diagnostic$ess_bulk[active],na.rm=TRUE) else NA_real_)
}
write.table(do.call(rbind,rows),file.path(root,'chain-stability.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('Compared',nrow(repeats),'independent chain pairs\n')

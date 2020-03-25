
library(coxme)
library(data.table)
library(parallel)
extract_coxme_table <- function (mod){
    beta <- mod$coefficients #$fixed is not needed
    nvar <- length(beta)
    nfrail <- nrow(mod$var) - nvar
    se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
    z<- round(beta/se, 2)
    p<- signif(1 - pchisq((beta/se)^2, 1), 2)
    table=data.frame(cbind(beta,se,z,p))
    return(table)
}
eval.coxme = function(dat, dtx){
    tryCatch(
        {
            dtx$col = dat
            dtx = dtx[!is.na(col)]
            dtx = dtx[dataset %in% names(which(table(dtx$dataset) > 10))] 
            aa =  coxme(Surv(Surival, Event) ~ col + (1|dataset), dtx)
            extract_coxme_table(aa)
        },error=  function(e) rep(NA,4))
}
pfs.response.genes = mclapply(seq(ncol(expression.merge.match)), function(tt) eval.coxme(expression.merge.match[,tt], dtx=pfs.response.dt), mc.cores=45)

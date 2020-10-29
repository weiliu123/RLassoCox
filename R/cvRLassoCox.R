cvRLassoCox <- function(x, y, globalGraph=NULL, nfolds = 10, Gamma=0.3, 
                        DEBUG=TRUE, standardize=TRUE,  ...){
    if(standardize){
        x <- scale(x, center = TRUE, scale = TRUE)
    }
    commonmRNA <- intersect(colnames(x),V(globalGraph)$name)
    x <- x[, commonmRNA]  
    # calculate the p-value of each gene using Cox PH model
    if(DEBUG) cat('Calculating Cox p-value...')
    colnames(y) <- c("status", "time")
    Survdata.mRNA <- data.frame(x, y,check.names = TRUE)
    geneCoxZP <- matrix(NA,nrow=ncol(x),ncol=2)
    rownames(geneCoxZP) <- colnames(x)
    colnames(geneCoxZP) <- c("ZScore","p-value")
    for(i in seq(from = 1, to = ncol(x))){
        res.coxph <- coxph(as.formula(paste("Surv(time, status)~", 
                        colnames(Survdata.mRNA)[i])), Survdata.mRNA)
        geneCoxZP[i,] <- summary(res.coxph)$coefficients[c(4,5)]
    }
    if(DEBUG) cat('Done\n')
    if(DEBUG) cat('Performing directed random walk...')
    p0 <- rep(0, length = vcount(globalGraph))
    names(p0) <- V(globalGraph)$name
    p0[rownames(geneCoxZP)] <- -log(geneCoxZP[, "p-value"])  
    g.adj <- as_adjacency_matrix(globalGraph, type = "both", sparse = FALSE)
    PT <- rw(W = g.adj, p0 = p0, gamma = Gamma) 
    rownames(PT) <- names(p0)
    PT <- PT[commonmRNA,]
    if(DEBUG) cat('Done\n')
    selGene <- rownames(geneCoxZP)[which(geneCoxZP[, "p-value"] < 0.05)]
    PT <- PT[selGene]
    x <- x[, selGene]
    TPWeight <- 1/order(PT)
    if(DEBUG) cat('Performing cv.glmnet...')
    res <- cv.glmnet(x, y = Surv(time = y[,"time"], event = y[,"status"]),
                family = "cox", nfolds = nfolds, alpha = 1,  # lasso
                penalty.factor = TPWeight, ...)
    if(DEBUG) cat('Done\n') 
    mod <- list(glmnetRes = res, PT = PT)
    class(mod) <- "cvRLassoCox"
    return(mod)
}

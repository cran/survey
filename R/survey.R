library(survey)
library(R.utils)
library(foreign)
library(ff)

ffsvymean <- function (x, design, na.rm = FALSE, deff = FALSE, ...) 
{
    usedff <- FALSE
    
    if (inherits(x, "formula")) {
        usedff <- TRUE
        mf <- model.frame(x, design$variables, na.action = na.pass)
        xx <- lapply(attr(terms(x), "variables")[-1], function(tt) as.ff(model.matrix(eval(bquote(~0 + 
                                                                                                .(tt))), mf)))
        cols <- sapply(lapply(xx, function(x){x[]}), NCOL)
        x <- as.ff(matrix(nrow = NROW(xx[[1]]), ncol = sum(cols)))
        scols <- c(0, cumsum(cols))
        for (i in 1:length(xx)) {
            x[, scols[i] + 1:cols[i]] <- unlist(lapply(xx[i], function(x){x[]}))
        }
        colnames(x) <- do.call("c", lapply(lapply(xx, function(x){x[]}), colnames))
    }
    else {
        if (typeof(x) %in% c("expression", "symbol")) 
            x <- eval(x, design$variables)
        else if (is.data.frame(x) && any(sapply(x, is.factor))) {
            xx <- lapply(x, function(xi) {
                if (is.factor(xi)) 
                    0 + (outer(xi, levels(xi), "=="))
                else xi
            })
            cols <- sapply(xx, NCOL)
            scols <- c(0, cumsum(cols))
            cn <- character(sum(cols))
            for (i in 1:length(xx)) cn[scols[i] + 1:cols[i]] <- paste(names(x)[i], 
                                                                      levels(x[[i]]), sep = "")
            x <- matrix(nrow = NROW(xx[[1]]), ncol = sum(cols))
            for (i in 1:length(xx)) {
                x[, scols[i] + 1:cols[i]] <- xx[[i]]
            }
            colnames(x) <- cn
        }
    }
    if (usedff == TRUE){
        x <- as.matrix(unlist(x[]))
    }else{
        x <- as.matrix(x)
    }
    
    if (na.rm) {
        nas <- rowSums(is.na(x))
        design <- design[nas == 0, ]
        if (length(nas) > length(design$prob)) 
            x <- x[nas == 0, , drop = FALSE]
        else x[nas > 0, ] <- 0
    }
    pweights <- 1/design$prob
    psum <- sum(pweights)
    average <- colSums(x * pweights/psum)
    x <- sweep(x, 2, average)
    v <- svyrecvar(x * pweights/psum, design$cluster, design$strata, 
                   design$fpc, postStrata = design$postStrata)
    attr(average, "var") <- v
    attr(average, "statistic") <- "mean"
    class(average) <- "svystat"
    if (is.character(deff) || deff) {
        nobs <- sum(weights(design) != 0)
        if (deff == "replace") {
            vsrs <- svyvar(x, design, na.rm = na.rm)/(nobs)
        }
        else {
            if (psum < nobs) {
                vsrs <- NA * v
                warning("Sample size greater than population size: are weights correctly scaled?")
            }
            else {
                vsrs <- svyvar(x, design, na.rm = na.rm) * (psum - 
                                                                nobs)/(psum * nobs)
            }
        }
        attr(average, "deff") <- v/vsrs
    }
    return(average)
}

mydata <- 
    read.dta( 
        "http://www.ats.ucla.edu/stat/books/sop/momsag.dta" , 
        convert.factors = FALSE 
    )


mydesign <- 
    svydesign( 
        ids = ~1 , 
        data = mydata , 
        weights = ~weight1 , 
        fpc = ~birth 
    )



origsvymean <- function (x, design, na.rm = FALSE, deff = FALSE, ...) 
{
    if (inherits(x, "formula")) {
        mf <- model.frame(x, design$variables, na.action = na.pass)
        xx <- lapply(attr(terms(x), "variables")[-1], function(tt) model.matrix(eval(bquote(~0 + 
                                                                                                .(tt))), mf))
        cols <- sapply(xx, NCOL)
        x <- matrix(nrow = NROW(xx[[1]]), ncol = sum(cols))
        scols <- c(0, cumsum(cols))
        for (i in 1:length(xx)) {
            x[, scols[i] + 1:cols[i]] <- xx[[i]]
        }
        colnames(x) <- do.call("c", lapply(xx, colnames))
    }
    else {
        if (typeof(x) %in% c("expression", "symbol")) 
            x <- eval(x, design$variables)
        else if (is.data.frame(x) && any(sapply(x, is.factor))) {
            xx <- lapply(x, function(xi) {
                if (is.factor(xi)) 
                    0 + (outer(xi, levels(xi), "=="))
                else xi
            })
            cols <- sapply(xx, NCOL)
            scols <- c(0, cumsum(cols))
            cn <- character(sum(cols))
            for (i in 1:length(xx)) cn[scols[i] + 1:cols[i]] <- paste(names(x)[i], 
                                                                      levels(x[[i]]), sep = "")
            x <- matrix(nrow = NROW(xx[[1]]), ncol = sum(cols))
            for (i in 1:length(xx)) {
                x[, scols[i] + 1:cols[i]] <- xx[[i]]
            }
            colnames(x) <- cn
        }
    }
    x <- as.matrix(x)
    if (na.rm) {
        nas <- rowSums(is.na(x))
        design <- design[nas == 0, ]
        if (length(nas) > length(design$prob)) 
            x <- x[nas == 0, , drop = FALSE]
        else x[nas > 0, ] <- 0
    }
    pweights <- 1/design$prob
    psum <- sum(pweights)
    average <- colSums(x * pweights/psum)
    x <- sweep(x, 2, average)
    v <- svyrecvar(x * pweights/psum, design$cluster, design$strata, 
                   design$fpc, postStrata = design$postStrata)
    attr(average, "var") <- v
    attr(average, "statistic") <- "mean"
    class(average) <- "svystat"
    if (is.character(deff) || deff) {
        nobs <- sum(weights(design) != 0)
        if (deff == "replace") {
            vsrs <- svyvar(x, design, na.rm = na.rm)/(nobs)
        }
        else {
            if (psum < nobs) {
                vsrs <- NA * v
                warning("Sample size greater than population size: are weights correctly scaled?")
            }
            else {
                vsrs <- svyvar(x, design, na.rm = na.rm) * (psum - 
                                                                nobs)/(psum * nobs)
            }
        }
        attr(average, "deff") <- v/vsrs
    }
    return(average)
}

reassignInPackage("svymean.survey.design2", pkgName="survey", ffsvymean); 
rslt <- svymean(~momsag, mydesign)

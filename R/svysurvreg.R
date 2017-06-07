svysurvreg<-function (formula, design, weights=NULL, subset = NULL, ...) 
{
    UseMethod("svysurvreg", design)
}

residuals.svysurvreg<-function(object, type = c("response", "deviance", "dfbeta", 
                                           "dfbetas", "working", "ldcase", "ldresp", "ldshape", "matrix"), 
                               rsigma = TRUE, collapse = FALSE, weighted = TRUE, ...) {
    NextMethod()
}


svysurvreg.survey.design<-
    function (formula, design, weights=NULL, subset=NULL, ...) 
{
    subset <- substitute(subset)
    subset <- eval(subset, model.frame(design), parent.frame())
    if (!is.null(subset)) 
        design <- design[subset, ]
    if (any(weights(design) < 0)) 
        stop("weights must be non-negative")
    data <- model.frame(design)
    g <- match.call()
    g$formula <- eval.parent(g$formula)
    g$design <- NULL
    g$var <- NULL
    if (is.null(g$weights)) 
        g$weights <- quote(.survey.prob.weights)
    else g$weights <- bquote(.survey.prob.weights * .(g$weights))
    g[[1]] <- quote(survreg)
    g$data <- quote(data)
    g$subset <- quote(.survey.prob.weights > 0)
    g$model <- TRUE
    data$.survey.prob.weights <- (1/design$prob)/mean(1/design$prob)
    if (!all(all.vars(formula) %in% names(data))) 
        stop("all variables must be in design= argument")
    g <- with(list(data = data), eval(g))
    g$call <- match.call()
    g$call[[1]] <- as.name(.Generic)
    g$printcall <- sys.call(-1)
    g$printcall[[1]] <- as.name(.Generic)
    class(g) <- c("svysurvreg", class(g))
    g$survey.design <- design
    nas <- g$na.action
    if (length(nas)) 
        design <- design[-nas, ]
    dbeta.subset <- resid(g, "dfbeta", weighted = TRUE)
    if (nrow(design) == NROW(dbeta.subset)) {
        dbeta <- as.matrix(dbeta.subset)
    }
    else {
        dbeta <- matrix(0, ncol = NCOL(dbeta.subset), nrow = nrow(design))
        dbeta[is.finite(design$prob), ] <- dbeta.subset
    }
    g$inv.info <- g$var
    if (inherits(design, "survey.design2")) 
        g$var <- svyrecvar(dbeta, design$cluster, design$strata, 
                           design$fpc, postStrata = design$postStrata)
    else if (inherits(design, "twophase")) 
        g$var <- twophasevar(dbeta, design)
    else if (inherits(design, "twophase2")) 
        g$var <- twophase2var(dbeta, design)
    else if (inherits(design, "pps")) 
        g$var <- ppsvar(dbeta, design)
    else g$var <- svyCprod(dbeta, design$strata, design$cluster[[1]], 
                           design$fpc, design$nPSU, design$certainty, design$postStrata)
    g$ll <- g$loglik
    g$loglik <- NA
    g$degf.resid <- degf(design) - length(coef(g)[!is.na(coef(g))]) + 
        1
    g
}

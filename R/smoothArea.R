svysmoothArea <-  function(formula,
                        domain,
                        design = NULL,
                        adj.mat = NULL,
                        X.domain = NULL,
                        direct.est = NULL,
                        domain.size = NULL,
                        transform = c("identity", "logit", "log"),
                        pc.u = 1,
                        pc.alpha = 0.01,
                        pc.u.phi = 0.5,
                        pc.alpha.phi = 2/3,
                        level = .95,
                        n.sample = 250,
                        var.tol = 1e-10,
                        return.samples = FALSE,...) {
    UseMethod("svysmoothArea",design)
}

svysmoothArea.default <-  function(formula,
                        domain,
                        design = NULL,
                        adj.mat = NULL,
                        X.domain = NULL,
                        direct.est = NULL,
                        domain.size = NULL,
                        transform = c("identity", "logit", "log"),
                        pc.u = 1,
                        pc.alpha = 0.01,
                        pc.u.phi = 0.5,
                        pc.alpha.phi = 2/3,
                        level = .95,
                        n.sample = 250,
                        var.tol = 1e-10,
                        return.samples = FALSE,...) {
   stop("smoothArea is only available for survey.design objects")
}
    
svysmoothArea.survey.design <-  function(formula,
                        domain,
                        design = NULL,
                        adj.mat = NULL,
                        X.domain = NULL,
                        direct.est = NULL,
                        domain.size = NULL,
                        transform = c("identity", "logit", "log"),
                        pc.u = 1,
                        pc.alpha = 0.01,
                        pc.u.phi = 0.5,
                        pc.alpha.phi = 2/3,
                        level = .95,
                        n.sample = 250,
                        var.tol = 1e-10,
                        return.samples = FALSE,...) {
  if (!requireNamespace("SUMMER", quietly = TRUE)) {
    stop(
      "Package \"SUMMER\" must be installed to use this function.",
      call. = FALSE
    )
  }
 rval<-SUMMER::smoothArea(formula = formula,
                            domain = domain,
                            design = design,
                            adj.mat = adj.mat,
                            X.domain = X.domain,
                            direct.est = direct.est,
                            domain.size = domain.size,
                            transform = transform,
                            pc.u = pc.u,
                            pc.alpha = pc.alpha,
                            pc.u.phi = pc.u.phi,
                            pc.alpha.phi = pc.alpha.phi,
                            level = level,
                            n.sample = n.sample,
                            var.tol = var.tol,
                           return.samples = return.samples)
     rval$call<-sys.call(-1)
     rval
}

## print.svysae <- function(x, ...) {
##   x_att <- attributes(x)

##   # print out call
##   cat("Call:\n")
##   print(x$call)
##   cat("\n")

##   # print estimation methods
##   cat("Methods used: ")
##   cat(x_att$method.names, sep = ", ")
##   cat("\n\n")

##   # print transform, if used
##   if (!is.null(x_att$transform)) {
##     cat("Transform: ")
##     cat(x_att$transform)
##     cat("\n\n")
##   }
## }

summary.svysae <- function(object, ...) {
  x<-object
  x_att <- attributes(x)

  # print out call
  cat("Call:\n")
  print(x$call)
  cat("\n")

  # print estimation methods
  cat("Methods used: ")
  cat(x_att$method.names, sep = ", ")
  cat("\n\n")

  # print transform, if used
  if (!is.null(x_att$transform)) {
    cat("Transform: ")
    cat(x_att$transform)
    cat("\n\n")
  }

  # print estimates
  for (i in x_att$method.names) {
    cat(i, "\n")
    print(x[[i]])
    cat("\n")
  }
}


## plot.svysae <- function(x, return_list = F, ...) {
##   combined_est <- do.call(rbind, x[attr(x, "method.names")])

##   # pull the first method (should be "direct.est")
##   ref_method <- attr(x, "method.names")[1]
##   n_methods <- length(attr(x, "method.names"))
##   m <- nrow(x[[ref_method]])
##   sorted_levels <- x[[ref_method]]$domain[order(x[[ref_method]]$mean)]
##   combined_est$domain <- factor(combined_est$domain,
##                                 levels = sorted_levels)
##   combined_est <- combined_est[order(combined_est$domain),]
##   # split across multiple plots for many estimates
##   if (m > 30) {
##     plot_breaks <- c(seq(0, m, by = 30), m)
##     plot_labels <- paste0(
##       "Areas ",
##       (seq_along(plot_breaks)[-1] - 2) * 30 + 1,
##       " to ",
##       pmin((seq_along(plot_breaks)[-1] - 1) * 30, m)
##     )
##     combined_est$plot <-
##       cut(1:m, breaks = plot_breaks, labels = plot_labels)
##   } else {
##     combined_est$plot <- 1
##   }
##   plot_list <- split(combined_est, combined_est$plot)
##   for (y in plot_list) {
##     n_domains <- length(unique(y$domain))
##     plot_dat <- data.frame(
##       domain = y$domain,
##       xval = match(y$domain, unique(y$domain)),
##       est = y$median,
##       lower = y$lower,
##       upper = y$upper,
##       method = as.factor(y$method)
##     )
##     plot_dat$xval <- plot_dat$xval + (as.numeric(plot_dat$method) - 1) *
##       .2 - (.2 * (n_methods - 1) / 2)
##     par(mar = c(9, 5, 4, 3))
##     plot(est ~ xval, data = plot_dat, xaxt = "n", xlab = "", pch = 16,
##          col = (2:(n_methods + 1))[plot_dat$method], cex = .5,
##          ylim = c(min(plot_dat$lower), max(plot_dat$upper)),
##          ylab = "Estimate",
##          cex.axis = .75)
##     segments(plot_dat$xval, plot_dat$lower, plot_dat$xval, plot_dat$upper,
##              col = (2:(n_methods + 1))[plot_dat$method], cex = .5)
##     axis(1, at = 1:n_domains, labels = unique(as.character(y$domain)), las = 2,
##          cex.axis = .65)
##     legend("topleft", inset = .01, legend = levels(plot_dat$method),
##            col=(2:(n_methods + 1)), lty=1, cex=0.8, box.lty=0)

##   }
## }

naa_longer<-function(naa, object,...) UseMethod("naa_longer",naa)
naa_shorter<-function(naa, object,...) UseMethod("naa_shorter",naa)

naa_longer.NULL<-function(naa, object,...) object
naa_shorter.NULL<-function(naa, object,...) object

naa_longer.default<-function(naa, object,...) stop("no default method (not psychic)")
naa_shorter.default<-function(naa, object,...) stop("no default method (not psychic)")

naa_longer.fail<-function(naa, object,...) stop("can't happen (na.fail)")
naa_shorter.fail<-function(naa, object,...) stop("can't happen (na.fail)")

naa_shorter.omit<-function(naa, object,...) object
naa_longer.omit<-function(naa,object,...){ ##from naresid.exclude
    if (length(naa) == 0 || !is.numeric(naa)) 
        stop("invalid argument 'naa'")
    if (is.null(object)) 
        return(object)
    n <- NROW(object)
    keep <- rep.int(NA, n + length(naa))
    keep[-naa] <- 1L:n
    if (is.matrix(object)) {
        object <- object[keep, , drop = FALSE]
        temp <- rownames(object)
        if (length(temp)) {
            temp[naa] <- names(naa)
            rownames(object) <- temp
        }
    }
    else if (is.array(object) && length(d <- dim(object)) > 2L) {
        object <- object[keep, , , drop = FALSE]
        temp <- (dn <- dimnames(object))[[1L]]
        if (!is.null(temp)) {
            temp[naa] <- names(naa)
            dimnames(object)[[1L]] <- temp
        }
    }
    else {
        object <- object[keep]
        temp <- names(object)
        if (length(temp)) {
            temp[naa] <- names(naa)
            names(object) <- temp
        }
    }
    object	
}

naa_longer.exclude<-function(naa,object,...) object
naa_shorter.exclude<-function(naa,object,...) {
    if (length(naa) == 0 || !is.numeric(naa)) 
        stop("invalid argument 'naa'")
    if (is.null(object)) 
        return(object)
    n <- NROW(object)
    keep <- (1:n)[-naa]
    if (is.matrix(object)) {
        object <- object[keep, , drop = FALSE]
        temp <- rownames(object)
     }
    else if (is.array(object) && length(d <- dim(object)) > 2L) {
        object <- object[keep, , , drop = FALSE]
        temp <- (dn <- dimnames(object))[[1L]]
    }
    else {
        object <- object[keep]
        temp <- names(object)
     }
    object		
}





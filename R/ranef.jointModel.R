`ranef.jointModel` <-
function (object, postVar = FALSE, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    out <- as.matrix(object$EB$post.b)
    rownames(out) <- names(object$y$logT)
    if (postVar) {
        n <- nrow(out)
        ncz <- ncol(out)
        vars <- vector("list", n)
        for (i in 1:n) {
            vars[[i]] <- matrix(object$EB$post.vb[i, ], ncz, ncz)
        }
        names(vars) <- rownames(out)
        attr(out, "postVar") <- vars
    }
    out
}


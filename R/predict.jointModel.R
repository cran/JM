predict.jointModel <-
function (object, newdata, se.fit = FALSE, level = 0.95, 
    returnData = FALSE, ...) {
    TermsX <- delete.response(object$termsYx)
    mf <- model.frame(TermsX, data = newdata)
    form <- reformulate(attr(TermsX, "term.labels"))
    X <- model.matrix(form, data = mf)
    out <- c(X %*% object$coefficients$betas)
    names(out) <- row.names(newdata)
    if (se.fit) {
        V <- vcov(object)
        ind <- head(grep("Y.", colnames(V), fixed = TRUE), -1)
        se.fit <- sqrt(diag(X %*% tcrossprod(V[ind, ind], X)))
        alpha <- 1 - level
        low <- out + qnorm(alpha/2) * se.fit
        up <- out + qnorm(1-alpha/2) * se.fit
        names(se.fit) <- names(low) <- names(up) <- row.names(newdata)
        out <- list(pred = out, se.fit = se.fit, low = low, upp = up)
    }
    if (returnData)
        cbind(newdata, if (is.list(out)) do.call(cbind, out) else out)
    else
        out
}


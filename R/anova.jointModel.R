anova.jointModel <-
function (object, object2, test = TRUE, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (!inherits(object2, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (object$method != object2$method)
        stop("You compare joint models with different survival submodels.\n")
    L0 <- logLik(object)
    L1 <- logLik(object2)
    nb0 <- attr(L0, "df")
    nb1 <- attr(L1, "df")
    df <- nb1 - nb0
    if (df < 0)
        stop("'object' should be nested in 'object2'.\n")
    out <- list(nam0 = deparse(substitute(object)), L0 = L0, aic0 = AIC(object), 
        bic0 = AIC(object, k = log(attr(L0, "nobs"))), nam1 = deparse(substitute(object2)), L1 = L1, aic1 = AIC(object2), 
        bic1 = AIC(object2, k = log(attr(L1, "nobs"))), df = df, test = test)
    if (test) {
        LRT <- - 2 * (L0 - L1)
        attributes(LRT) <- NULL
        if (LRT < 0)
            warning("either the two models are not nested or the model represented by 'object2' fell on a local maxima.\n")
        out$LRT <- LRT
        out$p.value <- pchisq(LRT, df, lower.tail = FALSE)
    }
    class(out) <- "aov.jointModel"
    out
}


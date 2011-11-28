cumHaz <-
function (object) {
    method <- object$method
    timeVar <- object$timeVar
    interFact <- object$interFact
    parameterization <- object$parameterization
    derivForm <- object$derivForm
    indFixed <- derivForm$indFixed
    indRandom <- derivForm$indRandom
    LongFormat <- object$LongFormat
    lag <- object$y$lag
    TermsX <- object$termsYx
    TermsZ <- object$termsYz
    TermsX.deriv <- object$termsYx.deriv
    TermsZ.deriv <- object$termsYz.deriv
    formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
    formYz <- object$formYz
    id <- object$id
    W <- object$x$W
    WintF.vl <- object$x$WintF.vl
    WintF.sl <- object$x$WintF.sl
    wk <- gaussKronrod(object$control$GKk)$wk
    sk <- gaussKronrod(object$control$GKk)$sk
    Time0 <- object$times
    last.t <- if (LongFormat) tapply(exp(object$y$logT), object$x$idT, tail, 1) else exp(object$y$logT)
    tt <- mapply(function (x, y) if (!x[1]) c(x[-1], y) else c(x, y),
        split(Time0, id), last.t)
    id <- rep(seq_along(tt), sapply(tt, length))
    Time1 <- unlist(tt, use.names = FALSE)
    P <- Time1 / 2
    P1 <- Time1 / 2
    st <- outer(P, sk) + P1
    id.GK <- rep(seq_along(Time1), each = object$control$GKk)
    id.T <- rep(id, each = object$control$GKk)
    data.id <- object$data.id
    data2 <- data.id[id.T, ]
    data2[[timeVar]] <- pmax(c(t(st)) - lag, 0)
    betas <- object$coefficients$betas
    gammas <- object$coefficients$gammas
    alpha <- object$coefficients$alpha
    Dalpha <- object$coefficients$Dalpha
    b <- ranef(object)
    if (parameterization %in% c("value", "both")) {
        mfX <- model.frame(delete.response(TermsX), data = data2)
        mfZ <- model.frame(TermsZ, data = data2)
        Xs <- model.matrix(formYx, mfX)
        Zs <- model.matrix(formYz, mfZ)
        Ws.intF.vl <- WintF.vl[id.T, , drop = FALSE]
    }
    if (parameterization %in% c("slope", "both")) {
        mfX.deriv <- model.frame(TermsX.deriv, data = data2)
        mfZ.deriv <- model.frame(TermsZ.deriv, data = data2)
        Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
        Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
        Ws.intF.sl <- WintF.sl[id.T, , drop = FALSE]
    }
    if (parameterization %in% c("value", "both"))
        Ys <- as.vector(Xs %*% betas + rowSums(Zs * b[id.T, , drop = FALSE]))
    if (parameterization %in% c("slope", "both"))
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
            rowSums(Zs.deriv * b[id.T, indRandom, drop = FALSE])
    tt <- switch(parameterization,
        "value" = c(Ws.intF.vl %*% alpha) * Ys, 
        "slope" = c(Ws.intF.sl %*% Dalpha) * Ys.deriv,
        "both" = c(Ws.intF.vl %*% alpha) * Ys + 
            c(Ws.intF.sl %*% Dalpha) * Ys.deriv)
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, nrow(data.id))
    if (length(eta.tw) < length(id))
        eta.tw <- eta.tw[id]
    if (method == "weibull-PH-GH") {
        sigma.t <- object$coefficients$sigma.t
        Vi <- exp(log(sigma.t) + (sigma.t - 1) * log(c(t(st))) + tt)
        exp(eta.tw) * P * tapply(wk * Vi, id.GK, sum)    
    } else if (method == "weibull-AFT-GH") {
        sigma.t <- object$coefficients$sigma.t
        Vi <- exp(eta.tw) * P * tapply(wk * tt, id.GK, sum)
        Vi^sigma.t
    } else if (method == "spline-PH-GH") {
        gammas.bs <- object$coefficients$gammas.bs
        W2s <- if (length(kn <- object$control$knots) == 1) {
            splineDesign(unlist(kn, use.names = FALSE), c(t(st)), 
                ord = object$control$ord, outer.ok = TRUE)
        } else {
            strt <- object$y$strata
            strt.s <- rep(strt[id], each = object$control$GKk)
            split.Time <- split(c(t(st)), strt.s)
            w2s <- mapply(function (k, t) splineDesign(k, t, 
                ord = object$control$ord, outer.ok = TRUE), 
                kn, split.Time, SIMPLIFY = FALSE)
            w2s <- mapply(function (w2s, ind) {
                out <- matrix(0, length(st), ncol(w2s))
                out[strt.s == ind, ] <- w2s
                out
            }, w2s, levels(strt), SIMPLIFY = FALSE)
            do.call(cbind, w2s)
        }        
        Vi <- exp(c(W2s %*% gammas.bs) + tt)
        exp(eta.tw) * P * tapply(wk * Vi, id.GK, sum)# <-----------
    } else if (method == "piecewise-PH-GH") {
        xi <- object$coefficients$xi
        qs <- c(0, object$control$knots, max(exp(object$y$logT)) + 1)
        ind.K <- findInterval(c(t(st)), qs, rightmost.closed = TRUE)
        exp(eta.tw) * P * tapply(wk * xi[ind.K] * exp(tt), id.GK, sum)    
    }
}

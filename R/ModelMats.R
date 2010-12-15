ModelMats <-
function (time, ii) {
    if (method %in% c("weibull-AFT-GH", "weibull-PH-GH", "spline-PH-GH", "spline-PH-Laplace")) {
        id.GK <- rep(ii, each = 15)
        wk <- gaussKronrod(15)$wk
        sk <- gaussKronrod(15)$sk
        P <- time / 2
        st <- P * (sk + 1)
        data.id2 <- data.id[id.GK, ]
        data.id2[timeVar] <- pmax(st - lag, 0)
        mfX <- model.frame(TermsX, data = data.id2)
        mfZ <- model.frame(TermsZ, data = data.id2)
        out <- list(st = st, wk = wk, P = P)
        if (parameterization %in% c("value", "both")) {
            out$Xs <- model.matrix(formYx, mfX)
            out$Zs <- model.matrix(formYz, mfZ)
            out$Ws.intF.vl <- WintF.vl[id.GK, , drop = FALSE]
        }
        if (parameterization %in% c("slope", "both")) {
            out$Xs.deriv <- model.matrix(derivForm$fixed, mfX)
            out$Zs.deriv <- model.matrix(derivForm$random, mfZ)
            out$Ws.intF.sl <- WintF.sl[id.GK, , drop = FALSE]
        }      
    }
    if (method == "piecewise-PH-GH") {
        wk <- gaussKronrod(7)$wk
        sk <- gaussKronrod(7)$sk
        nk <- length(sk)
        qs <- c(0, sort(object$control$knots), max(survTimes, object$control$knots) + 1)
        ind <- findInterval(time, qs, rightmost.closed = TRUE)
        Tiq <- outer(time, qs, pmin)
        Lo <- Tiq[, 1:Q]
        Up <- Tiq[, 2:(Q+1)]
        T <- Up - Lo
        P <- T / 2
        P[P < sqrt(.Machine$double.eps)] <- as.numeric(NA)
        P1 <- (Up + Lo) / 2
        st <- rep(P, each = nk) * rep(sk, Q) + rep(P1, each = nk)
        data.id2 <- data.id[rep(ii, each = nk*Q), ]
        data.id2[timeVar] <- pmax(st - lag, 0)
        id.GK <- rep(ii, sum(!is.na(st)))
        mfX <- model.frame(TermsX, data = data.id2)
        mfZ <- model.frame(TermsZ, data = data.id2)
        out <- list(st = st, wk = wk, P = P, ind = ind)
        if (parameterization %in% c("value", "both")) {
            out$Xs <- model.matrix(formYx, mfX)
            out$Zs <- model.matrix(object$formYz, mfZ)
            out$Ws.intF.vl <- WintF.vl[id.GK, , drop = FALSE]
        }
        if (parameterization %in% c("slope", "both")) {
            out$Xs.deriv <- model.matrix(derivForm$fixed, mfX)
            out$Zs.deriv <- model.matrix(derivForm$random, mfZ)
            out$Ws.intF.sl <- WintF.sl[id.GK, , drop = FALSE]
        }
    }
    out
}


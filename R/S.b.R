S.b <-
function (t, b, i) {
    if (t == 0)
        return(1)
    log.survival <- if (method == "weibull-PH-GH") {
        id.i <- rep(i, each = 15)
        wk <- gaussKronrod(15)$wk
        sk <- gaussKronrod(15)$sk
        P <- t/2
        st <- P * (sk + 1)
        log.st <- log(st)
        data.id2 <- data.id[id.i, ]
        data.id2[timeVar] <- pmax(st - object$y$lag, 0)
        mf <- model.frame(TermsY, data = data.id2)
        if (parameterization %in% c("value", "both")) {
            Xs <- model.matrix(formYx, mf)
            Zs <- model.matrix(object$formYz, mf)
            Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * b[id.i, , drop = FALSE]))
        }
        if (parameterization %in% c("slope", "both")) {
            Xs.deriv <- model.matrix(derivForm$fixed, mf)
            Zs.deriv <- model.matrix(derivForm$random, mf)
            Ys.deriv <- as.vector(Xs.deriv %*% betas.new[indFixed] + 
                rowSums(Zs.deriv * b[id.i, indRandom, drop = FALSE]))
        }
        tt <- switch(parameterization, "value" = alpha.new * Ys, 
            "slope" = Dalpha.new * Ys.deriv, "both" = alpha.new * Ys + Dalpha.new * Ys.deriv)
        eta.tw <- as.vector(W %*% gammas.new)
        Vi <- exp(log(sigma.t.new) + (sigma.t.new - 1) * log.st + tt)
        - exp(eta.tw[i]) * P * sum(wk * Vi)
    } else if (method == "weibull-AFT-GH") {
        id.i <- rep(i, each = 15)
        wk <- gaussKronrod(15)$wk
        sk <- gaussKronrod(15)$sk
        P <- t/2
        st <- P * (sk + 1)
        data.id2 <- data.id[id.i, ]
        data.id2[timeVar] <- pmax(st - object$y$lag, 0)
        mf <- model.frame(TermsY, data = data.id2)
        if (parameterization %in% c("value", "both")) {
            Xs <- model.matrix(formYx, mf)
            Zs <- model.matrix(object$formYz, mf)
            Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * b[id.i, , drop = FALSE]))
        }
        if (parameterization %in% c("slope", "both")) {
            Xs.deriv <- model.matrix(derivForm$fixed, mf)
            Zs.deriv <- model.matrix(derivForm$random, mf)
            Ys.deriv <- as.vector(Xs.deriv %*% betas.new[indFixed] + 
                rowSums(Zs.deriv * b[id.i, indRandom, drop = FALSE]))
        }
        tt <- switch(parameterization, "value" = alpha.new * Ys, 
            "slope" = Dalpha.new * Ys.deriv, "both" = alpha.new * Ys + Dalpha.new * Ys.deriv)
        eta.tw <- as.vector(W %*% gammas.new)
        Vi <- exp(eta.tw[i]) * P * sum(wk * exp(tt))
        - Vi^sigma.t.new
    } else if (method == "spline-PH-GH") {
        id.i <- rep(i, each = 15)
        wk <- gaussKronrod(15)$wk
        sk <- gaussKronrod(15)$sk
        P <- t/2
        st <- P * (sk + 1)
        log.st <- log(st)
        data.id2 <- data.id[id.i, ]
        data.id2[timeVar] <- pmax(st - object$y$lag, 0)
        mf <- model.frame(TermsY, data = data.id2)
        if (parameterization %in% c("value", "both")) {
            Xs <- model.matrix(formYx, mf)
            Zs <- model.matrix(object$formYz, mf)
            Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * b[id.i, , drop = FALSE]))
        }
        if (parameterization %in% c("slope", "both")) {
            Xs.deriv <- model.matrix(derivForm$fixed, mf)
            Zs.deriv <- model.matrix(derivForm$random, mf)
            Ys.deriv <- as.vector(Xs.deriv %*% betas.new[indFixed] + 
                rowSums(Zs.deriv * b[id.i, indRandom, drop = FALSE]))
        }
        tt <- switch(parameterization, "value" = alpha.new * Ys, 
            "slope" = Dalpha.new * Ys.deriv, "both" = alpha.new * Ys + Dalpha.new * Ys.deriv)
        eta.tw <- if (!is.null(W)) as.vector(W[i, , drop = FALSE] %*% gammas.new) else 0
        W2s <- if (length(kn <- object$control$knots) == 1) {
            splineDesign(unlist(kn, use.names = FALSE), st, ord = object$control$ord, outer.ok = TRUE)
        } else {
            strt.i <- strt[i]
            w2s <- lapply(kn, function (kn) splineDesign(kn, st, ord = object$control$ord, outer.ok = TRUE))
            ll <- match(strt.i, names(w2s))
            w2s[-ll] <- lapply(w2s[-ll], function (m) {m[, ] <- 0; m})
            do.call(cbind, w2s)
        }
        Vi <- exp(c(W2s %*% gammas.bs.new) + tt)
        - exp(eta.tw) * P * sum(wk * Vi)
    } else if (method == "piecewise-PH-GH") {
        wk <- gaussKronrod(7)$wk
        sk <- gaussKronrod(7)$sk
        nk <- length(sk)
        qs <- c(0, sort(object$control$knots), max(survTimes, object$control$knots) + 1)
        ind <- findInterval(t, qs, rightmost.closed = TRUE)
        Tiq <- outer(t, qs, pmin)
        Lo <- Tiq[, 1:Q]
        Up <- Tiq[, 2:(Q+1)]
        T <- Up - Lo
        P <- T / 2
        P[P < sqrt(.Machine$double.eps)] <- as.numeric(NA)
        P1 <- (Up + Lo) / 2
        st <- rep(P, each = nk) * rep(sk, Q) + rep(P1, each = nk)
        data.id2 <- data.id[rep(i, each = nk*Q), ]
        data.id2[timeVar] <- pmax(st - object$y$lag, 0)
        mf <- model.frame(TermsY, data = data.id2)
        if (parameterization %in% c("value", "both")) {
            Xs <- model.matrix(formYx, mf)
            Zs <- model.matrix(object$formYz, mf)
            Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * rep(b[i, ], each = nrow(Zs))))
        }
        if (parameterization %in% c("slope", "both")) {
            Xs.deriv <- model.matrix(derivForm$fixed, mf)
            Zs.deriv <- model.matrix(derivForm$random, mf)
            Ys.deriv <- as.vector(Xs.deriv %*% betas.new[indFixed] + 
                rowSums(Zs.deriv * rep(b[i, indRandom], each = nrow(Zs))))
        }
        tt <- switch(parameterization, "value" = alpha.new * Ys, 
            "slope" = Dalpha.new * Ys.deriv, "both" = alpha.new * Ys + Dalpha.new * Ys.deriv)
        P <- P[!is.na(P)]
        ind.K <- rep(seq_len(ind), each = nk)
        wk <- rep(wk, ind)
        wkP <- wk * rep(P, each = nk)
        eta.tw <- if (!is.null(W)) as.vector(W[i, , drop = FALSE] %*% gammas.new) else 0 
        - exp(eta.tw) * sum(xi.new[ind.K] * wkP * exp(tt))
    }
    exp(log.survival)
}


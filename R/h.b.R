h.b <-
function (t, b, i) {
    log.hazard <- if (method == "weibull-PH-GH") {
        data.id2 <- data.id[i, ]
        data.id2[timeVar] <- pmax(t - object$y$lag, 0)
        mf <- model.frame(TermsY, data = data.id2)
        Xs <- model.matrix(object$formYx, mf)
        Zs <- model.matrix(object$formYz, mf)
        eta.tw <- as.vector(W %*% gammas.new)
        Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * b[rep(i, nrow(Zs)), , drop = FALSE]))
        log(sigma.t.new) + (sigma.t.new - 1) * log(t) + eta.tw[i] + alpha.new * Ys
    } else if (method == "weibull-AFT-GH") {
        id.i <- rep(i, each = 15)
        wk <- gaussKronrod(15)$wk
        sk <- gaussKronrod(15)$sk
        P <- t/2
        st <- P * (sk + 1)
        data.id2 <- data.id[id.i, ]
        data.id2[timeVar] <- pmax(st - object$y$lag, 0)
        mf <- model.frame(TermsY, data = data.id2)
        Xs <- model.matrix(object$formYx, mf)
        Zs <- model.matrix(object$formYz, mf)
        Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * b[id.i, , drop = FALSE]))
        data.id3 <- data.id[i, ]
        data.id3[timeVar] <- pmax(t - object$y$lag, 0)
        mf <- model.frame(TermsY, data = data.id3)
        X <- model.matrix(object$formYx, mf)
        Z <- model.matrix(object$formYz, mf)
        Y <- as.vector(X %*% betas.new + rowSums(Z * b[rep(i, nrow(Z)), , drop = FALSE]))
        eta.tw <- as.vector(W %*% gammas.new)
        Vi <- exp(eta.tw[i]) * P * sum(wk * exp(alpha.new * Ys))
        log(sigma.t.new) + (sigma.t.new - 1) * log(Vi) + eta.tw[i] + alpha.new * Y
    } else if (method == "piecewise-PH-GH") {
        qs <- c(0, sort(object$control$knots), max(survTimes) + 1)
        ind <- findInterval(t, qs, rightmost.closed = TRUE)
        data.id2 <- data.id[i, ]
        data.id2[timeVar] <- pmax(t - object$y$lag, 0)
        mf <- model.frame(TermsY, data = data.id2)
        Xs <- model.matrix(object$formYx, mf)
        Zs <- model.matrix(object$formYz, mf)
        eta.tw <- if (!is.null(W)) as.vector(W[i, , drop = FALSE] %*% gammas.new) else 0 
        Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * b[rep(i, nrow(Zs)), ]))
        log(xi[ind]) + eta.tw + alpha.new * Ys
    }
    exp(log.hazard)
}


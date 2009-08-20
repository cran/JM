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
        data.id2[timeVar] <- st
        mf <- model.frame(TermsY, data = data.id2)
        Xs <- model.matrix(object$formYx, mf)
        Zs <- model.matrix(object$formYz, mf)
        eta.tw <- as.vector(W %*% gammas.new)
        Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * b[id.i, , drop = FALSE]))
        Vi <- exp(log(sigma.t.new) + (sigma.t.new - 1) * log.st + alpha.new * Ys)
        - exp(eta.tw[i]) * P * sum(wk * Vi)
    } else if (method == "weibull-AFT-GH") {
        id.i <- rep(i, each = 15)
        wk <- gaussKronrod(15)$wk
        sk <- gaussKronrod(15)$sk
        P <- t/2
        st <- P * (sk + 1)
        data.id2 <- data.id[id.i, ]
        data.id2[timeVar] <- st
        mf <- model.frame(TermsY, data = data.id2)
        Xs <- model.matrix(object$formYx, mf)
        Zs <- model.matrix(object$formYz, mf)
        eta.tw <- as.vector(W %*% gammas.new)
        Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * b[id.i, , drop = FALSE]))
        Vi <- exp(eta.tw[i]) * P * sum(wk * exp(alpha.new * Ys))
        - Vi^sigma.t.new
    } else if (method == "piecewise-PH-GH") {
        wk <- gaussKronrod(7)$wk
        sk <- gaussKronrod(7)$sk
        nk <- length(sk)
        qs <- c(0, sort(object$control$knots), max(survTimes) + 1)
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
        data.id2[timeVar] <- st       
        mf <- model.frame(TermsY, data = data.id2)
        Xs <- model.matrix(object$formYx, mf)
        Zs <- model.matrix(object$formYz, mf)
        P <- P[!is.na(P)]
        ind.K <- rep(seq_len(ind), each = nk)
        wk <- rep(wk, ind)
        wkP <- wk * rep(P, each = nk)
        eta.tw <- if (!is.null(W)) as.vector(W[i, , drop = FALSE] %*% gammas.new) else 0 
        Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * rep(b[i, ], each = nrow(Zs))))
        - exp(eta.tw) * sum(xi.new[ind.K] * wkP * exp(alpha.new * Ys))
    }
    exp(log.survival)
}


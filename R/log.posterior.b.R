log.posterior.b <-
function (b, time, method, ii) {
    id.i <- id %in% ii
    X.i <- X[id.i, , drop = FALSE]
    Z.i <- Z[id.i, , drop = FALSE]
    mu.y <- as.vector(X.i %*% betas.new) + rowSums(Z.i * rep(b, each = nrow(Z.i)))
    logNorm <- dnorm(y[id.i], mu.y, sigma.new, TRUE)
    log.p.yb <- sum(logNorm)
    log.p.b <- dmvnorm(b, rep(0, ncz), D.new, TRUE)
    log.survival <- if (time[ii] == 0) 0 else {
        if (method == "weibull-PH-GH") {
            id.GK <- rep(ii, each = 15)
            wk <- gaussKronrod(15)$wk
            sk <- gaussKronrod(15)$sk
            P <- time[ii]/2
            st <- P * (sk + 1)
            log.st <- log(st)
            data.id2 <- data.id[id.GK, ]
            data.id2[timeVar] <- st
            mf <- model.frame(TermsY, data = data.id2)
            Xs <- model.matrix(object$formYx, mf)
            Zs <- model.matrix(object$formYz, mf)
            eta.tw <- as.vector(W[ii, , drop = FALSE] %*% gammas.new)
            Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * rep(b, each = nrow(Zs))))
            Vi <- exp(log(sigma.t.new) + (sigma.t.new - 1) * log.st + alpha.new * Ys)
            - exp(eta.tw) * P * sum(wk * Vi)
        } else if (method == "weibull-AFT-GH") {
            id.GK <- rep(ii, each = 15)
            wk <- gaussKronrod(15)$wk
            sk <- gaussKronrod(15)$sk
            P <- time[ii]/2
            st <- P * (sk + 1)
            log.st <- log(st)
            data.id2 <- data.id[id.GK, ]
            data.id2[timeVar] <- st
            mf <- model.frame(TermsY, data = data.id2)
            Xs <- model.matrix(object$formYx, mf)
            Zs <- model.matrix(object$formYz, mf)
            eta.tw <- as.vector(W[ii, , drop = FALSE] %*% gammas.new)
            Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * rep(b, each = nrow(Zs))))
            Vi <- exp(eta.tw) * P * sum(wk * exp(alpha.new * Ys))
            - Vi^sigma.t.new
        } else if (method == "piecewise-PH-GH") {
            wk <- gaussKronrod(7)$wk
            sk <- gaussKronrod(7)$sk
            nk <- length(sk)
            qs <- c(0, sort(object$control$knots), max(survTimes) + 1)
            ind <- findInterval(time[ii], qs, rightmost.closed = TRUE)
            Tiq <- outer(time[ii], qs, pmin)
            Lo <- Tiq[, 1:Q]
            Up <- Tiq[, 2:(Q+1)]
            T <- Up - Lo
            P <- T / 2
            P[P < sqrt(.Machine$double.eps)] <- as.numeric(NA)
            P1 <- (Up + Lo) / 2
            st <- rep(P, each = nk) * rep(sk, Q) + rep(P1, each = nk)
            data.id2 <- data.id[rep(ii, each = nk*Q), ]
            data.id2[timeVar] <- st       
            mf <- model.frame(TermsY, data = data.id2)
            Xs <- model.matrix(object$formYx, mf)
            Zs <- model.matrix(object$formYz, mf)
            P <- P[!is.na(P)]
            ind.K <- rep(seq_len(ind), each = nk)
            wk <- rep(wk, ind)
            wkP <- wk * rep(P, each = nk)
            eta.tw <- if (!is.null(W)) as.vector(W[ii, , drop = FALSE] %*% gammas.new) else 0 
            Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * rep(b, each = nrow(Zs))))
            - exp(eta.tw) * sum(xi.new[ind.K] * wkP * exp(alpha.new * Ys))
        }
    }
    log.p.yb + log.survival + log.p.b
}


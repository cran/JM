Score.chGH <-
function (thetas) {
    betas <- thetas[1:ncx]
    sigma <- exp(thetas[ncx + 1])
    gammas <- thetas[seq(ncx + 2, ncx + 1 + ncww)]
    gammas[1:nk] <- cumsum(c(gammas[1], exp(gammas[2:nk])))
    alpha <- thetas[ncx + ncww + 2]
    D <- thetas[seq(ncx + ncww + 3, length(thetas))]
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- as.vector(WW %*% gammas)
    Y <- eta.yxT + Ztime.b
    eta.t <- eta.tw + alpha * Y
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    sc <- as.vector(S %*% diff(gammas[1:nk]))
    ew <- - exp(eta.t)
    log.p.tb <- d * (log(sc) + eta.t + ew - logT) + (1 - d) * ew
    log.p.b <- if (ncz == 1) {
        dnorm(b, sd = sqrt(D), log = TRUE)
    } else {
        if (diag.D) {
            rowSums(dnorm(b, sd = rep(sqrt(D), each = k), log = TRUE))
        } else {
            dmvnorm(b, rep(0, ncz), D, TRUE)
        }
    }
    p.ytb <- exp((log.p.yb + log.p.tb) + rep(log.p.b, each = n)); dimnames(p.ytb) <- NULL
    p.yt <- c(p.ytb %*% wGH)
    p.byt <- p.ytb / p.yt
    post.b <- p.byt %*% (b * wGH)
    post.vb <- if (ncz == 1) {
            c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
    } else {
            (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
    }
    Zb <- if (ncz == 1) post.b[id] else rowSums(Z * post.b[id, ], na.rm = TRUE)
    mu <- y - eta.yx
    tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)
    ew1 <- - exp(eta.t)
    ew2 <- ew1 * Ztime.b
    out <- (d + c((ew1 * p.byt) %*% wGH)) * WW
    out[, 1:nk] <- out[, 1:nk] + d * SS / sc
    out <- colSums(out, na.rm = TRUE)
    out[1:nk] <- c(out[1:nk] %*% jacobian(thetas[seq(ncx + 2, ncx + 1 + nk)]))
    sc.alpha <- sum(((d * Y + eta.yxT * ew1 + ew2) * p.byt) %*% wGH, na.rm = TRUE)
    score.t <- - c(out, sc.alpha)
    score.y <- - c(crossprod(X, mu - Zb) / sigma^2 + colSums((d + c((ew1 * p.byt) %*% wGH)) * (alpha * Xtime), na.rm = TRUE),
        sigma * (- N / sigma + drop(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / sigma^3))
    score.b <- if (diag.D) {
        svD <- 1 / D
        svD2 <- svD^2
        cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
        dim(cS.postVB) <- c(ncz, ncz)
        D * 0.5 * (n * svD - diag(cS.postVB) * svD2 - colSums(as.matrix(post.b^2), na.rm = TRUE) * svD2)
    } else {
        svD <- solve(D)
        dD <- deriv.D(D)
        ndD <- length(dD)
        D1 <- sapply(dD, function (x) sum(svD * x))
        D2 <- t(sapply(dD, function (x) c(svD %*% x %*% svD)))
        cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
        out <- numeric(ndD)
        for (i in seq_along(dD)) {
            D.mat <- D2[i, ]
            dim(D.mat) <- c(ncz, ncz)
            out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) + sum((post.b %*% D.mat) * post.b, na.rm = TRUE)   
        }
        J <- jacobian2(attr(D, "L"), ncz)
        drop(0.5 * (n * D1 - out) %*% J)
    }
    c(score.y, score.t, score.b)
}


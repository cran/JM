`Score.phGH` <-
function (thetas) {
    betas <- thetas[1:ncx]
    sigma <- exp(thetas[ncx + 1])
    gammas <- thetas[seq(ncx + 2, ncx + 1 + ncww)]
    alpha <- thetas[ncx + ncww + 2]
    D <- thetas[seq(ncx + ncww + 3, length(thetas))]
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.yxT2 <- as.vector(Xtime2 %*% betas)
    Y <- eta.yxT + Ztime.b
    Y2 <- eta.yxT2 + Ztime2.b
    eta.t <- if (is.null(WW)) alpha * Y else as.vector(WW %*% gammas) + alpha * Y
    eta.t2 <- if (is.null(WW)) alpha * Y2 else as.vector(WW %*% gammas)[indT] + alpha * Y2    
    lambda0T <- d * lambda0[ind.T0]
    lambda0T[is.na(lambda0T)] <- 0
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    ew <- exp(eta.t2)
    log.p.tb <- d * (log(lambda0T) + eta.t)
    log.p.tb[ind0, ] <- 0
    log.p.tb[ind.lenN0, ] <- log.p.tb[ind.lenN0, ] - rowsum(lambda0[ind.lambda] * ew, indT)
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
    WX <- if (is.null(WW)) alpha * Xtime else cbind(WW, alpha * Xtime)
    WX2 <- if (is.null(WW)) alpha * Xtime2 else cbind(WW[indT, , drop = FALSE], alpha * Xtime2)
    nn <- ncol(WX)
    res <- numeric(nn)
    for (i in 1:nn) {
        p <- rowsum(lambda0[ind.lambda] * WX2[, i] * ew, indT)
        res[i] <- sum((p * p.byt[ind.lenN0, ]) %*% wGH)
    }
    sc <- crossprod(WX, d) - res
    sc.gammas <- if (is.null(WW)) NULL else sc[1:ncww]
    sc.betas <- if (is.null(WW)) sc else sc[seq(ncww + 1, nn)]
    p1 <- sum(d * c((Y * p.byt) %*% wGH))
    p2 <- rowsum(lambda0[ind.lambda] * Y2 * ew, indT)
    sc.alpha <- p1 - sum((p2 * p.byt[ind.lenN0, ]) %*% wGH, na.rm = TRUE)
    score.t <- - c(sc.gammas, sc.alpha)
    score.y <- - c(crossprod(X, mu - Zb) / sigma^2 + sc.betas, 
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


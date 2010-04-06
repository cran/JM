Score.weibullAFTGH <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    sigma.t <- if (is.null(scaleWB)) exp(thetas$log.sigma.t) else scaleWB
    D <- thetas$D
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- as.vector(WW %*% gammas)
    exp.eta.tw <- exp(eta.tw)
    Y <- eta.yxT + Ztime.b
    Ys <- as.vector(Xs %*% betas) + Zsb
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Ys
    wk.exp.eta.s <- wk * exp(eta.s)
    exp.eta.tw <- exp(eta.tw)
    exp.eta.tw.P <- exp(eta.tw) * P
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    Vi <- exp.eta.tw * P * rowsum(wk.exp.eta.s, id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
    log.survival <- - Vi^sigma.t
    log.p.tb <- d * log.hazard + log.survival
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
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    Vii <- d * (sigma.t - 1) / Vi - sigma.t * Vi^(sigma.t - 1)
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
        ki <- Vii * exp.eta.tw.P * rowsum(wk.exp.eta.s * alpha * Xs[, i], id.GK, reorder = FALSE)
        kii <- c((p.byt * ki) %*% wGH)
        sc2[i] <- - sum(d * alpha * Xtime[, i] + kii, na.rm = TRUE)
    }    
    score.y <- c(sc1 + sc2, - sigma * (- N / sigma + drop(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / sigma^3))
    Vi <- exp.eta.tw * P * rowsum(wk.exp.eta.s, id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    Vii <- d * (sigma.t - 1) / Vi - sigma.t * Vi^(sigma.t - 1)
    scgammas <- - colSums(WW * (d + c((p.byt * Vii * Vi) %*% wGH)), na.rm = TRUE)
    scalpha <- - sum((p.byt * (d * Y + Vii * exp.eta.tw * P * rowsum(wk.exp.eta.s * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    score.t <- if (is.null(scaleWB)) {
        scsigmat <- - sigma.t * sum((p.byt * (d / sigma.t + (d - Vi^sigma.t) * log(Vi))) %*% wGH, na.rm = TRUE)
        c(scgammas, scalpha, scsigmat)
    } else {
        c(scgammas, scalpha)
    }
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


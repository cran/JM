Score.splineGH <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    gammas.bs <- thetas$gammas.bs
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    D <- thetas$D
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.tw1 <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, n)
    eta.tw2 <- as.vector(W2 %*% gammas.bs)
    if (parameterization %in% c("value", "both")) {
        Y <- as.vector(Xtime %*% betas) + Ztime.b
        Ys <- as.vector(Xs %*% betas) + Zsb
        eta.t <- eta.tw2 + eta.tw1 + alpha * Y
        eta.s <- alpha * Ys
    }
    if (parameterization %in% c("slope", "both")) {
        Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + Ztime.b.deriv
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
        eta.t <- if (parameterization == "both") eta.t + Dalpha * Y.deriv else eta.tw2 + eta.tw1 + Dalpha * Y.deriv
        eta.s <- if (parameterization == "both") eta.s + Dalpha * Ys.deriv else Dalpha * Ys.deriv
    }    
    eta.ws <- as.vector(W2s %*% gammas.bs)
    exp.eta.tw.P <- exp(eta.tw1) * P
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    log.hazard <- eta.t
    Int <- wk * exp(eta.ws + eta.s)
    log.survival <- - exp.eta.tw.P * rowsum(Int, id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
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
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
        ki <- exp.eta.tw.P * switch(parameterization,
            "value" = rowsum(Int * alpha * Xs[, i], id.GK, reorder = FALSE),
            "slope" = {ii <- match(i, indFixed); 
                if (is.na(ii)) 0 else rowsum(Int * Dalpha * Xs.deriv[, ii], id.GK, reorder = FALSE)},
            "both" = {ii <- match(i, indFixed);
                rowsum(Int * (alpha * Xs[, i] + Dalpha * if (is.na(ii)) 0 else Xs.deriv[, ii]), id.GK, reorder = FALSE)}
        )
        kii <- c((p.byt * ki) %*% wGH)
        sc2[i] <- switch(parameterization,
            "value" = - sum(d * alpha * Xtime[, i] - kii, na.rm = TRUE),
            "slope" = {ii <- match(i, indFixed); 
                if (is.na(ii)) 0 else - sum(d * Dalpha * Xtime.deriv[, ii] - kii, na.rm = TRUE)},
            "both" = {ii <- match(i, indFixed);
                - sum(d * (alpha * Xtime[, i] + Dalpha * if (is.na(ii)) 0 else Xtime.deriv[, ii]) - kii, na.rm = TRUE)}
        )        
        
    }    
    score.y <- c(sc1 + sc2, - sigma * (- N / sigma + drop(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / sigma^3))
    scgammas1 <- if (!is.null(W1)) - colSums(W1 * (d + c((p.byt * log.survival) %*% wGH)), na.rm = TRUE) else NULL    
    scgammas2 <- numeric(nk)
    for (i in 1:nk) {
        kk <- exp.eta.tw.P * rowsum(Int * W2s[, i], id.GK, reorder = FALSE)
        scgammas2[i] <- - sum(W2[, i] * d - c((p.byt * kk) %*% wGH))
    }
    scalpha <- if (parameterization %in% c("value", "both")) {
        - sum((p.byt * (d * Y - exp.eta.tw.P * rowsum(Int * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)    
    } else NULL
    scalpha.D <- if (parameterization %in% c("slope", "both")) {
        - sum((p.byt * (d * Y.deriv - exp.eta.tw.P * rowsum(Int * Ys.deriv, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    } else NULL
    score.t <- c(scgammas1, scalpha, scalpha.D, scgammas2)
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


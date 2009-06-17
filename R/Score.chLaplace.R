Score.chLaplace <-
function (thetas, b) {
    betas <- thetas[1:ncx]
    sigma <- exp(thetas[ncx + 1])
    gammas <- thetas[seq(ncx + 2, ncx + 1 + ncww)]
    gammas[1:nk] <- cumsum(c(gammas[1], exp(gammas[2:nk])))
    alpha <- thetas[ncx + ncww + 2]
    D <- thetas[seq(ncx + ncww + 3, length(thetas))]
    D <- if (diag.D) exp(D) else chol.transf(D)
    environment(update.bCH) <- environment(fn.b) <- environment(gr.b) <- environment()
    environment(logsurvCH) <- environment(ScsurvCH) <- environment()
    new.b <- update.bCH(b, matrix(0, n, ncz * ncz), betas, sigma, c(gammas, alpha), D)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- as.vector(WW %*% gammas)
    sc <- as.vector(S %*% diff(gammas[1:nk]))
    for (i in 1:n) {
        eta.twi <- eta.tw[i]
        Ztime.i <- Ztime[i, ]
        bb <- attr(new.b, "b")[i, ]
        hes.b <- attr(new.b, "hes.b")[i, ]
        var.b <- attr(new.b, "vb")[i, ]            
        dim(var.b) <- dim(hes.b) <- c(ncz, ncz)
        Ztime.b <- sum(Ztime.i * bb)
        eta.ti <- eta.twi + alpha * (eta.yxT[i] + Ztime.b)
        exp.eta.ti <- exp(eta.ti)
        trc1 <- - (alpha^3 * exp.eta.ti) * outer.Ztime[[i]]
        trc2 <- colSums(Ztime[i, ] * var.b)
        K <- var.b %*% trc1
        tr.var.b.trc1 <- - 0.5 * sum(diag(K))
        trc.y1[i, ] <- tr.var.b.trc1 * trc2
        trc.y2[i, ] <- - 0.5 * sum(- K * t(K)) * c(trc2 %o% trc2)
        L <- alpha * c(trc2 %o% trc2)
        M <- c(trc2 %o% colSums(Ztime.i * (K %*% var.b)))
        trc.y3[i, ] <- tr.var.b.trc1 * (L + M)
        ZtSZ <- c(crossprod(Ztime.i, solve(hes.b, Ztime.i)))
        P <- alpha * exp.eta.ti * ZtSZ - 1 / alpha
        trc.t1[i] <- tr.var.b.trc1 * P
        Q <- exp.eta.ti * ZtSZ * (alpha * Ztime.b + 1) - (alpha * Ztime.b + 2) / alpha^2
        trc.t2[i] <- tr.var.b.trc1 * Q
    }
    b <- attr(new.b, "b")
    hes.b <- attr(new.b, "hes.b")
    b.hat <- b + trc.y1
    vb.hat <- attr(new.b, "vb") + trc.y2 + trc.y3
    Zb <- rowSums(Z * b.hat[id, ])
    btZtZb <- drop(crossprod(Zb))
    tr.tZZvarb <- sum(ZtZ * vb.hat)
    Ztime.b <- rowSums(Ztime * b)
    mu <- y - eta.yx
    Y <- eta.yxT + Ztime.b
    eta.t <- eta.tw + Y * alpha
    ew1 <- - exp(eta.t)
    ew2 <- ew1 - trc.t1
    ew3 <- ew1 * Ztime.b - trc.t2
    out <- (d + ew2) * WW
    out[, 1:nk] <- out[, 1:nk] + d * SS / sc
    out <- colSums(out, na.rm = TRUE)
    out[1:nk] <- c(out[1:nk] %*% jacobian(thetas[seq(ncx + 2, ncx + 1 + nk)]))
    score.t <- - c(out, sum(d * Y + eta.yxT * ew1 + ew3, na.rm = TRUE))
    score.y <- - c(crossprod(X, mu - Zb) / sigma^2 + colSums((d + ew2) * (alpha * Xtime)), 
        sigma * (- N / sigma + drop(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / sigma^3))
    score.b <- if (diag.D) {
        svD <- 1 / D
        svD2 <- svD^2
        cS.postVB <- colSums(as.matrix(vb.hat), na.rm = TRUE)
        dim(cS.postVB) <- c(ncz, ncz)
        D * 0.5 * (n * svD - diag(cS.postVB) * svD2 - colSums(as.matrix(b.hat^2), na.rm = TRUE) * svD2)
    } else {
        svD <- solve(D)
        dD <- deriv.D(D)
        ndD <- length(dD)
        D1 <- sapply(dD, function (x) sum(svD * x))
        D2 <- t(sapply(dD, function (x) c(svD %*% x %*% svD)))
        cS.postVB <- colSums(as.matrix(vb.hat), na.rm = TRUE)
        out <- numeric(ndD)
        for (i in seq_along(dD)) {
            D.mat <- D2[i, ]
            dim(D.mat) <- c(ncz, ncz)
            out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) + sum((b.hat %*% D.mat) * b.hat, na.rm = TRUE)   
        }
        J <- jacobian2(attr(D, "L"), ncz)
        drop(0.5 * (n * D1 - out) %*% J)
    }
    c(score.y, score.t, score.b)
}


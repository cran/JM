gr.survSplinePH <-
function (thetas) {
    gammas <- thetas[1:ncww]
    gammas.bs <- gammas[1:nk]
    gammas <- gammas[-(1:nk)]
    alpha <- thetas[ncww + 1]
    eta.tw1 <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, n)
    eta.tw2 <- as.vector(W2 %*% gammas.bs)
    eta.t <- eta.tw2 + eta.tw1 + alpha * Y
    eta.s <- alpha * Ys
    eta.ws <- as.vector(W2s %*% gammas.bs)
    exp.eta.tw.P <- exp(eta.tw1) * P
    Int <- wk * exp(eta.ws + eta.s)
    scgammas1 <- if (!is.null(W1)) {
        ki <- exp.eta.tw.P * rowsum(Int, id.GK, reorder = FALSE)
        kii <- c((p.byt * ki) %*% wGH)
        - colSums(W1 * (d - kii), na.rm = TRUE)
    } else 
        NULL
    scgammas2 <- numeric(nk)
    for (i in 1:nk) {
        kk <- exp.eta.tw.P * rowsum(Int * W2s[, i], id.GK, reorder = FALSE)
        scgammas2[i] <- - sum(W2[, i] * d - c((p.byt * kk) %*% wGH))
    }
    scalpha <- - sum((p.byt * (d * Y - exp.eta.tw.P * rowsum(Int * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)    
    c(scgammas2, scgammas1, scalpha)
}


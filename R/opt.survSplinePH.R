opt.survSplinePH <-
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
    log.hazard <- eta.t
    log.survival <- - exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + eta.s), id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- d * log.hazard + log.survival    
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}


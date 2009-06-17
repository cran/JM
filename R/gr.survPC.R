gr.survPC <-
function (thetas) {
    gammas <- thetas[seq_len(ncww)]
    alpha <- thetas[ncww + 1]
    xi <- exp(thetas[seq(ncww + 2, length(thetas))])
    eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else rep(0, n)
    eta.t <- eta.tw + alpha * Y
    exp.eta.s <- exp(alpha * Ys)
    exp.eta.tw <- exp(eta.tw)
    Int <- xi[ind.K] * wkP * exp.eta.s
    Int2 <- wkP * exp.eta.s
    scgammas <- if (!is.null(WW))
        - colSums(WW * (d - c((p.byt * (exp.eta.tw * rowsum(Int, id.GK, reorder = FALSE))) %*% wGH)), na.rm = TRUE)
    else 
        NULL
    scalpha <- - sum((p.byt * (d * Y - exp.eta.tw * rowsum(Int * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    scxi <- numeric(Q)
    for (i in 1:Q) {
        i1 <- ind.D == i
        i2 <- ind.K == i
        i3 <- ind.D >= i
        ki <- c((p.byt[i3, ] * (exp.eta.tw[i3] * rowsum(Int2[i2, ], id.GK[i2], reorder = FALSE))) %*% wGH)
        kk <- numeric(n); kk[i3] <- ki
        scxi[i] <- - xi[i] * sum((d * i1)/xi[i] - kk)
    }
    c(scgammas, scalpha, scxi)
}


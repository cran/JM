gr.survPC <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    xi <- exp(thetas$log.xi)
    eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else rep(0, n)
    eta.t <- switch(parameterization, "value" = eta.tw + alpha * Y, 
        "slope" = eta.tw + Dalpha * Y.deriv, "both" = eta.tw + alpha * Y + Dalpha * Y.deriv)    
    exp.eta.s <- exp(switch(parameterization, "value" = alpha * Ys, "slope" = Dalpha * Ys.deriv, 
        "both" = alpha * Ys + Dalpha * Ys.deriv))
    exp.eta.tw <- exp(eta.tw)
    Int <-  wkP * exp.eta.s
    Int2 <- xi[ind.K] * Int
    scgammas <- if (!is.null(WW)) {
        - colSums(WW * (d - c((p.byt * (exp.eta.tw * rowsum(Int2, id.GK, reorder = FALSE))) %*% wGH)), na.rm = TRUE)
    } else NULL
    scalpha <- if (parameterization %in% c("value", "both")) {
        - sum((p.byt * (d * Y - exp.eta.tw * rowsum(Int2 * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    } else NULL
    scalpha.D <- if (parameterization %in% c("slope", "both")) {
        - sum((p.byt * (d * Y.deriv - exp.eta.tw * rowsum(Int2 * Ys.deriv, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    } else NULL
    scxi <- numeric(Q)
    for (i in 1:Q) {
        i1 <- ind.D == i
        i2 <- ind.K == i
        i3 <- ind.D >= i
        ki <- c((p.byt[i3, ] * (exp.eta.tw[i3] * rowsum(Int[i2, ], id.GK[i2], reorder = FALSE))) %*% wGH)
        kk <- numeric(n); kk[i3] <- ki
        scxi[i] <- - xi[i] * sum((d * i1)/xi[i] - kk)
    }
    c(scgammas, scalpha, scalpha.D, scxi)
}


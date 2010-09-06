gr.survSplinePH <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    gammas.bs <- thetas$gammas.bs
    eta.tw1 <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, n)
    eta.tw2 <- as.vector(W2 %*% gammas.bs)
    eta.t <- switch(parameterization, "value" = eta.tw2 + eta.tw1 + alpha * Y, 
        "slope" = eta.tw2 + eta.tw1 + Dalpha * Y.deriv, "both" = eta.tw2 + eta.tw1 + alpha * Y + Dalpha * Y.deriv)    
    eta.s <- switch(parameterization, "value" = alpha * Ys, "slope" = Dalpha * Ys.deriv, 
        "both" = alpha * Ys + Dalpha * Ys.deriv)
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
    scalpha <- if (parameterization %in% c("value", "both")) {
        - sum((p.byt * (d * Y - exp.eta.tw.P * rowsum(Int * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    } else NULL
    scalpha.D <- if (parameterization %in% c("slope", "both")) {
        - sum((p.byt * (d * Y.deriv - exp.eta.tw.P * rowsum(Int * Ys.deriv, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    } else NULL
    c(scgammas1, scalpha, scalpha.D, scgammas2)
}


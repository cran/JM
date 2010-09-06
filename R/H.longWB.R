H.longWB <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    if (parameterization %in% c("value", "both")) {
        Ys <- as.vector(Xs %*% betas) + Zsb
        eta.s <- alpha * Ys
    }
    if (parameterization %in% c("slope", "both")) {
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
        eta.s <- if (parameterization == "both") eta.s + Dalpha * Ys.deriv else Dalpha * Ys.deriv
    }
    exp.eta.tw.P <- exp(eta.tw) * P
    H1 <- XtX / sigma^2
    Int <- wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s) 
    H2 <- matrix(0, ncx, ncx)
    for (i in 1:ncx) {
        for (j in i:ncx) {
            XX <- if (parameterization == "value") {
                alpha^2 * Xs[, i] * Xs[, j]
            } else if (parameterization == "slope") {
                if (i %in% indFixed && j %in% indFixed) {
                    ii <- match(i, indFixed)
                    jj <- match(j, indFixed)
                    Dalpha^2 * Xs.deriv[, ii] * Xs.deriv[, jj]
                } else
                    0
            } else {
                if (i %in% indFixed && j %in% indFixed) {
                    ii <- match(i, indFixed)
                    jj <- match(j, indFixed)
                    (alpha * Xs[, i] + Dalpha * Xs.deriv[, ii]) * (alpha * Xs[, j] + Dalpha * Xs.deriv[, jj])
                } else if (i %in% indFixed && !j %in% indFixed) {
                    ii <- match(i, indFixed)
                    (alpha * Xs[, i] + Dalpha * Xs.deriv[, ii]) * (alpha * Xs[, j])
                } else if (!i %in% indFixed && j %in% indFixed) {
                    jj <- match(j, indFixed)
                    (alpha * Xs[, i]) * (alpha * Xs[, j] + Dalpha * Xs.deriv[, jj])
                } else {
                    alpha^2 * Xs[, i] * Xs[, j]
                }
            }
            ki <- exp.eta.tw.P * rowsum(Int * XX, id.GK, reorder = FALSE)
            kii <- c((p.byt * ki) %*% wGH)
            H2[i, j] <- sum(kii, na.rm = TRUE)
        }
    }
    H2[lower.tri(H2)] <- t(H2)[lower.tri(H2)]
    H1 + H2
}


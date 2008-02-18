`gr.survCH` <-
function (thetas) {
    gammas <- thetas[1:ncww]
    alpha <- thetas[ncww + 1]
    gammas[1:nk] <- cumsum(c(gammas[1], exp(gammas[2:nk])))
    eta <- c(WW %*% gammas) + Y * alpha
    sc <- as.vector(S %*% diff(gammas[1:nk]))
    ew1 <- - exp(eta)
    ew2 <- ew1 * Ztime.b
    out <- (d + c((ew1 * p.byt) %*% wGH)) * WW
    out[, 1:nk] <- out[, 1:nk] + d * SS / sc
    out <- colSums(out, na.rm = TRUE)
    out[1:nk] <- c(out[1:nk] %*% jacobian(thetas[1:nk]))
    sc.alpha <- sum(((d * Y + eta.yxT * ew1 + ew2) * p.byt) %*% wGH, na.rm = TRUE)
    out <- c(out, sc.alpha)
    -out    
}


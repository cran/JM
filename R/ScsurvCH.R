`ScsurvCH` <-
function (thetas) {
    gamas <- thetas[-(ncww + 1)]
    alpha <- thetas[ncww + 1]
    gamas[1:nk] <- cumsum(c(gamas[1], exp(gamas[2:nk])))
    Y <- eta.yxT + rowSums(Ztime * b)
    eta <- as.vector(WW %*% gamas + Y * alpha)
    sc <- as.vector(S %*% diff(gamas[1:nk]))
    ew1 <- - exp(eta)
    ew2 <- ew1 - trc.t1
    ew3 <- ew1 * rowSums(Ztime * b) - trc.t2
    out <- (d + ew2) * WW
    out[, 1:nk] <- out[, 1:nk] + d * SS / sc
    out <- colSums(out, na.rm = TRUE)
    out[1:nk] <- c(out[1:nk] %*% jacobian(thetas[1:nk]))
    out <- c(out, sum(d * Y + eta.yxT * ew1 + ew3, na.rm = TRUE))
    -out
}


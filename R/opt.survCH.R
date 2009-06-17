opt.survCH <-
function (thetas) {
    gammas <- thetas[1:ncww]
    alpha <- thetas[ncww + 1]
    gammas[1:nk] <- cumsum(c(gammas[1], exp(gammas[2:nk])))
    eta <- c(WW %*% gammas) + Y * alpha
    sc <- as.vector(S %*% diff(gammas[1:nk]))
    log.p.tb <- d * (log(sc) + eta - logT) - exp(eta)
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}


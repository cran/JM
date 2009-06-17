opt.survPC <-
function (thetas) {
    gammas <- thetas[seq_len(ncww)]
    alpha <- thetas[ncww + 1]
    xi <- exp(thetas[seq(ncww + 2, length(thetas))])
    eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else 0
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Ys
    log.hazard <- log(xi[ind.D]) + eta.t
    log.survival <- - exp(eta.tw) * rowsum(xi[ind.K] * wkP * exp(eta.s), id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- d * log.hazard + log.survival    
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}


`opt.survPH` <-
function (thetas) {
     if (is.null(WW)) {
        eta.t <- thetas * Y
        eta.t2 <- thetas * Y2
    } else {
        nth <- length(thetas)
        eta.t <- as.vector(WW %*% thetas[-nth]) + thetas[nth] * Y
        eta.t2 <- as.vector(WW %*% thetas[-nth])[indT] + thetas[nth] * Y2
    }
    ew <- exp(eta.t2)
    log.p.tb <- d * (log(lambda0T) + eta.t)
    log.p.tb[ind0, ] <- 0
    log.p.tb[ind.lenN0, ] <- log.p.tb[ind.lenN0, ] - rowsum(lambda0[ind.lambda] * ew, indT)
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}


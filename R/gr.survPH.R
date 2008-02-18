`gr.survPH` <-
function (WW, Y, Y2, ew) {
    sc.gammas <- if (!is.null(WW)) {
        nw <- ncol(WW)
        res <- numeric(nw)
        for (i in 1:nw) {
            p <- rowsum(lambda0. * WW[indT, i] * ew, indT)
            res[i] <- sum((p * p.byt.) %*% wGH)
        }
        crossprod(WW, d) - res
    } else
        NULL
    p1 <- sum(d * c((Y * p.byt) %*% wGH))
    p2 <- rowsum(lambda0. * Y2 * ew, indT)
    sc.alpha <- p1 - sum((p2 * p.byt.) %*% wGH, na.rm = TRUE)
    - c(sc.gammas, sc.alpha)
}


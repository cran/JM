LogLik.chLaplace <-
function (thetas, b) {
    betas <- thetas[1:ncx]
    sigma <- exp(thetas[ncx + 1])
    gammas <- thetas[seq(ncx + 2, ncx + 1 + ncww)]
    gammas[1:nk] <- cumsum(c(gammas[1], exp(gammas[2:nk])))
    alpha <- thetas[ncx + ncww + 2]
    D <- thetas[seq(ncx + ncww + 3, length(thetas))]
    D <- if (diag.D) exp(D) else chol.transf(D)
    environment(fn.b) <- environment(gr.b) <- environment()
    environment(logsurvCH) <- environment(update.bCH) <- environment()
    - as.vector(update.bCH(b, hes.b, betas, sigma, c(gammas, alpha), D, FALSE))
}


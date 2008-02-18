`dmvnorm` <-
function (x, mu, Sigma, log = FALSE) {
    if (!is.matrix(x))
        x <- rbind(x)
    p <- length(mu)
    ed <- eigen(Sigma, symmetric = TRUE)
    ev <- ed$values
    evec <- ed$vectors
    if (!all(ev >= -1e-06 * abs(ev[1]))) 
        stop("'Sigma' is not positive definite")
    ss <- x - rep(mu, each = nrow(x))
    inv.Sigma <- evec %*% (t(evec) / ev)
    quad <- 0.5 * rowSums((ss %*% inv.Sigma) * ss)
    fact <- - 0.5 * (p * log(2 * pi) + sum(log(ev)))
    if (log)
        as.vector(fact - quad)
    else
        as.vector(exp(fact - quad))
}


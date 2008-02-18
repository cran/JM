`initial.surv` <-
function (logT, d, X, method, control) {
    if (method == "ph-GH") {
        dat <- data.frame(Time = exp(logT), d = d)
        cph <- coxph(Surv(Time, d) ~ X, data = dat)
        coefs <- coef(cph)
        nk <- length(coefs)
        lambda0 <- basehaz(cph, FALSE)
        out <- list(alpha = coefs[nk], lambda0 = lambda0[, 1])
        if (nk > 1)
            out$gammas <- coefs[-nk]
        out
    } else if (method == "weibull-GH") {
        dat <- data.frame(Time = exp(logT), d = d)
        init.fit <- survreg(Surv(Time, d) ~ X, data = dat)
        coefs <- init.fit$coef
        nk <- length(coefs)
        list(gammas = coefs[-nk], alpha = coefs[nk], sigma.t = init.fit$scale)    
    } else if (method == "ch-Laplace" || method == "ch-GH")  {
        dat <- data.frame(Time = exp(logT), d = d)
        init.fit <- survreg(Surv(Time, d) ~ X, data = dat)
        coefs <- - coef(init.fit) / init.fit$scale
        min.x <- min(logT)
        max.x <- max(logT)
        kn <- if (is.null(control$knots)) {
            kk <- seq(0, 1, length.out = control$lng.in.kn + 2)[-c(1, control$lng.in.kn + 2)]
            quantile(logT[d == 1], kk, names = FALSE)
        } else {
            control$knots
        }
        kn <- sort(c(rep(c(min.x, max.x), control$ord), kn))
        W <- splineDesign(kn, logT, ord = control$ord)
        nk <- ncol(W)
        nx <- NCOL(X)
        logH <- coefs[1] + logT / init.fit$scale
        coefs <- c(as.vector(lm.fit(W, logH)$coefficients), coefs[-1])
        list(gammas = c(sort(coefs[1:nk]), if (nx > 1) coefs[seq(nk + 1, nk + nx - 1)] else NULL), 
            alpha = coefs[nk + nx])
    }
}


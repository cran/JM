initial.surv <-
function (Time, d, W, WintF.vl, WintF.sl, id, times, method, 
        parameterization, long = NULL, long.deriv = NULL, extra = NULL) {
    old <- options(warn = (-1))
    on.exit(options(old))
    if (!is.null(long)) {
        long.id <- tapply(long, id, tail, 1)
        if (parameterization == "value") 
            longD.id <- NULL
    }
    if (!is.null(long.deriv)) {
        longD.id <- tapply(long.deriv, id, tail, 1)
        if (parameterization == "slope") 
            long.id <- NULL
    }
    WW <- cbind(W, long.id, longD.id)
    if (method %in% c("Cox-PH-GH", "weibull-PH-GH", "piecewise-PH-GH", "spline-PH-GH", "spline-PH-Laplace")) {
        DD <- data.frame(id = id, Time = Time[id], d = d[id], times = times)
        if (!is.null(long)) {
            DD$long <- long * WintF.vl[id, , drop = FALSE]
            k <- ncol(DD$long)
        }
        if (!is.null(long.deriv)) {
            DD$longD <- long.deriv * WintF.sl[id, , drop = FALSE]
            l <- ncol(DD$longD)
        }
        dW <- as.data.frame(W[id, , drop = FALSE], row.names = row.names(DD))
        if (ncol(dW)) {
            names(dW) <- paste("W", seq_along(dW), sep = "")
            DD <- cbind(DD, dW)
        }
        DD$start <- DD$times
        DD$stop <- unlist(lapply(split(DD[c("id", "start", "Time")], DD$id), 
            function (d) c(d$start[-1], d$Time[1])))
        DD$event <- ave(DD$d, DD$id, FUN = function(x) {
            if (length(x) == 1) {
                x
            } else {
                x[seq(length(x) - 1)] <- 0
                x
            }
        })
        baseCovs <- if (ncol(dW))
            paste("+", paste(names(dW), collapse = " + "))
        else 
            NULL
        form <- switch(parameterization,
            "value" = paste("Surv(start, stop, event) ~", "long", baseCovs),
            "slope" = paste("Surv(start, stop, event) ~", "longD", baseCovs),
            "both" = paste("Surv(start, stop, event) ~", "long + longD", baseCovs)
        )
        form <- as.formula(form)
        cph <- coxph(form, data = DD)
        coefs <- cph$coefficients
        out <- switch(parameterization,
            "value" = list(alpha = coefs[1:k], gammas = coefs[-(1:k)]), 
            "slope" = list(Dalpha = coefs[1:l], gammas = coefs[-(1:l)]),
            "both" = list(alpha = coefs[1:k], Dalpha = coefs[(k+1):(k+l)], 
                gammas = coefs[-(1:(k+l))])
        )
        if (method == "Cox-PH-GH") {
            out$lambda0 <- basehaz(cph, FALSE)$hazard
        }
        if (method == "weibull-PH-GH") {
            dat <- data.frame(Time = Time, d = d)
            init.fit <- survreg(Surv(Time, d) ~ WW, data = dat)
            coefs <- - init.fit$coef / init.fit$scale
            out$gammas <- c(coefs[1], out$gammas)
            out$sigma.t <- 1 / init.fit$scale
        }
        if (method == "piecewise-PH-GH") {
            dat <- data.frame(Time = Time, d = d)
            cph. <- coxph(Surv(Time, d) ~ WW, data = dat, x = TRUE)
            init.fit <- piecewiseExp.ph(cph., knots = extra$control$knots)
            coefs <- init.fit$coef
            out$xi <- exp(coefs[grep("xi", names(coefs))])
        }
        if (method == "spline-PH-GH" || method == "spline-PH-Laplace") {
            dat <- data.frame(Time = Time, d = d)
            init.fit <- survreg(Surv(Time, d) ~ WW, data = dat)
            xi <- 1 / init.fit$scale
            phi <- exp(coefs[1])
            logh <- log(phi * xi * Time^(xi - 1))
            out$gammas.bs <- as.vector(lm.fit(extra$W2, logh)$coefficients)
            out
        }
    }
    if (method == "weibull-AFT-GH") {
        dat <- data.frame(Time = Time, d = d)
        if (!is.null(long.id)) {
            long.id <- c(long.id) * WintF.vl
            k <- ncol(WintF.vl)
        }
        if (!is.null(longD.id)) {
            longD.id <- c(longD.id) * WintF.sl
            l <- ncol(WintF.sl)
        }         
        WW <- cbind(W, long.id, longD.id)
        init.fit <- survreg(Surv(Time, d) ~ WW, data = dat)
        coefs <- - init.fit$coef
        nk <- if (is.null(W)) 1 else ncol(W) + 1
        out <- switch(parameterization, 
            "value" = list(gammas = coefs[1:nk], alpha = coefs[-(1:nk)], sigma.t = 1 / init.fit$scale),
            "slope" = list(gammas = coefs[1:nk], Dalpha = coefs[-(1:nk)], sigma.t = 1 / init.fit$scale),
            "both" = list(gammas = coefs[1:nk], alpha = coefs[seq(nk+1, nk+k)], 
                Dalpha = coefs[-seq(1, nk+k)], sigma.t = 1 / init.fit$scale)
        )
    }
    if (method == "ch-Laplace") {
        dat <- data.frame(Time = Time, d = d)
        init.fit <- survreg(Surv(Time, d) ~ WW, data = dat)
        coefs <- - coef(init.fit) / init.fit$scale
        min.x <- min(logT)
        max.x <- max(logT)
        kn <- if (is.null(extra$control$knots)) {
            kk <- seq(0, 1, length.out = extra$control$lng.in.kn + 2)[-c(1, extra$control$lng.in.kn + 2)]
            quantile(log(Time)[d == 1], kk, names = FALSE)
        } else {
            extra$control$knots
        }
        kn <- sort(c(rep(c(min.x, max.x), extra$control$ord), kn))
        W <- splineDesign(kn, log(Time), ord = extra$control$ord)
        nk <- ncol(W)
        nx <- NCOL(X)
        logH <- coefs[1] + logT / init.fit$scale
        coefs <- c(as.vector(lm.fit(W, logH)$coefficients), coefs[-1])
        out <- list(gammas = c(sort(coefs[1:nk]), if (nx > 1) coefs[seq(nk + 1, nk + nx - 1)] else NULL), 
            alpha = coefs[nk + nx])
    }
    out
}


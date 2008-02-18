`predict.jointModel` <-
function (object, newdata, b = NULL, id.var = "id", process = c("Longitudinal", "Event"), 
    type = c("Marginal", "Subject"), scale = c("survival", "cumulative-Hazard", 
    "log-cumulative-Hazard", "expected-future-lifetime"), survTimes = NULL, M = 1000, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (!is.data.frame(newdata))
        stop("'newdata' must be a data.frame.\n")
    process <- match.arg(process)
    type <- match.arg(type)
    scale <- match.arg(scale)
    method <- object$method
    TermsY <- delete.response(object$termsY)
    TermsT <- delete.response(object$termsT)
    if (process == "Longitudinal") {
        mf <- model.frame(TermsY, data = newdata)
        X <- model.matrix(TermsY, mf)
        betas <- object$coefficients$betas
        if (type == "Marginal") {
            c(X %*% betas)
        } else {
            if (any(!is.character(id.var) | length(id.var) > 1 | !id.var %in% names(newdata)))
                stop("'id.var' must be a character string corresponding to one of the columns of 'newdata'.\n")
            if (is.null(b)) {
                id <- newdata[[id.var]]
                b <- ranef(object)
                if (!all(id %in% rownames(b)))
                    stop("'b' was not supplied or some id's do not much the row names of 'ranef(object)'.\n")
                b <- b[rownames(b) %in% id, , drop = FALSE]
            }
            id <- as.vector(unclass(factor(id)))
            Z <- model.matrix(object$formYz, mf)
            c(X %*% betas) + rowSums(Z * b[id, , drop = FALSE])
        }
    } else {
        if (any(!is.character(id.var) | length(id.var) > 1 | !id.var %in% names(newdata)))
            stop("'id.var' must be a character string corresponding to one of the columns of 'newdata'.\n")
        id <- as.vector(unclass(factor(newdata[[id.var]])))        
        newdata.id <- newdata[tapply(rownames(newdata), id, tail, n = 1), ]
        n <- nrow(newdata.id)
        # components of the event process
        mf <- model.frame(TermsT, data = newdata.id)
        W1 <- model.matrix(TermsT, mf)[, -1, drop = FALSE]
        if (ncol(W1) == 0)
            W1 <- NULL
        gammas <- object$coefficients$gammas
        alpha <- object$coefficients$alpha
        if (type == "Marginal") {
            # event times at which the survival function will be computed
            if (is.null(survTimes)) 
                survTimes <- seq(min(object$y$logT), max(object$y$logT), length.out = 50)
            nst <- length(survTimes)
            # simulate random effects
            D <- object$coefficients$D
            diag.D <- (ncz <- ncol(D)) == 1 & nrow(D) > 1
            b <- mvrnorm(M, rep(0, ncz), if (diag.D) diag(c(D)) else D)
            # components of the longitudinal process
            mf <- model.frame(TermsY, data = newdata.id)
            Xtime <- model.matrix(object$formYx, mf)
            Ztime <- model.matrix(object$formYz, mf)        
            betas <- object$coefficients$betas
            Y <- NA + NA
            fitT <- if (object$method == "ph-GH") {
                # event process; marginal predictions; PH model
                eta <- if (is.null(W1)) Y * alpha else c(W1 %*% gammas) + Y * alpha
                switch(scale,
                    "survival" = NA,
                    "cumulative-Hazard" = NA,
                    "log-cumulative-Hazard" = NA)
            } else if (object$method == "weibull-GH") {
                # event process; marginal predictions; Weibull model
                WW <- if (is.null(W1)) as.matrix(rep(1, nst)) else cbind(1, W1)
                eta <- c(WW %*% gammas) + Y * alpha
                logT.mat <- matrix(survTimes, nrow = n, ncol = nst, byrow = TRUE)
                switch(scale,
                    "survival" = exp(- exp((logT.mat - eta) / object$coefficients$sigma.t)),
                    "cumulative-Hazard" = exp((logT.mat - eta) / object$coefficients$sigma.t),
                    "log-cumulative-Hazard" = (logT.mat - eta) / object$coefficients$sigma.t)
            } else {
                # event process; marginal predictions; CH models
                W2 <- splineDesign(object$knots, survTimes, ord = object$control$ord)
                eta <- apply(W2, 1, function (x) {
                    w <- matrix(x, n, length(x), TRUE)
                    WW <- if (is.null(W1)) w else cbind(w, W1)
                    c(WW %*% gammas) + Y * alpha
                })
                switch(scale,
                    "survival" = exp(- exp(eta)),
                    "cumulative-Hazard" = exp(eta),
                    "log-cumulative-Hazard" = eta)
            }
            colMeans(fitT)
        } else {
            # censored times
            mf <- model.frame(object$termsT, data = newdata.id)
            Y <- model.extract(mf, "response")
            Time <- Y[, 1]
            d <- Y[, 2]
            if (any(d == 1))
                stop("'newdata' should only contain subjects that did not exhibit the event.\n")
            # components of the longitudinal process
            mf <- model.frame(object$termsY, data = newdata)
            y <- model.response(mf, "numeric")
            X <- model.matrix(object$formYx, mf)
            Z <- model.matrix(object$formYz, mf)
            mf <- model.frame(object$termsY, data = newdata.id)
            Xtime <- model.matrix(object$formYx, mf)
            Ztime <- model.matrix(object$formYz, mf)        
            betas <- object$coefficients$betas
            sigma <- object$coefficients$sigma
            # predict random effects and calculate Y
            b <- predict.b(method, y, X, Xtime, Z, Ztime, betas, sigma, Time, W1, gammas, alpha, 
                sigma.t = object$coefficients$sigma.t, D = object$coefficients$D, id, control = object$control, 
                knots = object$knots)
            Y <- c(Xtime %*% betas) + rowSums(Ztime * b)
            if (object$method == "ph-GH") {
                # event process; subject-specific predictions; PH model
                NA
                NA
                NA
            } else if (object$method == "weibull-GH") {
                # event process; subject-specific predictions; Weibull model
                WW <- attr(b, "WW")
                sigma.t <- object$coefficients$sigma.t
                eta <- c(WW %*% gammas) + Y * alpha
                S0 <- exp(- exp((log(Time) - eta) / sigma.t))
                f.t <- function (u, eta.i, t0) {
                    u. <- u + t0
                    w <- (log(u.) - eta.i) / sigma.t
                    u * exp(w - exp(w)) / (u. * sigma.t)
                }
                switch(scale,
                    "survival" = S0,
                    "cumulative-Hazard" = - log(S0),
                    "log-cumulative-Hazard" = log(- log(S0)),
                    "expected-future-lifetime" = sapply(1:length(Time), function (i) {
                        integrate(f.t, 0, Inf, eta.i = eta[i], t0 = Time[i])$value
                    }) / S0)
            } else {
                # event process; subject-specific predictions; CH models
                WW <- attr(b, "WW")
                nk <- ncol(WW) - ncol(W1)
                eta <- c(WW %*% gammas) + Y * alpha
                S0 <- exp(- exp(eta))
                f.t <- function (u, eta.i, t0) {
                    u. <- u + t0
                    S <- splineDesign(object$knots[-c(1, length(object$knots))], log(u.), 
                        ord = object$control$ord - 1, outer.ok = TRUE)
                    S <- object$control$ord * S / rep(diff(object$knots, lag = object$control$ord + 1), 
                        each = nrow(S))
                    sc <- as.vector(S %*% diff(gammas[1:nk]))
                    ew <- - exp(eta.i)
                    u * (log(sc) + eta.i + ew - log(u.))
                }
                
                f.t(1, eta[1], Time[1])
                
                sapply(1:length(Time), function (i) {
                    integrate(f.t, 0, Inf, eta.i = eta[i], t0 = Time[i])$value
                }) / S0
                
                
                
                switch(scale,
                    "survival" = S0,
                    "cumulative-Hazard" = exp(eta),
                    "log-cumulative-Hazard" = eta,
                    "expected-future-lifetime" = sapply(1:length(Time), function (i) {
                        integrate(f.t, 0, Inf, eta.i = eta[i], t0 = Time[i])$value
                    }) / S0)
            }            
        }
    }
}


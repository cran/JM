fitted.jointModel <-
function (object, process = c("Longitudinal", "Event"), type = c("Marginal", "Subject", "EventTime"), 
    scale = c("survival", "cumulative-Hazard", "log-cumulative-Hazard"), M = 200, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    process <- match.arg(process)
    type <- match.arg(type)
    scale <- match.arg(scale)
    method <- object$method
    if (process == "Longitudinal") {
        fitY <- c(object$x$X %*% object$coefficients$betas)
        names(fitY) <- names(object$y$y)
        if (type == "Marginal") fitY else if (type == "Subject") fitY + object$EB$Zb else {
            fitYEvent <- if (method == "Cox-PH-GH") {
                c(object$x$Xtime2 %*% object$coefficients$betas + object$EB$Ztime2b)
            } else {
                c(object$x$Xtime %*% object$coefficients$betas + object$EB$Ztimeb)
            }
            names(fitYEvent) <- names(object$y$logT)
            fitYEvent
        }
    } else {
        W1 <- object$x$W
        Y <- if (type == "Marginal") {
            D <- object$coefficients$D
            diag.D <- (ncz <- ncol(D)) == 1 & nrow(D) > 1
            b <- mvrnorm(M, rep(0, ncz), if (diag.D) diag(c(D)) else D)
            if (method == "Cox-PH-GH") {
                Zb <- object$x$Ztime2 %*% t(b)
                c(object$x$Xtime2 %*% object$coefficients$betas) + Zb
            } else {
                Zb <- object$x$Ztime %*% t(b)
                c(object$x$Xtime %*% object$coefficients$betas) + Zb
            }
        } else {
            if (method == "Cox-PH-GH") {
                c(object$x$Xtime2 %*% object$coefficients$betas + object$EB$Ztime2b)
            } else {
                c(object$x$Xtime %*% object$coefficients$betas + object$EB$Ztimeb)
            }
        }
        gammas <- object$coefficients$gammas
        alpha <- object$coefficients$alpha
        logT <- object$y$logT
        fitT <- if (method == "Cox-PH-GH") {
            indT <- object$indexes$indT
            lambda0 <- object$coefficients$lambda0[, "basehaz"]
            eta.tw <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, )
            Ztime2b <- if (type == "Marginal") {
                object$x$Ztime2 %*% t(b)
            } else {
                rowSums(object$x$Ztime2 * object$EB$post.b[indT, ])
            }
            Y2 <- c(object$x$Xtime2 %*% object$coefficients$betas) + Ztime2b
            eta.s <- object$coefficients$alpha * Y2
            if (type == "Marginal") {
                S <- matrix(0, length(logT), ncol(eta.s))
                S[unique(indT), ] <- rowsum(lambda0[object$indexes$ind.L1] * exp(eta.s), indT, reorder = FALSE)
            } else {
                S <- numeric(length(logT))
                S[unique(indT)] <- tapply(lambda0[object$indexes$ind.L1] * exp(eta.s), indT, sum)
            }
            Haz <- exp(eta.tw) * S
            switch(scale,
                "survival" = exp(- Haz),
                "cumulative-Hazard" = Haz,
                "log-cumulative-Hazard" = log(Haz))
        } else if (method == "weibull-PH-GH") {
            WW <- if (is.null(W1)) as.matrix(rep(1, length(logT))) else cbind(1, W1)
            eta.tw <- as.vector(WW %*% gammas)
            sigma.t <- object$coefficients$sigma.t
            P <- object$x$P
            log.st <- log(object$x$st)
            wk <- rep(object$x$wk, length(logT))
            id.GK <- rep(seq_along(logT), each = object$control$GKk)
            Zsb <- if (type == "Marginal") {
                object$x$Zs %*% t(b)
            } else {
                rowSums(object$x$Zs * object$EB$post.b[id.GK, ])
            }
            Ys <- c(object$x$Xs %*% object$coefficients$betas) + Zsb
            eta.s <- object$coefficients$alpha * Ys
            Haz <- exp(eta.tw) * P * rowsum(wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s), id.GK, reorder = FALSE)
            switch(scale,
                "survival" = exp(- Haz),
                "cumulative-Hazard" = Haz,
                "log-cumulative-Hazard" = log(Haz))
        } else if (method == "weibull-AFT-GH") {
            WW <- if (is.null(W1)) as.matrix(rep(1, length(logT))) else cbind(1, W1)
            eta.tw <- as.vector(WW %*% gammas)
            sigma.t <- object$coefficients$sigma.t
            P <- object$x$P
            log.st <- log(object$x$st)
            wk <- rep(object$x$wk, length(logT))
            id.GK <- rep(seq_along(logT), each = object$control$GKk)
            Zsb <- if (type == "Marginal") {
                object$x$Zs %*% t(b)
            } else {
                rowSums(object$x$Zs * object$EB$post.b[id.GK, ])
            }
            Ys <- c(object$x$Xs %*% object$coefficients$betas) + Zsb
            eta.s <- object$coefficients$alpha * Ys
            Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL            
            Haz <- Vi^sigma.t
            switch(scale,
                "survival" = exp(- Haz),
                "cumulative-Hazard" = Haz,
                "log-cumulative-Hazard" = log(Haz))
        } else if (method == "spline-PH-GH") {
            eta.tw <- if (!is.null(W1)) as.vector(W1 %*% gammas) else 0
            P <- object$x$P
            wk <- rep(object$x$wk, length(logT))
            id.GK <- rep(seq_along(logT), each = object$control$GKk)
            Zsb <- if (type == "Marginal") {
                object$x$Zs %*% t(b)
            } else {
                rowSums(object$x$Zs * object$EB$post.b[id.GK, ])
            }
            Ys <- c(object$x$Xs %*% object$coefficients$betas) + Zsb
            eta.s <- object$coefficients$alpha * Ys
            Haz <- exp(eta.tw) * P * rowsum(wk * exp(c(object$x$W2s %*% object$coefficients$gammas.bs) + eta.s), id.GK, reorder = FALSE)
            switch(scale,
                "survival" = exp(- Haz),
                "cumulative-Hazard" = Haz,
                "log-cumulative-Hazard" = log(Haz))            
        } else if (method == "piecewise-PH-GH") {
            WW <- W1
            eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else 0
            xi <- object$coefficients$xi
            nk <- object$control$GKk
            id.GK <- object$x$id.GK
            ind.K <- rep(unlist(lapply(object$y$ind.D, seq_len)), each = nk)
            wk <- unlist(lapply(object$y$ind.D, function (n) rep(object$x$wk, n)))
            P <- c(t(object$x$P))
            wkP <- wk * rep(P[!is.na(P)], each = nk)
            Zsb <- if (type == "Marginal") {
                object$x$Zs %*% t(b)
            } else {
                rowSums(object$x$Zs * object$EB$post.b[id.GK, ])
            }
            Ys <- c(object$x$Xs %*% object$coefficients$betas) + Zsb
            eta.s <- object$coefficients$alpha * Ys
            Haz <- exp(eta.tw) * rowsum(xi[ind.K] * wkP * exp(eta.s), id.GK, reorder = FALSE)
            switch(scale,
                "survival" = exp(- Haz),
                "cumulative-Hazard" = Haz,
                "log-cumulative-Hazard" = log(Haz))
        } else {
            W2 <- splineDesign(object$knots, logT, ord = object$control$ord)
            WW <- if (is.null(W1)) W2 else cbind(W2, W1)
            eta <- c(WW %*% gammas) + Y * alpha
             switch(scale,
                "survival" = exp(- exp(eta)),
                "cumulative-Hazard" = exp(eta),
                "log-cumulative-Hazard" = eta)
        }
        fitT <- if (type == "Marginal") rowMeans(fitT) else c(fitT)
        names(fitT) <- names(logT)
        fitT
    }
}


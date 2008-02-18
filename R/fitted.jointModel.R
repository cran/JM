`fitted.jointModel` <-
function (object, process = c("Longitudinal", "Event"), type = c("Marginal", "Subject"), 
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
        if (type == "Marginal") fitY else fitY + object$EB$Zb
    } else {
        W1 <- object$x$W
        Y <- if (type == "Marginal") {
            D <- object$coefficients$D
            diag.D <- (ncz <- ncol(D)) == 1 & nrow(D) > 1
            b <- mvrnorm(M, rep(0, ncz), if (diag.D) diag(c(D)) else D)
            if (method == "ph-GH") {
                Zb <- object$x$Ztime2 %*% t(b)
                c(object$x$Xtime2 %*% object$coefficients$betas) + Zb
            } else {
                Zb <- object$x$Ztime %*% t(b)
                c(object$x$Xtime %*% object$coefficients$betas) + Zb            
            }
        } else {
            if (method == "ph-GH") {
                c(object$x$Xtime2 %*% object$coefficients$betas + object$EB$Ztime2b)
            } else {
                c(object$x$Xtime %*% object$coefficients$betas + object$EB$Ztimeb)
            }
        }
        gammas <- object$coefficients$gammas
        alpha <- object$coefficients$alpha
        logT <- object$y$logT
        fitT <- if (method == "ph-GH") {
            indT <- object$indexes$indT
            ind.lambda <- object$indexes$ind.lambda
            lambda0 <- object$coefficients$lambda0[, "basehaz"]
            eta.t2 <- if (is.null(W1)) alpha * Y else as.vector(W1 %*% gammas)[indT] + alpha * Y
            ew <- exp(eta.t2)
            S <- exp(- rowsum(lambda0[ind.lambda] * ew, indT))
            switch(scale,
                "survival" = S,
                "cumulative-Hazard" = - log(S),
                "log-cumulative-Hazard" = log(- log(S)))
        } else if (method == "weibull-GH") {
            WW <- if (is.null(W1)) as.matrix(rep(1, length(logT))) else cbind(1, W1)
            eta <- c(WW %*% gammas) + Y * alpha
            switch(scale,
                "survival" = exp(- exp((logT - eta) / object$coefficients$sigma.t)),
                "cumulative-Hazard" = exp((logT - eta) / object$coefficients$sigma.t),
                "log-cumulative-Hazard" = (logT - eta) / object$coefficients$sigma.t)
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


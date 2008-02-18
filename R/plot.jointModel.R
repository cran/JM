`plot.jointModel` <-
function (x, which = 1:4, caption = c("Residuals vs Fitted", "Normal Q-Q", "Marginal Survival", 
    "Marginal Cumulative Hazard", "Marginal log Cumulative Hazard", "Baseline Hazard", "Cumulative Baseline Hazard", 
    "Subject-specific Survival", "Subject-specific Cumulative Hazard", "Subject-specific log Cumulative Hazard"), 
    survTimes = NULL, main = "", ask = prod(par("mfcol")) < length(which) && dev.interactive(), ..., ids = NULL, 
    add.smooth = getOption("add.smooth"), add.qqline = TRUE, add.KM = FALSE, cex.caption = 1) {
    if (!inherits(x, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (!is.numeric(which) || any(which < 1) || any(which > 10)) 
        stop("'which' must be in 1:10.\n")
    show <- rep(FALSE, 10)
    show[which] <- TRUE
    method <- x$method
    if (any(show[6], show[7]) && method != "ph-GH") {
        show[6] <- show[7] <- FALSE
        warning("the baseline hazard and the cumulative baseline hazard are only plotted for the 'ph-GH' method.\n")
    }   
    if (any(show[c(3:5, 8:10)])) {
        if (is.null(ids))
            ids <- seq_along(x$y$logT)
        if (is.null(survTimes) || !is.numeric(survTimes))
            survTimes <- seq(min(exp(x$y$logT)), max(exp(x$y$logT)), length.out = 15)
        log.survTimes <- log(survTimes)
        nt <- length(survTimes)
        n <- x$n
        W1 <- x$x$W
        gammas <- x$coefficients$gammas
        alpha <- x$coefficients$alpha
        fitT <- if (method == "ph-GH") {
            lambda0 <- x$coefficients$lambda0[, "basehaz"]
            unqT <- x$coefficients$lambda0[, "time"]
            times <- lapply(survTimes, function (t) unqT[t >= unqT])
            ind.lambda <- unlist(sapply(sapply(times, length), function (x) {
                if (x > 0) seq(1, x) else NULL
            }), use.names = FALSE)
            indL <- rep(seq_along(survTimes), sapply(times, length))
            times <- unlist(times, use.names = FALSE)            
            indT <- rep(1:x$n, each = length(times))
            data.id2 <- x$data.id[indT, ]
            data.id2[x$timeVar] <- unlist(times, use.names = FALSE)
            mf <- model.frame(x$termsY, data = data.id2)
            Xtime2 <- model.matrix(x$formYx, mf)
            Ztime2 <- model.matrix(x$formYz, mf)
            Y <- c(Xtime2 %*% x$coefficients$betas + rowSums(Ztime2 * x$EB$post.b[indT, ]))
            eta <- if (is.null(W1)) alpha * Y else as.vector(W1 %*% gammas)[indT] + alpha * Y
            ew <- exp(eta)
            indL <- rep(indL, x$n)
            ind.lambda <- rep(ind.lambda, x$n)
            ff <- paste(indL, "\t", indT)
            H <- rowsum(lambda0[ind.lambda] * ew, ff, reorder = FALSE)            
            dim(H) <- c(nt, n)
            H <- t(H)
            list("survival" = exp(- H), "cumulative-Hazard" = H, "log-cumulative-Hazard" = log(H))
        } else if (method == "weibull-GH") {
            WW <- if (is.null(W1)) as.matrix(rep(1, n)) else cbind(1, W1)
            Y <- c(x$x$Xtime %*% x$coefficients$betas + x$EB$Ztimeb)
            eta <- c(WW %*% gammas) + Y * alpha
            logT.mat <- matrix(log.survTimes, nrow = n, ncol = nt, byrow = TRUE)
            sigma.t <- x$coefficients$sigma.t
            list("survival" = exp(- exp((logT.mat - eta) / sigma.t)),
                 "cumulative-Hazard" = exp((logT.mat - eta) / sigma.t),
                 "log-cumulative-Hazard" = (logT.mat - eta) / sigma.t)
        } else {
            W2 <- splineDesign(x$knots, log.survTimes, ord = x$control$ord)
            Y <- c(x$x$Xtime %*% x$coefficients$betas + x$EB$Ztimeb)
            eta <- apply(W2, 1, function (x) {
                w <- matrix(x, n, length(x), TRUE)
                WW <- if (is.null(W1)) w else cbind(w, W1)
                c(WW %*% gammas) + Y * alpha
            })
            list("survival" = exp(- exp(eta)), "cumulative-Hazard" = exp(eta), "log-cumulative-Hazard" = eta)
        }
    }
    one.fig <- prod(par("mfcol")) == 1   
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
        fitY <- fitted(x, process = "Longitudinal", type = "Subject")
        resY <- residuals(x, process = "Longitudinal", type = "Subject")
        plot(fitY, resY, xlab = "Fitted Values", ylab = "Residuals", main = main, ...)
        if (add.smooth) {
            abline(h = 0, lty = 3, col = "grey", lwd = 2)
            panel.smooth(fitY, resY, lwd = 2)
        }
        mtext(caption[1], 3, 0.25, cex = cex.caption)
    } 
    if (show[2]) {
        resY <- residuals(x, process = "Longitudinal", type = "stand-Subject")
        qqnorm(resY, ylab = "Standardized Residuals", main = main, ...)
        if (add.qqline)
            qqline(resY, lty = 3, col = "grey50")
        mtext(caption[2], 3, 0.25, cex = cex.caption)
    }
    if (show[3]) {
        yy <- colMeans(fitT[["survival"]])
        if (add.KM) {
            Time <- exp(x$y$logT)
            failure <- x$y$d
            sf <- survfit(Surv(Time, failure) ~ 1)
            plot(sf, xlab = "Time", ylab = "Survival", main = main, mark.time = FALSE)
            lines(survTimes, yy, ...)
        } else {
            plot(survTimes, yy, xlab = "Time", ylab = "Survival", main = main, ylim = c(0, 1), type = "l", ...)
        }
        mtext(caption[3], 3, 0.25, cex = cex.caption)
    }
    if (show[4]) {
        yy <- colMeans(fitT[["cumulative-Hazard"]])
        plot(survTimes, yy, xlab = "Time", ylab = "Cumulative Hazard", main = main, type = "l", ...)
        mtext(caption[4], 3, 0.25, cex = cex.caption)
    }
    if (show[5]) {
        yy <- colMeans(fitT[["log-cumulative-Hazard"]])
        plot(survTimes, yy, xlab = "Time", ylab = "log Cumulative Hazard", main = main, type = "l", ...)
        mtext(caption[5], 3, 0.25, cex = cex.caption)
    }
    if (show[6]) {
        lambda0 <- x$coefficients$lambda0
        plot(lambda0[, "time"], lambda0[, "basehaz"], xlab = "Time", ylab = "", main = main, ...)
        if (add.smooth) {
            panel.smooth(lambda0[, "time"], lambda0[, "basehaz"], lwd = 2)
        }
        mtext(caption[6], 3, 0.25, cex = cex.caption)
    }
    if (show[7]) {
        lambda0 <- x$coefficients$lambda0
        plot(lambda0[, "time"], cumsum(lambda0[, "basehaz"]), xlab = "Time", ylab = "", main = main, type = "s", ...)
        mtext(caption[7], 3, 0.25, cex = cex.caption)
    }
    if (show[8]) {
        yy <- t(fitT[["survival"]])
        matplot(survTimes, yy[, ids], type = "l", col = "black", lty = 1, 
            xlab = "Time", ylab = "Survival", main = main)
        mtext(caption[8], 3, 0.25, cex = cex.caption)
    }
    if (show[9]) {
        yy <- t(fitT[["cumulative-Hazard"]])
        matplot(survTimes, yy[, ids], type = "l", col = "black", lty = 1, 
            xlab = "Time", ylab = "Cumulative Hazard", main = main)
        mtext(caption[9], 3, 0.25, cex = cex.caption)
    }
    if (show[10]) {
        yy <- t(fitT[["log-cumulative-Hazard"]])
        matplot(survTimes, yy[, ids], type = "l", col = "black", lty = 1, 
            xlab = "Time", ylab = "log Cumulative Hazard", main = main)
        mtext(caption[10], 3, 0.25, cex = cex.caption)
    }
    invisible()
}


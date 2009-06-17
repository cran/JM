plot.jointModel <-
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
            survTimes <- seq(min(exp(x$y$logT)), max(exp(x$y$logT)), length.out = 31)
        log.survTimes <- log(survTimes)
        nt <- length(survTimes)
        n <- x$n
        W1 <- x$x$W
        gammas <- x$coefficients$gammas
        alpha <- x$coefficients$alpha
        fitT <- if (method == "ph-GH") {
            lambda0 <- x$coefficients$lambda0[, "basehaz"]
            unqT <- x$coefficients$lambda0[, "time"]
            T.mat <- matrix(exp(log.survTimes), nrow = n, ncol = nt, byrow = TRUE)
            eta.tw <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, n)
            Haz <- matrix(0, n, nt)
            for (i in 1:nt) {
                times <- lapply(T.mat[, i], function (t) unqT[t >= unqT])
                ind.len <- sapply(times, length)
                indT <- rep(1:nrow(x$data.id), ind.len)
                data.id2 <- x$data.id[indT, ]
                data.id2[x$timeVar] <- unlist(times, use.names = FALSE)
                mf <- model.frame(x$termsY, data = data.id2)
                Xtime2 <- model.matrix(x$formYx, mf)
                Ztime2 <- model.matrix(x$formYz, mf)
                nk <- as.vector(sapply(split(indT, indT), length))
                ind.L1 <- unlist(lapply(nk, seq, from = 1))
                Y2 <- c(Xtime2 %*% x$coefficients$betas + rowSums(Ztime2 * x$EB$post.b[indT, ]))
                eta.s <- alpha * Y2
                S <- numeric(n)
                S[unique(indT)] <- tapply(lambda0[ind.L1] * exp(eta.s), indT, sum)
                Haz[, i] <- exp(eta.tw) * S
            }
            list("survival" = exp(- Haz), "cumulative-Hazard" = Haz, "log-cumulative-Hazard" = log(Haz))
        } else if (method == "weibull-PH-GH") {
            T.mat <- matrix(exp(log.survTimes), nrow = n, ncol = nt, byrow = TRUE)
            WW <- if (is.null(W1)) as.matrix(rep(1, n)) else cbind(1, W1)
            eta.tw <- as.vector(WW %*% gammas)
            sigma.t <- x$coefficients$sigma.t
            b <- x$EB$post.b
            wk <- gaussKronrod()$wk
            sk <- gaussKronrod()$sk
            id.GK <- rep(seq_len(n), each = x$control$GKk)
            Haz <- matrix(0, n, nt)
            for (i in 1:nt) {
                P <- T.mat[, i] / 2
                st <- outer(P, sk + 1)
                data.id <- x$data.id[id.GK, ]
                data.id[x$timeVar] <- c(t(st))
                mf <- model.frame(x$termsY, data = data.id)
                Xs <- model.matrix(x$formYx, mf)
                Zs <- model.matrix(x$formYz, mf)
                log.st <- log(c(t(st)))
                Ys <- c(Xs %*% x$coefficients$betas) + rowSums(Zs * b[id.GK, , drop = FALSE])
                eta.s <- alpha * Ys
                Haz[, i] <- exp(eta.tw) * P * rowsum(wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s), id.GK, reorder = FALSE)
            }
            list("survival" = exp(- Haz),
                 "cumulative-Hazard" = Haz,
                 "log-cumulative-Hazard" = log(Haz))
        } else if (method == "weibull-AFT-GH") {
            T.mat <- matrix(exp(log.survTimes), nrow = n, ncol = nt, byrow = TRUE)
            WW <- if (is.null(W1)) as.matrix(rep(1, n)) else cbind(1, W1)
            eta.tw <- as.vector(WW %*% gammas)
            sigma.t <- x$coefficients$sigma.t
            b <- x$EB$post.b
            wk <- gaussKronrod()$wk
            sk <- gaussKronrod()$sk
            id.GK <- rep(seq_len(n), each = x$control$GKk)
            Haz <- matrix(0, n, nt)
            for (i in 1:nt) {
                P <- T.mat[, i] / 2
                st <- outer(P, sk + 1)
                data.id <- x$data.id[id.GK, ]
                data.id[x$timeVar] <- c(t(st))
                mf <- model.frame(x$termsY, data = data.id)
                Xs <- model.matrix(x$formYx, mf)
                Zs <- model.matrix(x$formYz, mf)
                log.st <- log(c(t(st)))
                Ys <- c(Xs %*% x$coefficients$betas) + rowSums(Zs * b[id.GK, , drop = FALSE])
                eta.s <- alpha * Ys
                Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
                Haz[, i] <- Vi^sigma.t
            }
            list("survival" = exp(- Haz),
                 "cumulative-Hazard" = Haz,
                 "log-cumulative-Hazard" = log(Haz))
        } else if (method == "piecewise-PH-GH") {
            T.mat <- matrix(exp(log.survTimes), nrow = n, ncol = nt, byrow = TRUE)
            Q <- x$x$Q
            qs <- c(0, x$control$knots, max(exp(x$y$logT)) + 1)
            WW <- W1
            eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else 0
            xi <- x$coefficients$xi
            b <- x$EB$post.b
            sk <- gaussKronrod(x$control$GKk)$sk
            nk <- length(sk)
            Haz <- matrix(0, n, nt)
            for (i in 1:nt) {
                ind.D <- findInterval(T.mat[, i], qs, rightmost.closed = TRUE)
                Tiq <- outer(T.mat[, i], qs, pmin)
                Lo <- Tiq[, 1:Q]
                Up <- Tiq[, 2:(Q+1)]
                T <- Up - Lo
                P <- T / 2
                P[P < x$control$tol3] <- as.numeric(NA)
                P1 <- (Up + Lo) / 2
                st <- matrix(0, n, nk*Q)
                skQ <- rep(sk, Q)
                for (ii in 1:n) {
                    st[ii, ] <- rep(P[ii, ], each = nk) * skQ + rep(P1[ii, ], each = nk)
                }
                data.id2 <- x$data.id[rep(1:n, each = nk*Q), ]
                data.id2[x$timeVar] <- c(t(st))
                mf <- model.frame(x$termsY, data = data.id2)
                Xs <- model.matrix(x$formYx, mf)
                Zs <- model.matrix(x$formYz, mf)
                id.GK <- rep(1:n, rowSums(!is.na(st)))
                Ys <- c(Xs %*% x$coefficients$betas) + rowSums(Zs * b[id.GK, , drop = FALSE])
                eta.s <- alpha * Ys
                ind.K <- rep(unlist(lapply(ind.D, seq_len)), each = nk)
                wk <- unlist(lapply(ind.D, function (n) rep(gaussKronrod(x$control$GKk)$wk, n)))
                P <- c(t(P))
                wkP <- wk * rep(P[!is.na(P)], each = nk)
                Haz[, i] <- exp(eta.tw) * rowsum(xi[ind.K] * wkP * exp(eta.s), id.GK, reorder = FALSE)
            }
            list("survival" = exp(- Haz),
                 "cumulative-Hazard" = Haz,
                 "log-cumulative-Hazard" = log(Haz))
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


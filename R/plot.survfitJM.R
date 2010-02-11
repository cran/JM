plot.survfitJM <-
function (x, estimator = c("both", "mean", "median"), which = NULL, fun = NULL, conf.int = FALSE,
        add.last.time.axis.tick = FALSE, include.y = FALSE, main = NULL, xlab = NULL, ylab = NULL, 
        lty = NULL, col = NULL, lwd = NULL, pch = NULL, ask = NULL, legend = FALSE, ...) {
    estimator <- match.arg(estimator)
    fun <- if (!is.null(fun)) match.fun(fun)
    if (is.null(which))
        which <- seq_along(x$summaries)
    if (conf.int && is.null(x$success.rate)) {
        warning("\na confidence interval can be included only when argument 'simulate' of survfitJM() was set to TRUE.")
        conf.int <- FALSE
    }
    if (is.null(ask))
        ask <- prod(par("mfcol")) < length(which)
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (is.null(main)) {
        main <- paste("Subject", names(x$summaries))
        names(main) <- names(x$summaries)
    }
    if (is.null(xlab))
        xlab <- rep("Time", length(which))
    if (is.null(ylab))
        ylab <- if (is.null(fun))
            rep(expression(paste("Pr(", T[i] >= u, " | ", T[i] > t, ", ", tilde(y)[i](t), ")", sep = " ")), length(which))
        else
            rep("", length(which))
    if (!is.null(x$success.rate)) {
    if (is.null(col))
        col <- switch(estimator, both = c(2, 3, 1, 1), mean = c(2, 1, 1), median = c(3, 1, 1))
    if (is.null(lty))
        lty <- switch(estimator, both = c(1, 1, 2, 2), mean = c(1, 2, 2), median = c(1, 2, 2))
    if (is.null(lwd))
        lwd <- switch(estimator, both = c(1, 1, 1, 1), mean = c(1, 1, 1), median = c(1, 1, 1))
    } else {
        col <- lty <- lwd <- 1
    }
    if (is.null(pch))
        pch <- 8
    for (i in seq_along(which)) {
        ii <- which[i]
        r <- x$summaries[[ii]]
        r <- if (!is.null(x$success.rate)) {
            rbind(cbind(c(0, x$last.time[ii]), matrix(1, 2, 4)), r)
        } else {
            rbind(cbind(c(0, x$last.time[ii]), matrix(1, 2, 1)), r)
        }
        if (!is.null(fun) && is.function(fun))
            r[, 2:ncol(r)] <- fun(r[, 2:ncol(r)])
        if (!is.null(x$success.rate) && estimator == "mean")
            r <- r[, -2]
        if (!is.null(x$success.rate) && estimator == "median")
            r <- r[, -3]
        if (!conf.int && !is.null(x$success.rate)) {
            exc <- c(ncol(r) - 1, ncol(r))
            r <- r[, -exc, drop = FALSE]
            col <- col[-exc]
            lty <- lty[-exc]
            lwd <- lwd[-exc]
        }
        if (include.y)
            layout(matrix(c(1, 2)))
        ylim <- if (is.null(fun)) c(0, 1) else { rr <- r[, -1, drop = FALSE]; range(rr[is.finite(rr)]) } 
        matplot(r[, 1], r[, -1, drop = FALSE], type = "l", col = col, lwd = lwd, lty = lty, ylim = ylim, 
            main = main[ii], xlab = xlab[i], ylab = ylab[i], ...)
        if (add.last.time.axis.tick)
            axis(1, at = round(x$last.time[ii], 1))
        if (legend) {
            lab <- switch(estimator, both = c("Mean", "Median", "Lower limit", "Upper limit"), 
                mean = c("Mean", "Lower limit", "Upper limit"), median = c("Median", "Lower limit", "Upper limit"))
            legend("left", lab, lwd = lwd, lty = lty, col = col, bty = "n", ...)
        }
        if (include.y)
            plot(x$obs.times[[ii]], x$y[[ii]], xlim = range(x$survTimes), ylim = x$ry,
                xlab = xlab[i], ylab = "Longitudinal Outcome", pch = pch, ...)
    }
    invisible()
}


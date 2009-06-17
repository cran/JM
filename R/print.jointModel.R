print.jointModel <-
function (x, digits = max(4, getOption("digits") - 4), ...) {
    if (!inherits(x, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("Variance Components:\n")
    names(x$coefficients$alpha) <- ""
    D <- x$coefficients$D
    ncz <- nrow(D)
    diag.D <- ncz != ncol(D)
    sds <- if (diag.D) sqrt(D) else sqrt(diag(D))
    if (ncz > 1) {
        if (diag.D) {
            dat <- round(c(D), digits)
            names(dat) <- rownames(D)
        } else {
            corrs <- cov2cor(D)
            corrs[upper.tri(corrs, TRUE)] <- 0
            mat <- round(cbind(sds, corrs[, -ncz]), digits)
            mat <- apply(mat, 2, sprintf, fmt = "% .4f")
            mat[mat == mat[1, 2]] <- ""
            mat[1, -1] <- abbreviate(colnames(mat)[-1], 6)
            colnames(mat) <- c(colnames(mat)[1], rep("", ncz - 1))
            dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
            names(dat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
            row.names(dat) <- dimnames(D)[[1]]
        }
    } else {
        dat <- data.frame("StdDev" = sds, row.names = rownames(D), 
            check.rows = FALSE, check.names = FALSE)
    }
    lis <- list("Random-Effects" = dat, "Residual" = c("Longitudinal" = x$coefficients$sigma, 
        "Event" = x$coefficients$sigma.t))
    print(lapply(lis, function (x) if (!is.numeric(x)) x else round(x, digits = digits)))
    cat("\nCoefficients:\n")
    gammas <- x$coefficients$gammas
    if (x$method %in% c("ch-GH", "ch-Laplace")) {
        ng <- length(gammas)
        nw <- ncol(x$x$W)
        gammas <- if (is.null(nw)) NULL else gammas[seq(ng - nw + 1, ng)]
    }
    gammas <- c(gammas, "Assoct" = as.vector(x$coefficients$alpha), if (x$method == "piecewise-PH-GH") x$coefficients$xi else NULL)
    print(lapply(list("Longitudinal Process" = x$coefficients$betas, "Event Process" = gammas), 
        round, digits = digits))
    cat("\nlog-Lik:", x$logLik)
    cat("\n\n")
    invisible(x)
}


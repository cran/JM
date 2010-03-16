jointModel <-
function (lmeObject, survObject, timeVar, method = c("weibull-AFT-GH", "weibull-PH-GH", "piecewise-PH-GH", 
    "Cox-PH-GH", "spline-PH-GH", "ch-Laplace"), lag = 0, init = NULL, control = list(), ...) {
    cl <- match.call()
    if (!inherits(lmeObject, "lme"))
        stop("\n'lmeObject' must inherit from class lme.")
    if (length(lmeObject$group) > 1)
        stop("\nnested random-effects are not allowed in lme().")
    if (!is.null(lmeObject$modelStruct$corStruct))
        warning("correlation structure in 'lmeObject' is ignored.\n")
    if (!is.null(lmeObject$modelStruct$varStruct))
        warning("variance structure in 'lmeObject' is ignored.\n")        
    if (!inherits(survObject, "coxph") && !inherits(survObject, "survreg"))
        stop("\n'survObject' must inherit from class coxph or class survreg.")
    if (is.null(survObject$x))
        stop("\nuse argument 'x = TRUE' in ", if (inherits(survObject, "coxph")) "'coxph()'." else "'survreg()'.")
    if (length(timeVar) != 1 || !is.character(timeVar))
        stop("\n'timeVar' must be a character string.")
    method <- match.arg(method)
    if (method == "Cox-PH-GH" && !inherits(survObject, "coxph"))
        stop("\nfor 'method = Cox-PH-GH', 'survObject' must inherit from class coxph.")
    # survival process
    formT <- formula(survObject)
    if (inherits(survObject, "coxph")) {
        W <- survObject$x
        Time <- survObject$y[, 1]
    } else {
        W <- survObject$x[, -1, drop = FALSE]
        Time <- exp(survObject$y[, 1])
    }
    nT <- length(Time)
    if (!length(W))
        W <- NULL
    d <- survObject$y[, 2]
    if (sum(d) < 5)
        warning("\nmore than 5 events are required.")
    # longitudinal process
    id <- as.vector(unclass(lmeObject$groups[[1]]))
    b <- data.matrix(ranef(lmeObject))
    dimnames(b) <- NULL
    nY <- nrow(b)
    if (nY != nT)
        stop("sample sizes in the longitudinal and event processes differ.\n")
    Terms <- lmeObject$terms
    data <- lmeObject$data[all.vars(Terms)]
    formYx <- formula(lmeObject)
    formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
    mf <- model.frame(Terms, data = data)
    y.long <- model.response(mf, "numeric")
    X <- model.matrix(formYx, mf)
    Z <- model.matrix(formYz, mf)
    data.id <- data[!duplicated(id), ]
    if (!timeVar %in% names(data))
        stop("\n'timeVar' does not correspond to one of the columns in the model.frame of 'lmeObject'.")
    # control values
    con <- list(only.EM = FALSE, iter.EM = 50, iter.qN = 150, optimizer = "optim", tol1 = 1e-03, tol2 = 1e-04, 
        tol3 = sqrt(.Machine$double.eps), numeriDeriv = "fd", eps.Hes = 1e-06, parscale = NULL, step.max = 0.1, 
        backtrackSteps = 2, knots = NULL, lng.in.kn = if (method == "piecewise-PH-GH") 6 else 5, ord = 4, 
        equal.strata.knots = TRUE, GHk = if (ncol(Z) < 3) 15 else 9, GKk = if (method == "piecewise-PH-GH") 7 else 15, verbose = FALSE)
    if (method == "Cox-PH-GH") {
        con$only.EM <- TRUE
        con$iter.EM <- 200
        con$GHk <- if (ncol(Z) == 1) 15 else if (ncol(Z) == 2) 11 else 9
    }
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0) 
        warning("unknown names in 'control': ", paste(noNms, collapse = ", "))
    if (method == "Cox-PH-GH" && !con$only.EM)
        stop("with method 'Cox-PH-GH' only the EM algorithm is used.\n")
    if (method == "Cox-PH-GH" && any(!is.na(match(c("iter.qN", "optimizer"), namc))))
        warning("method 'Cox-PH-GH' uses only the EM algorithm.\n")
    # extra design matrices for the longitudinal part
    data.id[timeVar] <- pmax(Time - lag, 0)
    mf <- model.frame(lmeObject$terms, data = data.id)
    Xtime <- model.matrix(formYx, mf)
    Ztime <- model.matrix(formYz, mf)
    long <- as.vector(c(Xtime %*% fixef(lmeObject)) + rowSums(Ztime * b))
    # response vectors and design matrices
    y <- list(y = y.long, logT = log(Time), d = d, lag = lag)
    x <- list(X = X, Xtime = Xtime, Z = Z, Ztime = Ztime, W = W)
    # extra design matrices for 'method = "weibull-AFT-GH"' and 'method = "weibull-PH-GH"'
    # extra design matrices for 'method = "spline-PH-GH"' and 'method = "spline-PH-Laplace"'
    if (method %in% c("weibull-AFT-GH", "weibull-PH-GH", "spline-PH-GH", "spline-PH-Laplace")) {
        wk <- gaussKronrod(con$GKk)$wk
        sk <- gaussKronrod(con$GKk)$sk
        P <- Time/2
        st <- outer(P, sk + 1)
        data.id2 <- data.id[rep(seq_len(nY), each = con$GKk), ]
        data.id2[timeVar] <- pmax(c(t(st)) - lag, 0)
        mf <- model.frame(lmeObject$terms, data = data.id2)
        Xs <- model.matrix(formYx, mf)
        Zs <- model.matrix(formYz, mf)
        x <- c(x, list(Xs = Xs, Zs = Zs, P = P, st = c(t(st)), wk = wk))
        if (method == "spline-PH-GH" || method == "spline-PH-Laplace") {
            strt <- if (is.null(survObject$strata)) gl(1, nT) else survObject$strata
            nstrt <- length(levels(strt))
            split.Time <- split(Time, strt)
            kn <- if (con$equal.strata.knots) {
                kk <- if (is.null(con$knots)) {
                    pp <- seq(0, 1, length.out = con$lng.in.kn + 2)[-c(1, con$lng.in.kn + 2)]
                    quantile(Time, pp, names = FALSE)
                } else {
                    con$knots
                }
                rr <- rep(list(sort(c(rep(range(Time, st), con$ord), kk))), nstrt)
                names(rr) <- names(split.Time)
                rr
            } else {
                lapply(split.Time, function (t) {
                    kk <- if (is.null(con$knots)) {
                        pp <- seq(0, 1, length.out = con$lng.in.kn + 2)[-c(1, con$lng.in.kn + 2)]
                        quantile(t, pp, names = FALSE)
                    } else {
                        con$knots
                    }
                    sort(c(rep(range(Time, st), con$ord), kk))
                })
            }
            con$knots <- kn
            W2 <- mapply(function (k, t) splineDesign(k, t, ord = con$ord), kn, split.Time, SIMPLIFY = FALSE)
            W2 <- mapply(function (w2, ind) {
                out <- matrix(0, nT, ncol(w2))
                out[strt == ind, ] <- w2
                out
            }, W2, levels(strt), SIMPLIFY = FALSE)
            W2 <- do.call(cbind, W2)
            split.Time <- split(c(t(st)), rep(strt, each = con$GKk))
            W2s <- mapply(function (k, t) splineDesign(k, t, ord = con$ord), kn, split.Time, SIMPLIFY = FALSE)
            W2s <- mapply(function (w2s, ind) {
                out <- matrix(0, nT * con$GKk, ncol(w2s))
                out[strt == ind, ] <- w2s
                out
            }, W2s, levels(strt), SIMPLIFY = FALSE)
            W2s <- do.call(cbind, W2s)
            #W2 <- splineDesign(kn, Time, ord = con$ord)
            #W2s <- splineDesign(kn, c(t(st)), ord = con$ord)
            y <- c(y, list(strata = strt))
            x <- c(x, list(W2 = W2, W2s = W2s))
        }
    }
    # extra design matrices for 'method = "piecewise-PH-GH"'
    if (method == "piecewise-PH-GH") {
        wk <- gaussKronrod(con$GKk)$wk
        sk <- gaussKronrod(con$GKk)$sk
        nk <- length(sk)
        if (is.null(con$knots) || !is.numeric(con$knots)) {
            Q <- con$lng.in.kn + 1
            qs <- quantile(Time, seq(0, 1, len = Q + 1), names = FALSE)[-c(1, Q + 1)] + 1e-06
            con$knots <- qs
            qs <- c(0, qs, max(Time) + 1)
        } else {
            qs <- c(0, sort(con$knots), max(Time) + 1)
            Q <- length(qs) - 1
        }
        ind <- findInterval(Time, qs, rightmost.closed = TRUE)
        D <- matrix(0, length(ind), Q)
        D[cbind(seq_along(ind), ind)] <- 1
        D <- D * d
        Tiq <- outer(Time, qs, pmin)
        Lo <- Tiq[, 1:Q]
        Up <- Tiq[, 2:(Q+1)]
        T <- Up - Lo
        P <- T / 2
        P[P < con$tol3] <- as.numeric(NA)
        P1 <- (Up + Lo) / 2
        st <- matrix(0, nY, nk*Q)
        skQ <- rep(sk, Q)
        for (i in seq_len(nY)) {
            st[i, ] <- rep(P[i, ], each = nk) * skQ + rep(P1[i, ], each = nk)
        }
        data.id2 <- data.id[rep(seq_len(nY), each = nk*Q), ]
        data.id2[timeVar] <- pmax(c(t(st)) - lag, 0)
        mf <- model.frame(lmeObject$terms, data = data.id2)
        Xs <- model.matrix(formYx, mf)
        Zs <- model.matrix(formYz, mf)
        y <- c(y, list(ind.D = ind))
        id.GK <- rep(1:nY, rowSums(!is.na(st)))
        P <- c(t(P))
        x <- c(x, list(Xs = Xs, Zs = Zs, P = P[!is.na(P)], st = st[!is.na(st)], wk = wk, id.GK = id.GK, Q = Q))
    }
    # extra design matrices for 'method = "Cox-PH-GH"' with event times prior to observed time for the ith subject
    if (method == "Cox-PH-GH") {
        unqT <- sort(unique(Time[d == 1]))
        times <- lapply(Time, function (t) unqT[t >= unqT])
        ind.len <- sapply(times, length)
        indT <- rep(1:nrow(data.id), ind.len)
        data.id2 <- data.id[indT, ]
        data.id2[timeVar] <- pmax(unlist(times, use.names = FALSE) - lag, 0)
        mf <- model.frame(lmeObject$terms, data = data.id2)
        Xtime2 <- model.matrix(formYx, mf)
        Ztime2 <- model.matrix(formYz, mf)
        x <- c(x, list(Xtime2 = Xtime2, Ztime2 = Ztime2, indT = indT))
    }
    # initial values
    VC <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", lmeObject$sigma^2)[[1]]
    con$inv.chol.VC <- solve(chol(solve(VC)))
    con$det.inv.chol.VC <- det(con$inv.chol.VC)
    if (all(VC[upper.tri(VC)] == 0))
        VC <- diag(VC)
    init.surv <- if (method == "piecewise-PH-GH") {
        initial.surv(T, D, if (is.null(W)) as.matrix(long) else cbind(W, long), method, con, 
            extra = list(Time = Time, d = d, y = y$y, id = id, times = data[[timeVar]]))
    } else {
        initial.surv(log(Time), d, if (is.null(W)) as.matrix(long) else cbind(W, long), method, con, 
            extra = list(y = y$y, id = id, times = data[[timeVar]], W2 = x$W2))
    }
    if (method == "Cox-PH-GH" && length(init.surv$lambda0) < length(unqT))
        init.surv$lambda0 <- basehaz(survObject)$hazard
    initial.values <- c(list(betas = fixef(lmeObject), sigma = lmeObject$sigma, D = VC), init.surv)
    if (!is.null(init)) {
        nams1 <- names(init)
        nams2 <- names(initial.values)
        if (!is.list(init) || length(noNms <- nams1[!nams1 %in% nams2])) {
            warning("unknown names in 'init': ", paste(noNms, collapse = ", "))
        } else {
            initial.values[nams1] <- init
        }
    }
    # joint model fit
    out <- switch(method,
        "Cox-PH-GH" = phGH.fit(x, y, id, initial.values, con),
        "weibull-AFT-GH" = weibullAFTGH.fit(x, y, id, initial.values, con),
        "weibull-PH-GH" = weibullPHGH.fit(x, y, id, initial.values, con),
        "piecewise-PH-GH" = piecewisePHGH.fit(x, y, id, initial.values, con),
        "spline-PH-GH" = splinePHGH.fit(x, y, id, initial.values, con),
        "ch-Laplace" = chLaplace.fit(x, y, id, initial.values, b, con))
    # check for problems with the Hessian at convergence
    H <- out$Hessian
    if (any(is.na(H) | !is.finite(H))) {
        warning("infinite or missing values in Hessian at convergence.\n")
    } else {
        ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        if (!all(ev >= -1e-06 * abs(ev[1]))) 
            warning("Hessian matrix at convergence is not positive definite.\n")
    }
    out$x <- x
    out$y <- y
    out$times <- data[[timeVar]]
    out$data <- data
    out$data.id <- data.id
    out$method <- method
    out$termsY <- lmeObject$terms
    out$termsT <- survObject$terms
    out$formYx <- formYx
    out$formYz <- formYz
    out$formT <- formT
    out$timeVar <- timeVar
    out$control <- con
    out$call <- cl
    class(out) <- "jointModel"
    out
}


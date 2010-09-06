splinePHGH.fit <-
function (x, y, id, initial.values, parameterization, derivForm, control) {
    # response vectors
    logT <- as.vector(y$logT)
    Time <- exp(logT)
    d <- as.vector(y$d)
    strata <- y$strata
    y <- as.vector(y$y)
    # design matrices
    X <- x$X
    Xtime <- x$Xtime
    Xs <- x$Xs
    Xtime.deriv <- x$Xtime.deriv
    Xs.deriv <- x$Xs.deriv    
    Z <- x$Z
    Ztime <- x$Ztime
    Zs <- x$Zs
    Ztime.deriv <- x$Ztime.deriv
    Zs.deriv <- x$Zs.deriv
    W1 <- x$W
    W2 <- x$W2
    W2s <- x$W2s
    WW <- if (is.null(W1)) W2 else cbind(W2, W1)
    X <- dropAttr(X); Z <- dropAttr(Z); WW <- dropAttr(WW)
    if (parameterization == "value") {
        Xtime <- dropAttr(Xtime); Ztime <- dropAttr(Ztime); Xs <- dropAttr(Xs); Zs <- dropAttr(Zs)
    } else if (parameterization == "slope") {
        Xtime.deriv <- dropAttr(Xtime.deriv); Ztime.deriv <- dropAttr(Ztime.deriv)
        Xs.deriv <- dropAttr(Xs.deriv); Zs.deriv <- dropAttr(Zs.deriv)
    } else {
        Xtime <- dropAttr(Xtime); Ztime <- dropAttr(Ztime); Xs <- dropAttr(Xs); Zs <- dropAttr(Zs)
        Xtime.deriv <- dropAttr(Xtime.deriv); Ztime.deriv <- dropAttr(Ztime.deriv)
        Xs.deriv <- dropAttr(Xs.deriv); Zs.deriv <- dropAttr(Zs.deriv)
    }
    indFixed <- derivForm$indFixed
    indRandom <- derivForm$indRandom
    # sample size settings
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(WW)
    nk <- ncol(W2)
    n <- length(logT)
    N <- length(y)
    ni <- as.vector(tapply(id, id, length))
    # crossproducts and others
    XtX <- crossprod(X)
    ZtZ <- lapply(split(Z, id), function (x) crossprod(matrix(x, ncol = ncz)))
    names(ZtZ) <- NULL
    ZtZ <- matrix(unlist(ZtZ), n, ncz * ncz, TRUE)
    outer.Ztime <- lapply(1:n, function (x) Ztime[x, ] %o% Ztime[x, ])
    # Gauss-Hermite quadrature rule components
    GH <- gauher(control$GHk)
    b <- as.matrix(expand.grid(lapply(1:ncz, function (k, u) u$x, u = GH)))
    k <- nrow(b)
    wGH <- as.matrix(expand.grid(lapply(1:ncz, function (k, u) u$w, u = GH)))
    wGH <- 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * b)) * control$det.inv.chol.VC
    b <- sqrt(2) * t(control$inv.chol.VC %*% t(b)); dimnames(b) <- NULL
    b2 <- if (ncz == 1) b * b else t(apply(b, 1, function (x) x %o% x))
    Ztb <- Z %*% t(b)    
    if (parameterization %in% c("value", "both")) {
        Ztime.b <- Ztime %*% t(b)
        Zsb <- Zs %*% t(b)
    }
    if (parameterization %in% c("slope", "both")) {
        if (length(indRandom) > 1 || indRandom) {
            Ztime.b.deriv <- Ztime.deriv %*% t(b[, indRandom, drop = FALSE])
            Zsb.deriv <- Zs.deriv %*% t(b[, indRandom, drop = FALSE])
        } else {
            Ztime.b.deriv <- matrix(0, nrow(Ztime.deriv), k)
            Zsb.deriv <- matrix(0, nrow(Zs.deriv), k)
        }
    }
    # Gauss-Kronrod rule
    wk <- rep(x$wk, length(logT))
    P <- as.vector(x$P)
    id.GK <- rep(seq_along(logT), each = control$GKk)
    # initial values
    betas <- as.vector(initial.values$betas)
    sigma <- initial.values$sigma
    gammas <- if (!is.null(W1)) as.vector(initial.values$gammas) else NULL
    gammas.bs <- as.vector(initial.values$gammas.bs)
    alpha <- as.vector(initial.values$alpha)
    Dalpha <- as.vector(initial.values$Dalpha)
    D <- initial.values$D
    diag.D <- !is.matrix(D)
    if (!diag.D) dimnames(D) <- NULL else names(D) <- NULL
    # fix environments for functions
    environment(opt.survSplinePH) <- environment(gr.survSplinePH) <- environment()
    environment(opt.longSplinePH) <- environment(gr.longSplinePH) <- environment(H.longSplinePH) <- environment()
    environment(LogLik.splineGH) <- environment(Score.splineGH) <- environment()
    old <- options(warn = (-1))
    on.exit(options(old))
    # EM iterations
    iter <- control$iter.EM
    Y.mat <- matrix(0, iter + 1, ncx + 1)
    T.mat <- matrix(0, iter + 1, switch(parameterization, "value" = , "slope" = ncww + 1, "both" = ncww + 2))
    B.mat <- if (diag.D) matrix(0, iter + 1, ncz) else matrix(0, iter + 1, ncz * ncz)
    lgLik <- numeric(iter + 1)
    conv <- FALSE
    for (it in 1:iter) {
        # save parameter values in matrix
        Y.mat[it, ] <- c(betas, sigma)        
        T.mat[it, ] <- switch(parameterization, "value" = c(gammas, alpha, gammas.bs),
            "slope" = c(gammas, Dalpha, gammas.bs), "both" = c(gammas, alpha, Dalpha, gammas.bs))
        B.mat[it, ] <- D
        
        # linear predictors
        eta.yx <- as.vector(X %*% betas)
        eta.tw1 <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, n)
        eta.tw2 <- as.vector(W2 %*% gammas.bs)
        eta.ws <- as.vector(W2s %*% gammas.bs)
        if (parameterization %in% c("value", "both")) {
            Y <- as.vector(Xtime %*% betas) + Ztime.b
            Ys <- as.vector(Xs %*% betas) + Zsb
            eta.t <- eta.tw2 + eta.tw1 + alpha * Y
            eta.s <- alpha * Ys
        }
        if (parameterization %in% c("slope", "both")) {
            Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + Ztime.b.deriv
            Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
            eta.t <- if (parameterization == "both") eta.t + Dalpha * Y.deriv else eta.tw2 + eta.tw1 + Dalpha * Y.deriv
            eta.s <- if (parameterization == "both") eta.s + Dalpha * Ys.deriv else Dalpha * Ys.deriv
        }

        
        # E-step
        mu.y <- eta.yx + Ztb
        logNorm <- dnorm(y, mu.y, sigma, TRUE)
        log.p.yb <- rowsum(logNorm, id, reorder = FALSE); dimnames(log.p.yb) <- NULL
        log.hazard <- eta.t        
        log.survival <- - exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + eta.s), id.GK, reorder = FALSE)
        dimnames(log.survival) <- NULL
        log.p.tb <- d * log.hazard + log.survival
        log.p.b <- if (ncz == 1) {
            dnorm(b, sd = sqrt(D), log = TRUE)
        } else {
            if (diag.D) {
                rowSums(dnorm(b, sd = rep(sqrt(D), each = k), log = TRUE))
            } else {
                dmvnorm(b, rep(0, ncz), D, TRUE)
            }
        }
        p.ytb <- exp((log.p.yb + log.p.tb) + rep(log.p.b, each = n))
        p.yt <- c(p.ytb %*% wGH)
        p.byt <- p.ytb / p.yt
        post.b <- p.byt %*% (b * wGH)
        post.vb <- if (ncz == 1) {
            c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
        } else {
            (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
        }
        
        # compute log-likelihood
        log.p.yt <- log(p.yt)
        lgLik[it] <- sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
        
        # print results if verbose
        if (control$verbose) {
            cat("\n\niter:", it, "\n")
            cat("log-likelihood:", lgLik[it], "\n")
            cat("betas:", round(betas, 4), "\n")
            cat("sigma:", round(sigma, 4), "\n")
            cat("gammas:", round(c(gammas.bs, gammas), 4), "\n")
            if (parameterization %in% c("value", "both"))
                cat("alpha:", round(alpha, 4), "\n")
            if (parameterization %in% c("slope", "both"))
                cat("Dalpha:", round(Dalpha, 4), "\n")
            cat("D:", if (!diag.D) round(D[lower.tri(D, TRUE)], 4) else round(D, 4), "\n")
        }
        
        # check convergence
        if (it > 5 && lgLik[it] > lgLik[it - 1]) {
            thets1 <- c(Y.mat[it - 1, ], T.mat[it - 1, ], B.mat[it - 1, ])
            thets2 <- c(Y.mat[it, ], T.mat[it, ], B.mat[it, ])
            check1 <- max(abs(thets2 - thets1) / (abs(thets1) + control$tol1)) < control$tol2
            check2 <- (lgLik[it] - lgLik[it - 1]) < control$tol3 * (abs(lgLik[it - 1]) + control$tol3)
            if (check1 || check2) {
                conv <- TRUE
                if (control$verbose)
                    cat("\n\nconverged!\ncalculating Hessian...\n")
                break
            }
        }
        if (iter == 0) break
        
        # M-step
        Zb <- rowSums(Z * post.b[id, ], na.rm = TRUE)
        mu <- y - eta.yx
        tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)                
        sigman <- sqrt(c(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / N)
        Dn <- matrix(colMeans(p.byt %*% (b2 * wGH), na.rm = TRUE), ncz, ncz)
        Dn <- if (diag.D) diag(Dn) else 0.5 * (Dn + t(Dn))
        Hbetas <- nearPD(H.longSplinePH(betas))
        scbetas <- gr.longSplinePH(betas)
        betasn <- betas - c(solve(Hbetas, scbetas))
        list.thetas <- list(gammas = gammas, alpha = alpha, Dalpha = Dalpha, gammas.bs = gammas.bs)
        list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
        thetas <- unlist(as.relistable(list.thetas))
        optz.surv <- optim(thetas, opt.survSplinePH, gr.survSplinePH, method = "BFGS", 
            control = list(maxit = if (it < 5) 20 else 4, 
                parscale = if (it < 5) rep(0.01, length(thetas)) else rep(0.1, length(thetas))))

        # update parameter values
        betas <- betasn
        sigma <- sigman
        D <- Dn
        thetasn <- relist(optz.surv$par, skeleton = list.thetas)
        gammas <- if (!is.null(W1)) thetasn$gammas else NULL
        alpha <- thetasn$alpha
        Dalpha <- thetasn$Dalpha
        gammas.bs <- thetasn$gammas.bs
    }
    list.thetas <- list(betas = betas, log.sigma = log(sigma), gammas = gammas, alpha = alpha, Dalpha = Dalpha, 
        gammas.bs = gammas.bs, D = if (diag.D) log(D) else chol.transf(D))
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    lgLik <- - LogLik.splineGH(thetas)
    # if not converged, start quasi-Newton iterations
    if (!conv && !control$only.EM) {
        if (is.null(control$parscale))
            control$parscale <- rep(0.01, length(thetas))
        if (control$verbose)
            cat("\n\nquasi-Newton iterations start.\n\n")
        out <- if (control$optimizer == "optim") {
            optim(thetas, LogLik.splineGH, Score.splineGH, method = "BFGS",
                control = list(maxit = control$iter.qN, parscale = control$parscale, 
                trace = 10 * control$verbose))
        } else {
            nlminb(thetas, LogLik.splineGH, Score.splineGH, scale = control$parscale, 
                control = list(iter.max = control$iter.qN, trace = 1 * control$verbose))
        }
        if ((conv <- out$convergence) == 0 || - out[[2]] > lgLik) {
            lgLik <- - out[[2]]            
            thetas <- relist(out$par, skeleton = list.thetas)
            betas <- thetas$betas
            sigma <- exp(thetas$log.sigma)
            gammas <- if (!is.null(W1)) thetas$gammas else NULL
            gammas.bs <- thetas$gammas.bs
            alpha <- thetas$alpha
            Dalpha <- thetas$Dalpha
            D <- thetas$D
            D <- if (diag.D) exp(D) else chol.transf(D)
            it <- it + if (control$optimizer == "optim") out$counts[1] else out$iterations
            # compute posterior moments for thetas after quasi-Newton
            eta.yx <- as.vector(X %*% betas)
            eta.tw1 <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, n)
            eta.tw2 <- as.vector(W2 %*% gammas.bs)
            exp.eta.tw <- exp(eta.tw1)
            eta.ws <- as.vector(W2s %*% gammas.bs)            
            if (parameterization %in% c("value", "both")) {
                Y <- as.vector(Xtime %*% betas) + Ztime.b
                Ys <- as.vector(Xs %*% betas) + Zsb
                eta.t <- eta.tw2 + eta.tw1 + alpha * Y
                eta.s <- alpha * Ys
            }
            if (parameterization %in% c("slope", "both")) {
                Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + Ztime.b.deriv
                Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
                eta.t <- if (parameterization == "both") eta.t + Dalpha * Y.deriv else eta.tw2 + eta.tw1 + Dalpha * Y.deriv
                eta.s <- if (parameterization == "both") eta.s + Dalpha * Ys.deriv else Dalpha * Ys.deriv
            }
            mu.y <- eta.yx + Ztb
            logNorm <- dnorm(y, mu.y, sigma, TRUE)
            log.p.yb <- rowsum(logNorm, id)
            log.hazard <- eta.t
            log.survival <- - exp.eta.tw * P * rowsum(wk * exp(eta.ws + eta.s), id.GK, reorder = FALSE)
            dimnames(log.survival) <- NULL
            log.p.tb <- d * log.hazard + log.survival            
            log.p.b <- if (ncz == 1) {
                dnorm(b, sd = sqrt(D), log = TRUE)
            } else {
                if (diag.D) {
                    rowSums(dnorm(b, sd = rep(sqrt(D), each = k), log = TRUE))
                } else {
                    dmvnorm(b, rep(0, ncz), D, TRUE)
                }
            }
            p.ytb <- exp((log.p.yb + log.p.tb) + rep(log.p.b, each = n))
            dimnames(p.ytb) <- NULL
            p.yt <- c(p.ytb %*% wGH)
            p.byt <- p.ytb / p.yt
            post.b <- p.byt %*% (b * wGH)
            post.vb <- if (ncz == 1) {
                c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
            } else {
                (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
            }
            Zb <- if (ncz == 1) post.b[id] else rowSums(Z * post.b[id, ], na.rm = TRUE)
            if (control$verbose)
                cat("\n\nconverged!\ncalculating Hessian...\n")
        }
    }
    # calculate Hessian matrix
    Hessian <- if (control$numeriDeriv == "fd") {
        fd.vec(unlist(thetas), Score.splineGH, eps = control$eps.Hes)
    } else { 
        cd.vec(unlist(thetas), Score.splineGH, eps = control$eps.Hes)
    }
    names(betas) <- names(initial.values$betas)
    if (!diag.D) dimnames(D) <- dimnames(initial.values$D) else names(D) <- names(initial.values$D)
    names(gammas) <- colnames(W1)
    names(gammas.bs) <- if (length(levels(strata)) == 1) paste("bs", 1:nk, sep = "") else {
        len.kn <- sapply(control$knots, length) - control$ord
        paste("bs", sapply(len.kn, seq_len), "(", rep(levels(strata), len.kn), ")", sep = "")
    }
    gg <- switch(parameterization, "value" = "alpha", "slope" = "alphaD", "both" = c("alpha", "alphaD"))
    nams <- c(paste("Y.", c(names(betas), "sigma"), sep = ""), 
        paste("T.", c(names(gammas), gg, names(gammas.bs)), sep = ""),
        paste("B.", if (!diag.D) paste("D", seq(1, ncz * (ncz + 1) / 2), sep = "") else names(D), sep = ""))
    dimnames(Hessian) <- list(nams, nams)
    colnames(post.b) <- colnames(x$Z)
    list(coefficients = list(betas = betas, sigma = sigma, gammas = gammas, alpha = alpha, Dalpha = Dalpha, gammas.bs = gammas.bs,
        D = as.matrix(D)), Hessian = Hessian, logLik = lgLik, EB = list(post.b = post.b, post.vb = post.vb, 
        Zb = if (iter == 0) rowSums(Z * post.b[id, ], na.rm = TRUE) else Zb,
        Ztimeb = if (parameterization %in% c("value", "both")) rowSums(Ztime * post.b) else NULL,
        Ztimeb.deriv = if (parameterization %in% c("slope", "both")) {
            if (indRandom) rowSums(Ztime.deriv * post.b[, indRandom, drop = FALSE]) else rep(0, nrow(Ztime.deriv))
        } else NULL), iters = it, convergence = conv, n = n, N = N, ni = ni, d = d, id = id)
}


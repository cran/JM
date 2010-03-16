weibullAFTGH.fit <-
function (x, y, id, initial.values, control) {
    # response vectors
    logT <- as.vector(y$logT)
    d <- as.vector(y$d)
    y <- as.vector(y$y)
    # design matrices
    X <- x$X
    Xtime <- x$Xtime
    Xs <- x$Xs
    Z <- x$Z
    Ztime <- x$Ztime
    Zs <- x$Zs
    W1 <- x$W
    WW <- if (is.null(W1)) as.matrix(rep(1, length(logT))) else cbind(1, W1)
    dimnames(X) <- dimnames(Xtime) <- dimnames(Xs) <- dimnames(Z) <- dimnames(Ztime) <- dimnames(Zs) <- dimnames(WW) <- NULL
    attr(X, "assign") <- attr(X, "contrasts") <- attr(Xtime, "assign") <- attr(Xtime, "contrasts") <- NULL
    attr(Xs, "assign") <- attr(Xs, "contrasts") <- attr(Zs, "assign") <- attr(Zs, "contrasts") <- attr(Z, "assign") <- attr(Ztime, "assign") <- NULL
    # sample size settings
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(WW)
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
    Ztime.b <- Ztime %*% t(b)
    Zsb <- Zs %*% t(b)
    # Gauss-Kronrod rule
    st <- x$st
    log.st <- log(st)
    wk <- rep(x$wk, length(logT))
    P <- as.vector(x$P)
    id.GK <- rep(seq_along(logT), each = control$GKk)
    # initial values
    betas <- as.vector(initial.values$betas)
    sigma <- initial.values$sigma
    gammas <- as.vector(initial.values$gammas)
    alpha <- as.vector(initial.values$alpha)
    sigma.t <- initial.values$sigma.t
    D <- initial.values$D
    diag.D <- !is.matrix(D)
    if (!diag.D) dimnames(D) <- NULL else names(D) <- NULL
    # fix environments for functions
    environment(opt.survAFTWB) <- environment(gr.survAFTWB) <- environment()
    environment(opt.longAFTWB) <- environment(gr.longAFTWB) <- environment()
    environment(LogLik.weibullAFTGH) <- environment(Score.weibullAFTGH) <- environment()
    old <- options(warn = (-1))
    on.exit(options(old))
    # EM iterations
    iter <- control$iter.EM
    Y.mat <- matrix(0, iter + 1, ncx + 1)
    T.mat <- matrix(0, iter + 1, ncww + 2)
    B.mat <- if (diag.D) matrix(0, iter + 1, ncz) else matrix(0, iter + 1, ncz * ncz)
    lgLik <- numeric(iter + 1)
    conv <- FALSE
    for (it in 1:iter) {
        # save parameter values in matrix
        Y.mat[it, ] <- c(betas, sigma)
        T.mat[it, ] <- c(gammas, alpha, sigma.t)
        B.mat[it,] <- D
        
        # linear predictors
        eta.yx <- as.vector(X %*% betas)
        eta.yxT <- as.vector(Xtime %*% betas)
        eta.tw <- as.vector(WW %*% gammas)
        Y <- eta.yxT + Ztime.b
        Ys <- as.vector(Xs %*% betas) + Zsb
        eta.t <- eta.tw + alpha * Y
        eta.s <- alpha * Ys
        
        # E-step
        mu.y <- eta.yx + Ztb
        logNorm <- dnorm(y, mu.y, sigma, TRUE)
        log.p.yb <- rowsum(logNorm, id, reorder = FALSE); dimnames(log.p.yb) <- NULL
        Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
        log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
        log.survival <- - Vi^sigma.t
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
            cat("gammas:", -round(gammas, 4), "\n")
            cat("alpha:", -round(alpha, 4), "\n")
            cat("sigma.t:", round(sigma.t, 4), "\n")
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
        Hbetas <- nearPD(fd.vec(betas, gr.longAFTWB))
        scbetas <- gr.longAFTWB(betas)
        betasn <- betas - c(solve(Hbetas, scbetas))
        thetas <- c(gammas, alpha, log(sigma.t))
        optz.surv <- optim(thetas, opt.survAFTWB, gr.survAFTWB, method = "BFGS", 
            control = list(maxit = if (it < 5) 20 else 5, 
                parscale = if (it < 10) rep(0.01, length(thetas)) else rep(0.1, length(thetas))))
        thetasn <- optz.surv$par

        # update parameter values
        betas <- betasn
        sigma <- sigman
        D <- Dn
        gammas <- thetasn[1:ncww]
        alpha <- thetasn[ncww + 1]
        sigma.t <- exp(thetasn[ncww + 2])
    }
    list.thetas <- list(betas = betas, log.sigma = log(sigma), gammas = gammas, alpha = alpha, 
        log.sigma.t = log(sigma.t), D = if (diag.D) log(D) else chol.transf(D))
    thetas <- unlist(as.relistable(list.thetas))
    lgLik <- - LogLik.weibullAFTGH(thetas)
    
    # if not converged, start quasi-Newton iterations
    if (!conv && !control$only.EM) {
        if (is.null(control$parscale))
            control$parscale <- rep(0.01, length(thetas))
        if (control$verbose)
            cat("\n\nquasi-Newton iterations start.\n\n")
        out <- if (control$optimizer == "optim") {
            optim(thetas, LogLik.weibullAFTGH, Score.weibullAFTGH, method = "BFGS",
                control = list(maxit = control$iter.qN, parscale = control$parscale, 
                trace = 10 * control$verbose))
        } else {
            nlminb(thetas, LogLik.weibullAFTGH, Score.weibullAFTGH, scale = control$parscale, 
                control = list(iter.max = control$iter.qN, trace = 1 * control$verbose))
        }
        if ((conv <- out$convergence) == 0 || - out[[2]] > lgLik) {
            lgLik <- - out[[2]]
            thetas <- relist(out$par, skeleton = list.thetas)
            betas <- thetas$betas
            sigma <- exp(thetas$log.sigma)
            gammas <- thetas$gammas
            alpha <- thetas$alpha
            sigma.t <- exp(thetas$log.sigma.t)
            D <- thetas$D
            D <- if (diag.D) exp(D) else chol.transf(D)
            it <- it + if (control$optimizer == "optim") out$counts[1] else out$iterations
            # compute posterior moments for thetas after quasi-Newton
            eta.yx <- as.vector(X %*% betas)
            eta.yxT <- as.vector(Xtime %*% betas)
            eta.tw <- as.vector(WW %*% gammas)
            exp.eta.tw <- exp(eta.tw)
            Y <- eta.yxT + Ztime.b
            Ys <- as.vector(Xs %*% betas) + Zsb
            eta.t <- eta.tw + alpha * Y
            eta.s <- alpha * Ys
            mu.y <- eta.yx + Ztb
            logNorm <- dnorm(y, mu.y, sigma, TRUE)
            log.p.yb <- rowsum(logNorm, id)
            Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
            log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
            log.survival <- - Vi^sigma.t
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
        }
    }
    # calculate Hessian matrix
    Hessian <- if (control$numeriDeriv == "fd") {
        fd.vec(unlist(thetas), Score.weibullAFTGH, eps = control$eps.Hes)
    } else { 
        cd.vec(unlist(thetas), Score.weibullAFTGH, eps = control$eps.Hes)
    }
    names(betas) <- names(initial.values$betas)
    if (!diag.D) dimnames(D) <- dimnames(initial.values$D) else names(D) <- names(initial.values$D)
    names(gammas) <- c("(Intercept)", colnames(W1))
    nams <- c(paste("Y.", c(names(betas), "sigma"), sep = ""), paste("T.", c(names(gammas), "alpha", "sigma.t"), sep = ""),
        paste("B.", if (!diag.D) paste("D", seq(1, ncz * (ncz + 1) / 2), sep = "") else names(D), sep = ""))
    dimnames(Hessian) <- list(nams, nams)
    colnames(post.b) <- colnames(x$Z)
    list(coefficients = list(betas = betas, sigma = sigma, gammas = gammas, alpha = alpha, sigma.t = sigma.t, 
        D = as.matrix(D)), Hessian = Hessian, logLik = lgLik, EB = list(post.b = post.b, post.vb = post.vb, 
        Zb = if (iter == 0) rowSums(Z * post.b[id, ], na.rm = TRUE) else Zb, 
        Ztimeb = rowSums(Ztime * post.b)), iters = it, convergence = conv, n = n, N = N, ni = ni, d = d, id = id)
}


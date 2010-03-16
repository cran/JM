phGH.fit <-
function (x, y, id, initial.values, control) {
    # response vectors
    logT <- as.vector(y$logT)
    d <- as.vector(y$d)
    y <- as.vector(y$y)
    # design matrices
    X <- x$X
    Xtime <- x$Xtime
    Xtime2 <- x$Xtime2
    Z <- x$Z
    Ztime <- x$Ztime
    Ztime2 <- x$Ztime2
    WW <- x$W
    dimnames(X) <- dimnames(Xtime) <- dimnames(Xtime2) <- NULL
    dimnames(Z) <- dimnames(Ztime) <- dimnames(Ztime2) <- dimnames(WW) <- NULL
    attr(X, "assign") <- attr(X, "contrasts") <- NULL
    attr(Xtime, "assign") <- attr(Xtime, "contrasts") <- NULL
    attr(Xtime2, "assign") <- attr(Xtime2, "contrasts") <- NULL    
    attr(Z, "assign") <- attr(Ztime, "assign") <- NULL
    # sample size settings
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- if (is.null(WW)) 0 else ncol(WW)
    n <- length(logT)
    N <- length(y)
    ni <- as.vector(tapply(id, id, length))
    # crossproducts and others
    XtX <- crossprod(X)
    ZtZ <- lapply(split(Z, id), function (x) crossprod(matrix(x, ncol = ncz)))
    names(ZtZ) <- NULL
    ZtZ <- matrix(unlist(ZtZ), n, ncz * ncz, TRUE)
    outer.Ztime <- lapply(1:n, function (x) Ztime[x, ] %o% Ztime[x, ])
    Time <- exp(logT)
    unqT <- sort(unique(Time[d == 1]))
    indT <- x$indT
    unq.indT <- unique(indT)
    nk <- as.vector(sapply(split(indT, indT), length))
    ind.L1 <- unlist(lapply(nk, seq, from = 1))
    ind.L2 <- colSums(d * outer(Time, unqT, "=="))
    ind.T0 <- match(Time, unqT)
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
    Ztime2.b <- Ztime2 %*% t(b)
    # initial values
    betas <- as.vector(initial.values$betas)
    sigma <- initial.values$sigma
    gammas <- as.vector(initial.values$gammas)
    alpha <- as.vector(initial.values$alpha)
    lambda0 <- initial.values$lambda0
    D <- initial.values$D
    diag.D <- !is.matrix(D)
    if (!diag.D) dimnames(D) <- NULL else names(D) <- NULL
    # fix environments for functions
    environment(opt.survPH) <- environment(gr.survPH) <- environment(gr.longPH) <- environment()
    environment(LogLik.phGH) <- environment(Score.phGH) <- environment()
    old <- options(warn = (-1))
    on.exit(options(old))
    # EM iterations
    iter <- control$iter.EM
    Y.mat <- matrix(0, iter + 1, ncx + 1)
    T.mat <- matrix(0, iter + 1, ncww + 1)
    B.mat <- if (diag.D) matrix(0, iter + 1, ncz) else matrix(0, iter + 1, ncz * ncz)
    lgLik <- numeric(iter + 1)
    conv <- FALSE
    for (it in 1:iter) {
        # save parameter values in matrix
        Y.mat[it, ] <- c(betas, sigma)
        T.mat[it, ] <- c(gammas, alpha)
        B.mat[it,] <- D
        
        # linear predictors
        eta.yx <- as.vector(X %*% betas)
        eta.yxT <- as.vector(Xtime %*% betas)
        eta.yxT2 <- as.vector(Xtime2 %*% betas)
        Y <- eta.yxT + Ztime.b
        Y2 <- eta.yxT2 + Ztime2.b
        eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else rep(0, n)
        eta.t <- eta.tw + alpha * Y
        eta.s <- alpha * Y2
        exp.eta.s <- exp(eta.s)
         
        # E-step
        mu.y <- eta.yx + Ztb
        logNorm <- dnorm(y, mu.y, sigma, TRUE)
        log.p.yb <- rowsum(logNorm, id); dimnames(log.p.yb) <- NULL
        log.lambda0T <- log(lambda0[ind.T0])
        log.lambda0T[is.na(log.lambda0T)] <- 0
        log.hazard <- log.lambda0T + eta.t
        S <- matrix(0, n, k)
        S[unq.indT, ] <- rowsum(lambda0[ind.L1] * exp.eta.s, indT, reorder = FALSE); dimnames(S) <- NULL
        log.survival <- - exp(eta.tw) * S
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
            if (!is.null(WW))
                cat("gammas:", round(gammas, 4), "\n")
            cat("alpha:", round(alpha, 4), "\n")
            cat("D:", if (!diag.D) round(D[lower.tri(D, TRUE)], 4) else round(D, 4), "\n")
        }
        
        # check convergence
        if (it > 5) {
            if (lgLik[it] > lgLik[it - 1]) {
                thets1 <- c(Y.mat[it - 1, ], T.mat[it - 1, ], B.mat[it - 1, ])
                thets2 <- c(Y.mat[it, ], T.mat[it, ], B.mat[it, ])
                check1 <- max(abs(thets2 - thets1) / (abs(thets1) + control$tol1)) < control$tol2
                check2 <- (lgLik[it] - lgLik[it - 1]) < control$tol3 * (abs(lgLik[it - 1]) + control$tol3)
                if (check1 || check2) {
                    conv <- TRUE
                    if (control$verbose)
                        cat("\n\nconverged!\n")
                    break
                }
            } else {
                lambda0 <- lambda0.old
                log.lambda0T <- log(lambda0[ind.T0])
                log.lambda0T[is.na(log.lambda0T)] <- 0
                log.hazard <- log.lambda0T + eta.t
                S <- matrix(0, n, k)
                S[unq.indT, ] <- rowsum(lambda0[ind.L1] * exp.eta.s, indT, reorder = FALSE)
                log.survival <- - exp(eta.tw) * S
                log.p.tb <- d * log.hazard + log.survival
                p.ytb <- exp((log.p.yb + log.p.tb) + rep(log.p.b, each = n))
                p.yt <- c(p.ytb %*% wGH)
                p.byt <- p.ytb / p.yt
                post.b <- p.byt %*% (b * wGH)
                post.vb <- if (ncz == 1) {
                    c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
                } else {
                    (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
                }
                log.p.yt <- log(p.yt)
                lgLik[it] <- sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
                if (control$verbose) {
                    cat("\n\niter:", it, "\n")
                    cat("log-likelihood:", lgLik[it], "\n")
                    cat("betas:", round(betas, 4), "\n")
                    cat("sigma:", round(sigma, 4), "\n")
                    if (!is.null(WW))
                        cat("gammas:", round(gammas, 4), "\n")
                    cat("alpha:", round(alpha, 4), "\n")
                    cat("D:", if (!diag.D) round(D[lower.tri(D, TRUE)], 4) else round(D, 4), "\n")
                }
            }
        }
        if (iter == 0) break
               
        # M-step
        if (it > 2) {
            Zb <- rowSums(Z * post.b[id, ], na.rm = TRUE)
            mu <- y - eta.yx
            tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)
            sigman <- sqrt(c(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / N)
            Dn <- matrix(colMeans(p.byt %*% (b2 * wGH), na.rm = TRUE), ncz, ncz)
            Dn <- if (diag.D) diag(Dn) else 0.5 * (Dn + t(Dn))
            Hbetas <- nearPD(fd.vec(betas, gr.longPH))
            scbetas <- gr.longPH(betas)
            betasn <- betas - c(solve(Hbetas, scbetas))
            thetas <- c(gammas, alpha)
            Hthetas <- nearPD(fd.vec(thetas, gr.survPH))
            scthetas <- gr.survPH(thetas)
            thetasn <- thetas - c(solve(Hthetas, scthetas))
        }
        ee <- c((exp.eta.s * exp(eta.tw[indT]) * p.byt[indT, ]) %*% wGH)
        lambda0n <- ind.L2 / as.vector(tapply(ee, ind.L1, sum, na.rm = TRUE))
        
        # update parameter values
        if (it > 2) {
            betas <- betasn
            sigma <- sigman
            D <- Dn
            gammas <- thetasn[seq_len(ncww)]
            alpha <- thetasn[ncww + 1]
        }
        lambda0.old <- lambda0
        lambda0 <- lambda0n
    }
    thetas <- c(betas, log(sigma), gammas, alpha, if (diag.D) log(D) else chol.transf(D))
    lgLik <- lgLik[it]
    # calculate Hessian matrix
    if (control$verbose) cat("\ncalculating Hessian...\n")
    Hessian <- if (control$numeriDeriv == "fd") {
        fd.vec(thetas, Score.phGH, eps = control$eps.Hes)
    } else { 
        cd.vec(thetas, Score.phGH, eps = control$eps.Hes)
    }
    names(betas) <- names(initial.values$betas)
    if (!diag.D) dimnames(D) <- dimnames(initial.values$D) else names(D) <- names(initial.values$D)
    names(gammas) <- colnames(x$W)
    nams <- c(paste("Y.", c(names(betas), "sigma"), sep = ""), paste("T.", c(names(gammas), "alpha"), sep = ""),
        paste("B.", if (!diag.D) paste("D", seq(1, ncz * (ncz + 1) / 2), sep = "") else names(D), sep = ""))
    dimnames(Hessian) <- list(nams, nams)
    colnames(post.b) <- colnames(x$Z)
    list(coefficients = list(betas = betas, sigma = sigma, gammas = gammas, alpha = alpha, 
        lambda0 = cbind("basehaz" = lambda0, "time" = unqT), D = as.matrix(D)), Hessian = Hessian, logLik = lgLik, 
        EB = list(post.b = post.b, post.vb = post.vb, Zb = if (iter == 0) rowSums(Z * post.b[id, ], na.rm = TRUE) else Zb, 
        Ztimeb = rowSums(Ztime * post.b), Ztime2b = rowSums(Ztime2 * post.b[indT, ])), indexes = list(indT = indT, ind.L1 = ind.L1), 
        iters = it, convergence = conv, n = n, N = N, ni = ni, d = d, id = id)
}


`phGH.fit` <-
function (x, y, id, indT, ind.len, initial.values, control) {
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
    ind.lambda <- unlist(sapply(ind.len, function (x) {
        if (x > 0) seq(1, x) else NULL
    }), use.names = FALSE)
    ind.T0 <- match(Time, unqT)
    ind.T0[d == 0] <- NA
    ind0 <- d == 0
    ind.lenN0 <- ind.len != 0
    indL <- colSums(d * outer(Time, unqT, "=="))
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
    environment(gr.survPH) <- environment(gr.longPH) <- environment()
    environment(HessSurvPH) <- environment(HessLongPH) <- environment()
    environment(LogLik.phGH) <- environment(Score.phGH) <- environment()
    old <- options(warn = (-1))
    on.exit(options(old))
    # EM iterations
    iter <- control$iter.EM
    Y.mat <- matrix(0, iter, ncx + 1)
    T.mat <- matrix(0, iter, ncww + 1)
    B.mat <- if (diag.D) matrix(0, iter, ncz) else matrix(0, iter, ncz * ncz)
    lgLik <- numeric(iter)
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
        eta.t <- if (is.null(WW)) alpha * Y else as.vector(WW %*% gammas) + alpha * Y
        eta.t2 <- if (is.null(WW)) alpha * Y2 else as.vector(WW %*% gammas)[indT] + alpha * Y2
        
        # compute lambda0T
        lambda0T <- lambda0[ind.T0]
        lambda0T[is.na(lambda0T)] <- 0
        
        # E-step
        mu.y <- eta.yx + Ztb
        logNorm <- dnorm(y, mu.y, sigma, TRUE)
        log.p.yb <- rowsum(logNorm, id); dimnames(log.p.yb) <- NULL
        ew <- exp(eta.t2)
        log.p.tb <- d * (log(lambda0T) + eta.t)
        log.p.tb[ind0, ] <- 0
        log.p.tb[ind.lenN0, ] <- log.p.tb[ind.lenN0, ] - rowsum(lambda0[ind.lambda] * ew, indT)
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
                
        # compute log-likelihood and check convergence
        log.p.yt <- log(p.yt)
        lgLik[it] <- sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
        if (it > 5) {
            if (lgLik[it] < lgLik[it - 1]) {
                betas <- Y.mat[it - 1, 1:ncx]
                sigma <- Y.mat[it - 1, ncx + 1]
                gammas <- T.mat[it - 1, 1:ncww]
                alpha <- T.mat[it - 1, ncww + 1]
                D <- B.mat[it - 1,  ]
                if (!diag.D) dim(D) <- c(ncz, ncz)
                lgLik[it] <- lgLik[it - 1]
                eta.t2 <- if (is.null(WW)) alpha * Y2 else as.vector(WW %*% gammas)[indT] + alpha * Y2
                ew <- exp(eta.t2)
                lambda0 <- lambda0.old
                break
            } else {
                thets1 <- c(Y.mat[it - 1, ], T.mat[it - 1, ], B.mat[it - 1, ])
                thets2 <- c(Y.mat[it, ], T.mat[it, ], B.mat[it, ])
                check1 <- max(abs(thets2 - thets1) / (abs(thets1) + control$tol1)) < control$tol2
                check2 <- (lgLik[it] - lgLik[it - 1]) < control$tol3 * (abs(lgLik[it - 1]) + control$tol3)
                if (check1 || check2) {
                    conv <- TRUE
                    if (control$verbose) cat("\n\nconverged!\n")
                    break
                }
            }
        }
        
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
        
        # M-step
        if (it > 1) {
            Zb <- rowSums(Z * post.b[id, ], na.rm = TRUE)
            mu <- y - eta.yx
            tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)
            sigman <- sqrt(c(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / N)
            Dn <- matrix(colMeans(p.byt %*% (b2 * wGH), na.rm = TRUE), ncz, ncz)
            Dn <- if (diag.D) diag(Dn) else 0.5 * (Dn + t(Dn))
            lambda0. <- lambda0[ind.lambda]
            p.byt. <- p.byt[ind.lenN0, ]
            Hbetas <- nearPD(HessLongPH(X, Xtime, Xtime2, ew))
            scbetas <- gr.longPH(X, Xtime, Xtime2, ew)
            betasn <- betas - c(solve(Hbetas, scbetas))
            thetas <- c(gammas, alpha)
            Hthetas <- nearPD(HessSurvPH(WW, Y, Y2, ew))
            scthetas <- gr.survPH(WW, Y, Y2, ew)
            thetasn <- thetas - c(solve(Hthetas, scthetas))
        }
        ee <- c((ew * p.byt[indT, ]) %*% wGH)        
        lambda0n <- indL / as.vector(sapply(split(ee, ind.lambda), sum))        
        
        # update parameter values
        if (it > 1) {
            betas <- betasn
            sigma <- sigman
            D <- Dn
            if (is.null(WW)) alpha <- thetasn else { gammas <- thetasn[1:ncww]; alpha <- thetasn[ncww + 1] }
        }
        lambda0.old <- lambda0
        lambda0 <- lambda0n
    }
    thetas <- c(betas, log(sigma), gammas, alpha, if (diag.D) log(D) else chol.transf(D))
    lgLik <- lgLik[it]
    # calculate Hessian matrix
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
        EB = list(post.b = post.b, post.vb = post.vb, Zb = Zb, Ztimeb = rowSums(Ztime * post.b), 
        Ztime2b = rowSums(Ztime2 * post.b[indT, ])), indexes = list(indT = indT, ind.lambda = ind.lambda, indL = indL), 
        iters = it, convergence = conv, n = n, N = N, ni = ni, d = d, id = id)
}


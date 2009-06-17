MI.fixed.times <-
function (time.points) {
    # indexes for missing data
    if (is.null(time.points))
        time.points <- sort(unique(obs.times))
    ind.miss <- id %in% which(ni < length(time.points))
    id.miss <- id[ind.miss]
    unq.id.miss <- unique(id.miss)
    ni.miss <- length(time.points) - ni
    ni.miss <- ni.miss[ni.miss > 0]
    id2.miss <- rep(seq_along(ni.miss), ni.miss)
    ni <- as.vector(tapply(id.miss, id.miss, length))
    id3.miss <- rep(seq_along(ni), ni)
    id.GK <- rep(seq_len(n) %in% unq.id.miss, each = object$control$GKk)
    # observed data corresponding to patients with
    # one or more missing values
    y.missO <- y[ind.miss]
    logT.missO <- logT[unq.id.miss]
    d.missO <- d[unq.id.miss]
    X.missO <- X[ind.miss, ]
    Z.missO <- Z[ind.miss, ]
    Xtime.missO <- Xtime[unq.id.miss, , drop = FALSE]
    Ztime.missO <- Ztime[unq.id.miss, , drop = FALSE]
    P.missO <- P[unq.id.miss]
    log.st.missO <- log.st[id.GK]
    Xs.missO <- Xs[id.GK, , drop = FALSE]
    Zs.missO <- Zs[id.GK, , drop = FALSE]
    WW.missO <- WW[unq.id.miss, ]
    # observed data corresponding to patients with
    # one or more missing values
    mis.times <- unlist(lapply(split(obs.times, id), function (x) time.points[!time.points %in% x]))
    dataM <- object$data[unq.id.miss, ]
    dataM <- dataM[id2.miss, ]
    dataM[object$timeVar] <- mis.times
    mf <- model.frame(object$termsY, data = dataM)
    X.missM <- model.matrix(object$formYx, mf)
    Z.missM <- model.matrix(object$formYz, mf)
    N.missM <- nrow(X.missM)
    n.missO <- nrow(Ztime.missO)
    # Estimated MLEs
    D <- object$coefficients$D
    diag.D <- ncz != ncol(D)
    thets <- c(object$coefficients$gammas, object$coefficients$alpha)
    thetas <- c(object$coefficients$betas, log(object$coefficients$sigma),
        if (object$method %in% c("weibull-PH-GH", "weibull-AFT-GH", "ph-GH")) thets else {
            thets[2:nk] <- log(diff(thets[1:nk]))
            thets
        }, if (object$method == "weibull-PH-GH" || object$method == "weibull-AFT-GH") log(object$coefficients$sigma.t) else NULL, 
        if (diag.D) log(D) else chol.transf(D))
    V.thetas <- vcov(object)
    EBs <- ranef(object, postVar = TRUE)
    Var <- attr(EBs, "postVar")[unq.id.miss]
    EBs <- proposed.b <- EBs[unq.id.miss, , drop = FALSE]
    # Fitted values for corresponding to Y_i^m
    fitted.valsM <- if (type == "Marginal" || type == "stand-Marginal") {
        as.vector(X.missM %*% object$coefficients$betas)
    } else {
        as.vector(X.missM %*% object$coefficients$betas + rowSums(Z.missM * EBs[id2.miss, , drop = FALSE]))
    }   
    current.b <- b.new <- EBs
    resid.valsM <- matrix(0, N.missM, M)
    environment(posterior.b) <- environment()
    old <- options(warn = (-1))
    on.exit(options(old))
    for (m in 1:M) {
        # Step1: simulate new parameter values from a multivariate normal
        thetas.new <- mvrnorm(1, thetas, V.thetas)
        betas.new <- thetas.new[1:ncx]
        sigma.new <- exp(thetas.new[ncx + 1])
        gammas.new <- thetas.new[seq(ncx + 2, ncx + 1 + ncww)]
        if (object$method == "ch-GH" || object$method == "ch-Laplace")
            gammas.new[1:nk] <- cumsum(c(gammas.new[1], exp(gammas.new[2:nk])))
        alpha.new <- thetas.new[ncx + ncww + 2]
        if (object$method == "weibull-PH-GH" || object$method == "weibull-AFT-GH") {
            sigma.t.new <- exp(thetas.new[ncx + ncww + 3])
            D.new <- thetas.new[seq(ncx + ncww + 4, length(thetas))]
            D.new <- if (diag.D) exp(D.new) else chol.transf(D.new)
        } else {   
            D.new <- thetas.new[seq(ncx + ncww + 3, length(thetas))]
            D.new <- if (diag.D) exp(D.new) else chol.transf(D.new)
        }
        # Step2: Simulate new values for the random effects
        eta.yx <- as.vector(X.missO %*% betas.new)
        eta.yxT <- as.vector(Xtime.missO %*% betas.new)
        eta.tw <- as.vector(WW.missO %*% gammas.new)
        dmvt.current <- dmvt.proposed <- numeric(n.missO)
        for (i in 1:n.missO) {
            proposed.b[i, ] <- rmvt(1, EBs[i, ], Var[[i]], 4)
            tt <- dmvt(rbind(current.b[i, ], proposed.b[i, ]), EBs[i, ], Var[[i]], 4, TRUE)
            dmvt.current[i] <- tt[1]
            dmvt.proposed[i] <- tt[2]
        }
        a <- pmin(exp(posterior.b(proposed.b) + dmvt.current - posterior.b(current.b) - dmvt.proposed), 1)
        ind <- runif(n.missO) <= a
        b.new[ind, ] <- proposed.b[ind, ]
        current.b <- b.new
        # Step3: Simulate new Y_i^m and calculate residuals
        mu <- as.vector(X.missM %*% betas.new + rowSums(Z.missM * b.new[id2.miss, , drop = FALSE]))
        y.new <- rnorm(N.missM, mu, sigma.new)
        resid.valsM[, m] <- y.new - fitted.valsM
    }
    mean.resid.valsM <- rowMeans(resid.valsM)
    if (type == "stand-Subject") {
        var.resid.valsM <- object$coefficients$sigma^2 + apply(resid.valsM, 1, var)
        mean.resid.valsM <- mean.resid.valsM / sqrt(var.resid.valsM)
    }
    if (type == "stand-Marginal") {
        mean.resid.valsM <- unlist(lapply(split(cbind(Z.missM, resid.valsM), id2.miss), function (x) {
            MM <- matrix(x, ncol = ncz + M)
            z <- MM[, 1:ncz, drop = FALSE]
            res <- MM[, -(1:ncz), drop = FALSE]
            V1 <- z %*% D %*% t(z)
            diag(V1) <- diag(V1) + object$coefficients$sigma^2
            rr <- res - rowMeans(res)
            V2 <- apply(rr, 2, function (y) y %o%y)
            V2 <- if (is.matrix(V2)) rowSums(V2) / (M - 1) else sum(V2) / (M - 1)
            dim(V2) <- c(nrow(rr), nrow(rr))
            solve(chol(V1 + V2)) %*% rowMeans(res)
            }))
        }
        resid.valsM <- apply(resid.valsM, 2, function (x) {
            if (type == "stand-Subject")
                x <- x / object$coefficients$sigma
            if (type == "stand-Marginal") {
                x <- unlist(lapply(split(cbind(Z.missM, x), id2.miss), function (y) {
                    M <- matrix(y, ncol = ncz + 1)
                    z <- M[, - (ncz + 1), drop = FALSE]
                    res <- M[, ncz + 1]
                    out <- z %*% D %*% t(z)
                    diag(out) <- diag(out) + object$coefficients$sigma^2
                    solve(chol(out)) %*% res
                }))
            }
            x
        })
        names(resid.vals) <- names(fitted.vals) <- names(y)
        names(fitted.valsM) <- names(mean.resid.valsM) <- rownames(resid.valsM) <- paste("m", 1:length(fitted.valsM), sep = "")
        list("fitted.values" = fitted.vals, "residuals" = resid.vals, "fitted.valsM" = fitted.valsM, 
             "mean.resid.valsM" = mean.resid.valsM, "resid.valsM" = resid.valsM, 
             "dataM" = if (return.data) dataM else NULL)
}


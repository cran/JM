survfitJM <-
function (object, newdata, idVar = "id", survTimes = NULL, 
            last.time = NULL, M = 200, CI.levels = c(0.025, 0.975), scale = 1.6) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    if (is.null(survTimes) || !is.numeric(survTimes))
        survTimes <- seq(min(exp(object$y$logT)), max(exp(object$y$logT)), length.out = 35)
    method <- object$method
    timeVar <- object$timeVar
    id <- as.numeric(unclass(newdata[[idVar]]))
    r <- rle(id)
    id <- rep(seq_along(r$values), r$lengths)
    TermsY <- object$termsY
    mfY <- model.frame(TermsY, data = newdata)
    y <- model.response(mfY)
    X <- model.matrix(object$formYx, mfY)
    Z <- model.matrix(object$formYz, mfY)
    TermsT <- object$termsT
    data.id <- newdata[!duplicated(id), ]
    mfT <- model.frame(delete.response(TermsT), data = data.id)
    formT <- as.character(object$formT)
    formT <- as.formula(paste(formT[1], formT[3], collapse = " "))
    W <- model.matrix(formT, mfT)
    obs.times <- split(mfY[[timeVar]], id)
    last.time <- if (is.null(last.time)) {
        tapply(mfY[[timeVar]], id, tail, n = 1)
    } else if (is.character(last.time) && length(last.time) == 1) {
        tapply(newdata[[last.time]], id, tail, n = 1)
    } else if (is.numeric(last.time) && length(last.time) == nrow(data.id)) {
        last.time
    } else {
        stop("\nnot appropriate value for 'last.time' argument.")
    }
    times.to.pred <- lapply(last.time, function (t) survTimes[survTimes > t])
    n <- object$n
    n.tp <- length(last.time)
    ncx <- ncol(X)
    ncww <- ncol(W)
    betas <- object$coefficients$betas
    sigma <- object$coefficients$sigma
    D <- object$coefficients$D
    diag.D <- (ncz <- ncol(D)) == 1 & nrow(D) > 1
    D <- if (diag.D) diag(c(D)) else D
    gammas <- object$coefficients$gammas
    alpha <- object$coefficients$alpha
    if (method == "weibull-PH-GH" || method == "weibull-AFT-GH") {
        sigma.t <- object$coefficients$sigma.t
        thetas <- c(betas, log(sigma), gammas, alpha, log(sigma.t), if (diag.D) log(D) else chol.transf(D))        
    } else if (method == "piecewise-PH-GH") {
        if (ncww == 1) {
            W <- NULL
            ncww <- 0
        } else {
            W <- W[, -1, drop = FALSE]
            ncww <- ncww - 1
        }
        Q <- object$x$Q
        xi <- object$coefficients$xi
        thetas <- c(betas, log(sigma), gammas, alpha, log(xi), if (diag.D) log(D) else chol.transf(D))
    } else {
        stop("\nsurvfitJM() is not yet available for this type of joint model.")
    }
    Var.thetas <- vcov(object)
    environment(log.posterior.b) <- environment(S.b) <- environment(h.b) <- environment()
    # calculate the Empirical Bayes estimates and their (scaled) variance
    modes.b <- matrix(0, n.tp, ncz)
    Vars.b <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
        betas.new <- betas
        sigma.new <- sigma
        D.new <- D
        gammas.new <- gammas
        alpha.new <- alpha
        if (method == "weibull-PH-GH" || method == "weibull-AFT-GH") {
            sigma.t.new <- sigma.t
        } else if (method == "piecewise-PH-GH") {
            xi.new <- xi
        }
        ff <- function (b, tt, mm, i) -log.posterior.b(b, time = tt, method = mm, ii = i)
        opt <- try(optim(rep(0, ncz), ff, tt = last.time, mm = method, i = i, method = "BFGS", hessian = TRUE), TRUE)
        if (inherits(opt, "try-error")) {
            gg <- function (b, tt, mm, i) cd(b, ff, tt = tt, mm = mm, i = i)
            opt <- optim(rep(0, ncz), ff, gg, tt = last.time, mm = method, i = i, method = "BFGS", hessian = TRUE)
        } 
        modes.b[i, ] <- opt$par
        Vars.b[[i]] <- scale * solve(opt$hessian)        
    }
    out <- vector("list", M)
    success.rate <- matrix(FALSE, M, n.tp)
    b.old <- b.new <- mvrnorm(n.tp, rep(0, ncz), D)
    if (n.tp == 1)
        dim(b.old) <- dim(b.new) <- c(1, ncz)    
    for (m in 1:M) {
        # Step 1: simulate new parameter values
        thetas.new <- mvrnorm(1, thetas, Var.thetas)
        betas.new <- thetas.new[1:ncx]
        sigma.new <- exp(thetas.new[ncx + 1])
        gammas.new <- if (method == "piecewise-PH-GH" && is.null(W)) NULL else thetas.new[seq(ncx + 2, ncx + 1 + ncww)]
        alpha.new <- thetas.new[ncx + ncww + 2]
        if (method == "weibull-PH-GH" || method == "weibull-AFT-GH") {
            sigma.t.new <- exp(thetas.new[ncx + ncww + 3])
            D.new <- thetas.new[seq(ncx + ncww + 4, length(thetas.new))]
        } else if (method == "piecewise-PH-GH") {
            xi.new <- exp(thetas.new[seq(ncx + ncww + 3, ncx + ncww + 2 + Q)])
            D.new <- thetas.new[seq(ncx + ncww + Q + 3, length(thetas.new))]
        }
        D.new <- if (diag.D) exp(D.new) else chol.transf(D.new)
        SS <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            # Step 2: simulate new random effects values
            proposed.b <- rmvt(1, modes.b[i, ], Vars.b[[i]], 4)
            dmvt.old <- dmvt(b.old[i, ], modes.b[i, ], Vars.b[[i]], 4, TRUE)
            dmvt.proposed <- dmvt(proposed.b, modes.b[i, ], Vars.b[[i]], 4, TRUE)
            a <- min(exp(log.posterior.b(proposed.b, last.time, method, ii = i) + dmvt.old - 
                    log.posterior.b(b.old[i, ], last.time, method, ii = i) - dmvt.proposed), 1)
            ind <- runif(1) <= a
            success.rate[m, i] <- ind
            if (ind)
                b.new[i, ] <- proposed.b
            # Step 3: compute Pr(T > t_k | T > t_{k - 1}; theta.new, b.new)
            S.last <- S.b(last.time[i], b.new, i)
            S.pred <- sapply(times.to.pred[[i]], S.b, b = b.new, i = i)
            SS[[i]] <- S.pred / S.last
        }
        b.old <- b.new
        out[[m]] <- SS
    }
    res <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
        rr <- sapply(out, "[[", i)
        if (!is.matrix(rr))
            rr <- rbind(rr)
        res[[i]] <- cbind(
            times = times.to.pred[[i]],
            "Mean" = rowMeans(rr, na.rm = TRUE),
            "Median" = apply(rr, 1, median, na.rm = TRUE),
            "Lower" = apply(rr, 1, quantile, probs = CI.levels[1], na.rm = TRUE),
            "Upper" = apply(rr, 1, quantile, probs = CI.levels[2], na.rm = TRUE)
        )
        rownames(res[[i]]) <- as.character(seq_len(NROW(res[[i]])))
    }
    y <- split(y, id)
    names(res) <- names(y) <- names(last.time) <- names(obs.times) <- unique(unclass(newdata[[idVar]]))
    res <- list(summaries = res, survTimes = survTimes, last.time = last.time, obs.times = obs.times, y = y, ry = range(object$y$y, na.rm = TRUE),
        full.results = out, success.rate = success.rate)
    class(res) <- "survfitJM"
    res
}


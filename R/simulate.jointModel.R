simulate.jointModel <-
function (object, nsim = 1, seed = NULL, max.time = NULL, censoring = "uniform", return.ranef = FALSE, ...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    if (is.null(max.time))
        max.time <- max(exp(object$y$logT)) + 100
    # function to compute the inverse survival function
    invS <- function (t, u, i) {
        TD <- function (v) {
            # function to compute the time-dependent 
            # part for patient i at time v
            dd <- Data[rep(i, length(v)), ]
            dd[[timeVar]] <- pmax(v - lag, 0)
            if (parameterization %in% c("value", "both")) {
                XX <- model.matrix(formYx, data = dd)
                ZZ <- model.matrix(formYz, data = dd)
                Y <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), , drop = FALSE]))
                out <- alpha * Y
            }
            if (parameterization %in% c("slope", "both")) {
                XX.deriv <- model.matrix(derivForm$fixed, data = dd)
                ZZ.deriv <- model.matrix(derivForm$random, data = dd)
                Y.deriv <- as.vector(XX.deriv %*% betas[derivForm$indFixed] + 
                    rowSums(ZZ.deriv * b[rep(i, nrow(ZZ.deriv)), derivForm$indRandom, drop = FALSE]))
                out <- if (parameterization == "both") out + Dalpha * Y.deriv else Dalpha * Y.deriv
            }
            out
        }
        h <- function (s) {
            TD.i <- TD(s)
            switch(method, 
                "weibull-PH-GH" = exp(log(sigma.t) + (sigma.t - 1) * log(s) + eta.t[i] + TD.i),
                "weibull-AFT-PH" = {
                    ff <- function (v) exp(TD(v))
                    V <- exp(eta.t[i]) * integrate(ff, lower = 0, upper = s)$value
                    exp(log(sigma.t) + (sigma.t - 1) * log(V) + eta.t[i] + TD.i)
                },
                "piecewise-PH-GH" = {
                    ind <- findInterval(s, object$control$knots, rightmost.closed = TRUE)
                    xi[ind] * exp(eta.t[i] + TD.i)
                },
                "spline-PH-GH" = {
                    W2 <- splineDesign(object$control$knots, s, ord = object$control$ord, outer.ok = TRUE)
                    exp(c(W2 %*% gammas.bs) + eta.t[i] + TD.i)
                }
            )
        }
        integrate(h, lower = 0, upper = t)$value + log(u)
    }
    # extract coefficients from the fitted joint model
    betas <- object$coefficients$betas
    sigma <- object$coefficients$sigma
    D <- object$coefficients$D
    gammas <- object$coefficients$gammas
    alpha <- object$coefficients$alpha
    Dalpha <- object$coefficients$Dalpha
    sigma.t <- object$coefficients$sigma.t
    xi <- object$coefficients$xi
    gammas.bs <- object$coefficients$gammas.bs
    # extract sample characteristics from the fitted joint model
    method <- object$method
    Data <- object$data.id
    formYx <- object$formYx
    formYz <- object$formYz
    formT <- object$formT    
    TermsT <- object$termsT
    timeVar <- object$timeVar
    n <- object$n
    ncz <- NCOL(D)
    lag <- object$y$lag
    times <- object$times
    parameterization <- object$parameterization
    derivForm <- object$derivForm
    # extract design matrices from the fitted joint model
    unq.times <- sort(unique(times))
    id <- rep(1:n, each = length(unq.times))
    DD <- Data[id, ]
    DD[[timeVar]] <- rep(unq.times, n)
    X <- model.matrix(formYx, data = DD)
    Z <- model.matrix(formYz, data = DD)
    formT <- if (!is.null(kk <- attr(TermsT, "specials")$strata)) {
        strt <- eval(attr(TermsT, "variables"), data.id)[[kk]]
        tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
        reformulate(attr(tt, "term.labels"))
    } else {
        tt <- attr(delete.response(TermsT), "term.labels")
        if (length(tt)) reformulate(tt) else reformulate("1")
    }
    W <- model.matrix(formT, Data)
    if (method %in% c("spline-PH-GH", "piecewise-PH-GH"))
        W <- W[, -1, drop = FALSE]
    # simulation
    val <- vector("list", nsim)
    ranef <- vector("list", nsim)
    for (ii in seq_len(nsim)) {
        # simulate random effects
        ranef[[ii]] <- b <- mvrnorm(n, rep(0, ncz), D)
        # simulate event times
        eta.t <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, n)
        u <- runif(n)
        trueTimes <- numeric(n)    
        for (i in 1:n) {
            Root <- try(uniroot(invS, interval = c(1e-05, max.time), u = u[i], i = i)$root, TRUE)
            while(inherits(Root, "try-error")) {
                b[i, ] <- c(mvrnorm(1, rep(0, ncz), D))
                Root <- try(uniroot(invS, interval = c(1e-05, max.time), u = u[i], i = i)$root, TRUE)
            }
            trueTimes[i] <- Root
        }
        # simulate longitudinal responses
        eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ]))
        y <- rnorm(n * length(unq.times), eta.y, sigma)
        if (!is.null(censoring)) {
            if (censoring == "uniform") {
                Ctimes <- runif(n, 0, max.time)
                Time <- pmin(trueTimes, Ctimes)
                event <- as.numeric(trueTimes <= Ctimes)
            } else if (is.numeric(censoring)) {
                Ctimes <- rep(censoring, length.out = n)
                Time <- pmin(trueTimes, Ctimes)
                event <- as.numeric(trueTimes <= Ctimes)
            } else {
                warning("not appropriate value for argument 'censoring'; uniform censoring is used instead.\n")
                Ctimes <- runif(n, 0, max.time)
                Time <- pmin(trueTimes, Ctimes)
                event <- as.numeric(trueTimes <= Ctimes)
            }
        } else {
            Time <- trueTimes
            event <- rep(1, length(Time))
        }
        ind.mat <- outer(unq.times, Time, "<=")
        ind <- c(ind.mat)
        y <- y[ind]
        id <- rep(1:n, colSums(ind.mat))
        #
        DD <- Data[id, ]
        DD[[timeVar]] <- X[ind, timeVar]
        DD$y <- y
        DD$id <- id
        aa <- all.vars(TermsT)
        DD[[aa[1]]] <- Time[id]
        DD[[aa[2]]] <- event[id]
        row.names(DD) <- seq_len(nrow(DD))
        val[[ii]] <- DD
    }
    if (return.ranef)
        attr(val, "ranef") <- ranef
    attr(val, "seed") <- RNGstate
    val
}


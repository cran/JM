`jointModel` <-
function (lmeObject, survObject, timeVar, method = c("weibull-GH", "ch-GH", "ph-GH", "ch-Laplace"), 
        init = NULL, control = list()) {
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
    data.id[timeVar] <- Time
    mf <- model.frame(lmeObject$terms, data = data.id)
    Xtime <- model.matrix(formYx, mf)
    Ztime <- model.matrix(formYz, mf)
    long <- as.vector(c(Xtime %*% fixef(lmeObject)) + rowSums(Ztime * b))
    # response vectors and design matrices
    y <- list(y = y.long, logT = log(Time), d = d)
    x <- list(X = X, Xtime = Xtime, Z = Z, Ztime = Ztime, W = W)
    # extra design matrices for 'method = "ph-GH"' with event times prior to observed time for the ith subject
    if (method == "ph-GH") {
        unqT <- sort(unique(Time[d == 1]))
        times <- lapply(Time, function (t) unqT[t >= unqT])
        ind.len <- sapply(times, length)
        indT <- rep(1:nrow(data.id), ind.len)
        data.id2 <- data.id[indT, ]
        data.id2[timeVar] <- unlist(times, use.names = FALSE)
        mf <- model.frame(lmeObject$terms, data = data.id2)
        Xtime2 <- model.matrix(formYx, mf)
        Ztime2 <- model.matrix(formYz, mf)
        x <- c(x, list(Xtime2 = Xtime2, Ztime2 = Ztime2))
    }
    # control values
    con <- list(only.EM = FALSE, iter.EM = 150, iter.qN = 100, optimizer = "optim", tol1 = 1e-03, tol2 = 1e-04, 
        tol3 = sqrt(.Machine$double.eps), numeriDeriv = "fd", eps.Hes = 1e-06, parscale = NULL, step.max = 0.1, 
        backtrackSteps = 2, knots = NULL, lng.in.kn = 3, ord = 4, GHk = if (ncol(Z) < 3) 21 else 11, verbose = FALSE)
    if (method == "ph-GH") {
        con$only.EM <- TRUE
        con$iter.EM <- 200
        con$GHk <- if (ncol(Z) < 2) 21 else 11
    }
    con[(namc <- names(control))] <- control
    if (method == "ph-GH" && !con$only.EM)
        stop("with method 'ph-GH' only the EM algorithm is used.\n")
    if (method == "ph-GH" && any(!is.na(match(c("iter.qN", "optimizer"), namc))))
        warning("method 'ph-GH' uses only the EM algorithm.\n")
    # initial values
    VC <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", lmeObject$sigma^2)[[1]]
    con$inv.chol.VC <- solve(chol(solve(VC)))
    con$det.inv.chol.VC <- det(con$inv.chol.VC)
    if (all(VC[upper.tri(VC)] == 0))
        VC <- diag(VC)
    init.surv <- initial.surv(log(Time), d, if (is.null(W)) long else cbind(W, long), method, con)
    initial.values <- c(list(betas = fixef(lmeObject), sigma = lmeObject$sigma, D = VC), init.surv)
    if (!is.null(init)) {
        lngths <- lapply(initial.values[(nam.init <- names(init))], length)
        if (!is.list(init) || !isTRUE(all.equal(lngths, lapply(init, length)))) {
            warning("'init' is not a list with elements numeric vectors of appropriate length; default starting values are used instead.\n")
        } else {
            initial.values[nam.init] <- init
        }
    }
    # joint model fit
    out <- switch(method,
        "ph-GH" = phGH.fit(x, y, id, indT, ind.len, initial.values, con),
        "weibull-GH" = weibullGH.fit(x, y, id, initial.values, con),
        "ch-GH" = chGH.fit(x, y, id, initial.values, con),
        "ch-Laplace" = chLaplace.fit(x, y, id, initial.values, b, con))
    # check if any problems with the Hessian at convergence
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


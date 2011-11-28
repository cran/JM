dynC.jointModel <-
function (object, dt, data, idVar = "id", times = NULL, nt = 10, 
        simulate = FALSE, estimator = c("median", "mean"), M = 10, validate = FALSE, B = 10, verbose = FALSE, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    estimator <- match.arg(estimator)
    if (is.null(times) || !is.numeric(times)) {
        times <- sort(unique(object$times))
        if ((nn <- length(times)) > nt)
            times <- seq(times[1], times[nn], length.out = nt)
    }
    CC <- function (tt, dt, object, data, idVar) {
        data[[idVar]] <- match(data[[idVar]], unique(data[[idVar]]))
        time <- exp(object$y$logT)
        event <- object$y$d
        id <- seq_len(object$n)
        ord.ind <- order(time)
        ord.times <- time[ord.ind]
        ord.event <- event[ord.ind]
        ord.evtimes <- ord.times[ord.event == 1]
        ord.id <- id[ord.ind]
        ord.evid <- ord.id[ord.event == 1]
        i <- ord.evid[ord.evtimes > tt & ord.evtimes <= tt + dt]
        j <- ord.id[ord.times > tt + dt | (ord.times == tt + dt & ord.event == 0)]
        combs <- expand.grid(list(i = i, j = j))        
        ND <- data[data[[idVar]] %in% c(i, j), ]
        ND <- ND[ND[[object$timeVar]] <= tt, ]
        if (nrow(combs)) { 
            pi.u.t <- survfitJM(object, newdata = ND, idVar = idVar, survTimes = tt + dt, simulate = simulate, M = M)
            pi.u.t <- if (simulate) {
                sapply(pi.u.t$summaries, function (x) if (estimator == "mean") x[2] else x[3])
            } else {
                sapply(pi.u.t$summaries, "[", 2)
            }
            sf <- survfit(Surv(time, event) ~ 1)
            ind1 <- findInterval(tt, sf$time)
            ind2 <- findInterval(tt + dt, sf$time)
            S1 <- if (tt == 0) 1 else sf$surv[ind1]
            S2 <- sf$surv[ind2]
            c("AUC(t)" = mean(pi.u.t[as.character(combs$i)] < pi.u.t[as.character(combs$j)]), w = (S1 - S2) * S2, N = nrow(combs))
        } else 
            c(NA, NA, NA)
    }
    res <- t(sapply(times, CC, dt = dt, object = object, data = data, idVar = idVar))
    res[, "w"] <- res[, "w"] / sum(res[, "w"], na.rm = TRUE)
    AUCt <- cbind(t = times, "t + dt" = times + dt, res)
    out <- list(AUCt = AUCt[complete.cases(AUCt), ], dynC = c(weighted.mean(res[, 1], res[, 2], na.rm = TRUE)))
    if (validate) {
        object. <- object
        out.Boot <- out.Orig <- vector("list", B)
        for (b in 1:B) {
            if (verbose)
                cat("\nBootstrap sample:", b)
            subjects.ind <- sample(unique(data[[idVar]]), replace = TRUE)
            lis.sub <- lapply(subjects.ind, function (i) which(data[[idVar]] %in% i))
            subjects.ind <- unlist(lis.sub)
            data.new <- data[subjects.ind, ]
            data.new[["id.new"]] <- rep(seq_len(object$n), sapply(lis.sub, length))
            data.new.id <- data.new[!duplicated(data.new$id.new), ]
            #
            fixef.form <- object$formYx
            ranef.form <- paste(paste(as.character(object$formYz), collapse = ""), "|", "id.new")
            surv.form <- object$formT
            environment(fixef.form) <- environment(ranef.form) <- environment(surv.form) <- environment()
            lme.new <- lme(fixed = object$formYx, random = as.formula(ranef.form), data = data.new)
            surv.new <- coxph(object$formT, data = data.new.id, x = TRUE)
            object.new <- jointModel(lme.new, surv.new, timeVar = object$timeVar, 
                method = object$method, init = object$coefficients)
            #
            object.$coefficients <- object.new$coefficients
            res.Orig <- t(sapply(times, CC, dt = dt, object = object., data = data, idVar = idVar))
            out.Orig[[b]] <- list(dynC = cbind(times = times, res.Orig), 
                mean = weighted.mean(res.Orig[, 1], res.Orig[, 2], na.rm = TRUE))
            res.Boot <- t(sapply(times, CC, dt = dt, object = object.new, data = data.new, idVar = "id.new"))
            out.Boot[[b]] <- list(dynC = cbind(times = times, res.Boot), 
                mean = weighted.mean(res.Boot[, 1], res.Boot[, 2], na.rm = TRUE))
        }
        Optimism.dynC <- rowMeans(sapply(out.Boot, function (x) x$dynC[, 2]) - sapply(out.Orig, function (x) x$dynC[, 2]))
        Optimism.mean <- mean(sapply(out.Boot, function (x) x$mean) - sapply(out.Orig, function (x) x$mean))
        out$dynC <- cbind(out$dynC[, 1:3], optimism = Optimism.dynC, corrected = out$dynC[, 3] - Optimism.dynC, N = out$dynC[, 4])
        rownames(out$dynC) <- rep("", nrow(out$dynC))
        out$mean <- c(value = out$mean, optimism = Optimism.mean, corrected = out$mean - Optimism.mean)
    } else {
        rownames(out$AUCt) <- rep("", nrow(out$AUCt))
        out
    }
    out$validate <- validate
    out$B <- B
    out$M <- M
    out$class <- class(object)
    class(out) <- "dynC"
    out
}

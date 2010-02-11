dynC.coxph <-
function (object, dt, data, idVar = "id", timeVar = "time", times = NULL, nt = 5, validate = FALSE, ...) {
    if (!inherits(object, "coxph"))
        stop("Use only with 'coxph' objects.\n")
    if (!(type <- attr(object$y, "type")) %in% c("counting", "right"))
        stop("Use only with counting process formulation.\n")
    if (is.null(times) || !is.numeric(times)) {
        times <- sort(unique(data[[timeVar]]))
        if ((nn <- length(times)) > nt)
            times <- seq(times[1], times[nn], length.out = nt)
    }
    CC <- function (tt, dt, object, data, idVar, timeVar) {
        data[[idVar]] <- match(data[[idVar]], unique(data[[idVar]]))
        Y <- unclass(object$y)
        if (type == "counting") {
            stop.name <- colnames(Y)[2]
            time <- as.vector(sapply(split(Y[, 2], data[[idVar]]), tail, n = 1))
            event <- as.vector(sapply(split(Y[, 3], data[[idVar]]), sum))
        } else {
            stop.name <- "Time"
            time <- Y[, 1]
            event <- Y[, 2]
        }
        id <- seq_along(time)
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
        ND <- ND[ND[[timeVar]] <= tt, ]
        if (nrow(combs)) { 
            pi.u.t <- sapply(split(ND, ND[[idVar]]), function (d) {
                if (type == "counting")
                    d[nrow(d), stop.name] <- max(time)
                sfit <- survfit(object, newdata = d[1, ], individual = if (type == "counting") TRUE else FALSE)
                ind <- findInterval(tt + dt, sfit$time)
                sfit$surv[ind]
            })
            sf <- survfit(Surv(time, event) ~ 1)
            ind1 <- findInterval(tt, sf$time)
            ind2 <- findInterval(tt + dt, sf$time)
            S1 <- if (tt == 0) 1 else sf$surv[ind1]
            S2 <- sf$surv[ind2]
            c("AUC(t)" = mean(pi.u.t[as.character(combs$i)] < pi.u.t[as.character(combs$j)]), w = (S1 - S2) * S2, N = nrow(combs))
        } else 
            c(NA, NA, NA)
    }
    res <- t(sapply(times, CC, dt = dt, object = object, data = data, idVar = idVar, timeVar = timeVar))
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
            ranef.form <- paste(paste(as.character(object$formYz), collapse = ""), "|", "id.new")
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
    out$class <- class(object)
    class(out) <- "dynC"
    out
}


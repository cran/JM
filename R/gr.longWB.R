gr.longWB <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    if (parameterization %in% c("value", "both")) {
        Ys <- as.vector(Xs %*% betas) + Zsb
        eta.s <- alpha * Ys
    }
    if (parameterization %in% c("slope", "both")) {
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
        eta.s <- if (parameterization == "both") eta.s + Dalpha * Ys.deriv else Dalpha * Ys.deriv
    }
    exp.eta.tw.P <- exp(eta.tw) * P
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    Int <- wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s)
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
        ki <- exp.eta.tw.P * switch(parameterization,
            "value" = rowsum(Int * alpha * Xs[, i], id.GK, reorder = FALSE),
            "slope" = {ii <- match(i, indFixed); if (is.na(ii)) 0 else rowsum(Int * Dalpha * Xs.deriv[, ii], id.GK, reorder = FALSE)},
            "both" = {ii <- match(i, indFixed); 
                rowsum(Int * (alpha * Xs[, i] + Dalpha * if (is.na(ii)) 0 else Xs.deriv[, ii]), id.GK, reorder = FALSE)} 
        )
        kii <- c((p.byt * ki) %*% wGH)
        sc2[i] <- switch(parameterization,
            "value" = - sum(d * alpha * Xtime[, i] - kii, na.rm = TRUE),
            "slope" = {ii <- match(i, indFixed); 
                if (is.na(ii)) 0 else - sum(d * Dalpha * Xtime.deriv[, ii] - kii, na.rm = TRUE)},
            "both" = {ii <- match(i, indFixed); 
                - sum(d * (alpha * Xtime[, i] + Dalpha * if (is.na(ii)) 0 else Xtime.deriv[, ii]) - kii, na.rm = TRUE)}
        )
    }    
    c(sc1 + sc2)
}


logsurvCH <-
function (alpha, bb) {
    eta.ti <- eta.twi + alpha * (eta.yxT[i] + sum(Ztime.i * bb))
    if (d[i]) log(sc[i]) + eta.ti - exp(eta.ti) - logT[i] else - exp(eta.ti)
}


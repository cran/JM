`gr.b` <-
function (bb) {
    gr.ybi <- - crossprod(Z.ind.i, Z.ind.i %*% bb - yi.eta.yxi) / sigma^2
    eta.ti <- eta.twi + alpha * (eta.yxT[i] + sum(Ztime.i * bb))
    gr.tbi <- alpha * (d[i] - exp(eta.ti) ) * cbind(Ztime.i)
    gr.bi <- if (diag.D) - bb / D else - solve(D, bb)
    - as.vector(gr.ybi + gr.tbi + gr.bi)
}


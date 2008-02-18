`residuals.jointModel` <-
function (object, process = c("Longitudinal", "Event"), 
    type = c("Marginal", "Subject", "stand-Marginal", "stand-Subject"), ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    process <- match.arg(process)
    type <- match.arg(type)
    if (process == "Longitudinal") {
        fits <- if (type == "Marginal" || type == "stand-Marginal") {
            fitted(object, process = "Longitudinal", type = "Marginal")
        } else {
            fitted(object, process = "Longitudinal", type = "Subject")
        }
        if (type == "Marginal" || type == "Subject") {
            object$y$y - fits
        } else {
            (object$y$y - fits) / object$coefficients$sigma
        }
    } else {
        cat("the residuals() method is not currently implemented for the Event process.\n")
    }
}


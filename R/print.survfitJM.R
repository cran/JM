print.survfitJM <-
function (x, ...) {
    if (!is.null(x$success.rate)) {
        cat("\nPrediction of Conditional Probabilities for Event based on", nrow(x$success.rate), "replications\n\n")
    } else {
        cat("\nPrediction of Conditional Probabilities for Events\n\n")
    }
    print(lapply(x$summaries, round, 4))
    invisible(x)
}


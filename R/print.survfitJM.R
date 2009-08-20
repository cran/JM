print.survfitJM <-
function (x, ...) {
        cat("\nPrediction of Conditional Probabilities for Event based on", nrow(x$success.rate), "replications\n\n")
    print(lapply(x$summaries, round, 4))
    invisible(x)
}


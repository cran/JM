print.dynC <-
function (x, ...) {
    if (!inherits(x, "dynC"))
        stop("Use only with 'dynC' objects.\n")
    if (x$class == "jointModel") 
        cat("\n\tDynamic Prediction for Joint Models\n")
    else
        cat("\n\tDynamic Prediction for Cox Models\n")
    cat("\nDynamic AUC:\n")
    print(round(x$AUCt, 4))
    cat("\nDynamic C index:", round(x$dynC, 4), "\n")
    cat("\n\n")
    invisible(x)
}

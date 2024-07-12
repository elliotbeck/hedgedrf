#' hedgedrf prediction
#'
#' @param object hedgedrf \code{hedgedrf} object.
#' @param data data New test data of class \code{data.frame} or
#' \code{gwaa.data} (GenABEL).
#' @param ... Additional arguments to pass to the \code{predict.ranger}
#' function.
#' @return The hedged random forest predictions. An object of class \code{matrix}.
#' @export
predict.hedgedrf <- function(object, data, ...) {
    # Random forest predictions on test data
    rf_predictions <- predict(
        object$rf_fit,
        data,
        predict.all = TRUE,
        ...
    )$predictions
    hedgedrf_predictions <- as.matrix(rf_predictions) %*% object$tree_weights
    return(hedgedrf_predictions)
}

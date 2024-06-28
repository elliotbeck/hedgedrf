#' hedgedrf
#'
#' @param  formula Object of class \code{formula} or \code{character}
#' describing the model to fit. Interaction terms supported only for
#' numerical variables.
#' @param data Training data of class \code{data.frame}, \code{matrix},
#' \code{dgCMatrix} (Matrix) or \code{gwaa.data} (GenABEL).
#' @param x Predictor data (independent variables), alternative interface to data with formula or
#' dependent.variable.name.
#' @param y Response vector (dependent variable), alternative interface to data with formula or
#' dependent.variable.name. For survival use a Surv() object or a matrix with time and status.
#' @param num_iter Number of iterations for the optimization algorithm.
#' @param kappa Amount of regularization to apply to the tree weights. 1 implies
#' no shorting, 2 implies no more than 50% shorting, etc.
#' @param ... Additional arguments to pass to the \code{ranger} function.
#'
#' @return An object of class \code{hedgedrf} containing the tree weights and
#' a ranger object. The tree weights can be used to construct a hedged random
#' forest with the \code{predict.hedgedrf} function. For more details about the
#' ranger object, see the ranger documentation.
#' @importFrom stats cov predict

#' @examples
#' rf <- hedgedrf(mpg ~ ., mtcars[1:26, ])
#' pred <- predict(rf, mtcars[27:32, ])
#' pred
#'
#' @export
#'
hedgedrf <- function(
    formula = NULL, data = NULL, x = NULL,
    y = NULL, num_iter = NULL, kappa = 2, ...) {
    if (kappa <= 0 || !is.numeric(kappa)) {
        stop("Kappa must be a positive real number")
    }
    if (is.factor(data[, all.vars(formula)[1]]) ||
        is.logical(data[, all.vars(formula)[1]])) {
        stop("Response variable must be numeric")
    }

    # Fit random forest
    rf_fit <- ranger::ranger(formula = formula, data = data, x = x, y = y, ...)

    # Random forest predictions on all trees
    rf_predictions_all <- predict(
        rf_fit,
        if (is.null(data)) x else data,
        predict.all = TRUE
    )$predictions

    # Calculate residuals
    if (is.null(data)) {
        rf_residuals <- y - rf_predictions_all
    } else {
        rf_residuals <- data[, all.vars(formula)[1]] -
            rf_predictions_all
    }

    # Calculate covariance matrix and mean vector of residuals
    mean_vector <- colMeans(rf_residuals)
    cov_matrix <- cov(rf_residuals)
    cov_matrix <- get_cov_qis(rf_residuals)

    # Get optimized tree weights
    w <- CVXR::Variable(ncol(cov_matrix))
    objective <- CVXR::square((t(w) %*% mean_vector)) +
        CVXR::quad_form(w, cov_matrix)
    constraints <- list(
        sum(w) == 1,
        sum(abs(w)) <= kappa
    )
    prob <- CVXR::Problem(CVXR::Minimize(objective), constraints)
    solution <- CVXR::solve(prob, num_iter = num_iter)

    hedgedrf <- structure(
        list(
            tree_weights = solution$getValue(w),
            rf_fit = rf_fit
        ),
        class = "hedgedrf"
    )
    return(hedgedrf)
}

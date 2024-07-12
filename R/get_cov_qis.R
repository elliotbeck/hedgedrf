#' Quadratic-inverse shrinkage
#'
#' Nonlinear shrinkage derived under Frobenius loss and its two cousins,
#' Inverse Stein’s loss and Minimum Variance loss, called quadratic-inverse
#' shrinkage (QIS). See Ledoit and Wolf (2022, Section 4.5).
#'
#' @param  data (n*p): raw data matrix of n iid observations on p random
#' variables
#' @param k If k < 0, then the algorithm demeans the data by default, and
#' adjusts the effective sample size accordingly. If the user inputs k = 0,
#' then no demeaning takes place; if user inputs k = 1, then it signifies that
#' the data data have already been demeaned.
#'
#' @return sigmahat (p*p): the QIS covariance matrix estimate. An object of
#' class \code{matrix}.
#'
get_cov_qis <- function(data, k = -1) {
    dim.data <- dim(data)
    n <- dim.data[1]
    p <- dim.data[2]
    if (k < 0) {
        # demean the data and set k = 1
        data <- scale(data, scale = FALSE)
        k <- 1
    }
    n_free <- n - k
    c <- p / n_free
    sample <- (t(data) %*% data) / n_free
    sample <- (t(sample) + sample) / 2
    spectral <- eigen(sample, symmetric = TRUE)
    lambda <- spectral$values[p:1]
    u <- spectral$vectors[, p:1]
    h <- min(c^2, 1 / c^2)^0.35 / p^0.35
    inv_lambda <- 1 / lambda[max(1, p - n_free + 1):p]
    l_j <- matrix(rep(inv_lambda, each = min(p, n_free)), nrow = min(p, n_free))
    l_j_i <- l_j - t(l_j)
    theta <- rowMeans(l_j * l_j_i / (l_j_i^2 + h^2 * l_j^2))
    h_theta <- rowMeans(l_j * (h * l_j) / (l_j_i^2 + h^2 * l_j^2))
    a_theta2 <- theta^2 + h_theta^2
    if (p <= n_free) {
        delta <- 1 / ((1 - c)^2 * inv_lambda + 2 * c * (1 - c) * inv_lambda *
            theta + c^2 * inv_lambda * a_theta2)
    } else {
        delta0 <- 1 / ((c - 1) * mean(inv_lambda))
        delta <- c(rep(delta0, p - n_free), 1 / (inv_lambda * a_theta2))
    }
    delta_qis <- delta * (sum(lambda) / sum(delta))
    sigma_hat <- u %*% diag(delta_qis) %*% t(u)
    return(sigma_hat)
}

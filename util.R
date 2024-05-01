library(mvtnorm)

# Function to check if a symmetric matrix is positive definite
## Input
### matrix: a symmetric matrix
### tol: tolerance
## Output
### A true or false value
check_positive_definite <- function(matrix, tol = 1e-8) {
  eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvalues <= tol)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# Function to generate the correlation matrix
## Input
### t: information time
## Output
### Correlation matrix
cr_function <- function(t){
  K <- length(t)
  cr <- diag(K)
  for (i in 1:K){
    for (j in i:K){
      cr[i, j] <- sqrt(t[i] / t[j])
    }
  }
  cr <- cr + t(cr) - diag(K)
  return(cr)
}

# Error spending function for the OBrien-Fleming Lan-DeMets design
## Input
### alpha: one-sided significance level
### t: information time
## Output
### Vector of cumulative error spent
esf_OBF_function <- function(alpha, t) {
  y <- 2 - 2 * pnorm(qnorm(1 - alpha / 2) / sqrt(t))
  return(y)
}

# Error spending function for the Pocock Lan-DeMets design
## Input
### alpha: one-sided significance level
### t: information time
## Output
### Vector of cumulative error spent
esf_POC_function <- function(alpha, t) {
  y <- alpha * log(1 + (exp(1) - 1) * t)
  return(y)
}

# Function to derive the stopping boundary using the cumulative error spent
## Input
### t: information time
### cumulative: cumulative error spent with the same length as t
## Output
### Vector of boundary values
solver_boundary_esf_function <- function(t, cumulative, steps = 1024){
  K <- length(t)
  cr <- cr_function(t)
  solver <- function(x, cumu, cr = cr, c_past){
    z <- c(c_past, x)
    return(1 - cumu - pmvnorm(upper = z, corr = cr[1:length(z), 1:length(z)],
                              algorithm = Miwa(steps = steps, checkCorr = FALSE))
    )
  }
  c_boundary <- rep(0,K)
  c_boundary[1] <- min(qnorm(1 - cumulative[1]), 10)
  if (K > 1){
    for (k in 2:K){
      a <- uniroot(solver, interval=c(0.001, 10), cumu = cumulative[k],
                   cr = cr, c_past = c_boundary[1:(k-1)])$root
      c_boundary[k] <- min(a, 10)
    }
  }
  return(c_boundary)
}

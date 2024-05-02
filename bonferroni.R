
alpha_BAG_delta <- function(R12, R135, R1246, delta, info, c_bound_1, c_bound_2, d_bound, steps = 512) {
  mean <- delta * c(sqrt(info), 0, 0, 0, 0)
  part1 <- pnorm(c_bound_1, mean = mean[1], sd = 1, lower.tail = F) -
    pmvnorm(lower = c(c_bound_1, -Inf, -Inf),
            upper = c(Inf, d_bound[1], d_bound[1]),
            mean = mean[c(1, 3, 5)], corr = R135,
            algorithm = Miwa(steps = steps, checkCorr = F), keepAttr = FALSE)
  part2 <- pmvnorm(lower = c(-Inf, c_bound_2),
                   upper = c(c_bound_1, Inf),
                   mean = mean[c(1, 2)], corr = R12,
                   algorithm = Miwa(steps = steps, checkCorr = F), keepAttr = FALSE) -
    pmvnorm(lower = c(-Inf, c_bound_2, -Inf, -Inf),
            upper = c(c_bound_1, Inf, d_bound[2], d_bound[2]),
            mean = mean[c(1, 2, 4, 6)], corr = R1246,
            algorithm = Miwa(steps = steps, checkCorr = F), keepAttr = FALSE)
  return(part1 + part2)
}

alpha_BAG <- function(R12, R135, R1246, info, c_bound_1, c_bound_2, d_bound_type, x) {
  # Secondary endpoint boundary
  if (d_bound_type == "OBF") {
    cumu_OBF <- esf_OBF_function(x / 2, info)
    d_bound <- solver_boundary_esf_function(info, cumu_OBF),
  } else if (d_bound_type == "POC") {
    cumu_POC <- esf_POC_function(x / 2, info)
    d_bound <- solver_boundary_esf_function(info, cumu_POC)
  } else {
    stop("Only 'OBF' and 'POC' are supported")
  }
  out <- optimize(alpha_BAG_delta, c(0, 10), maximum = T,
                  R12 = R12, R135 = R135, R1246 = R1246, info = info,
                  c_bound_1 = c_bound_1, c_bound_2 = c_bound_2, d_bound = d_bound)
  return(out)
}

search_alpha_BAG <- function(i, scen) {
  rho_xy <- scen$rho_xy[i]
  rho_xz <- scen$rho_xz[i]
  rho_yz <- scen$rho_yz[i]
  t <- scen$t[i]
  c_bound_1 <- scen$c_bound_1[i]
  c_bound_2 <- scen$c_bound_2[i]
  d_bound_type <- scen$d_bound_type[i]
  alpha <- scen$alpha[i]
  increment <- scen$increment[i]
  info <- c(t, 1)
  # Correlation structure
  R <- cbind(rbind(cr_function(info), rho_xy * cr_function(info), rho_xz * cr_function(info)),
             rbind(rho_xy * cr_function(info), cr_function(info), rho_yz * cr_function(info)),
             rbind(rho_xz * cr_function(info), rho_yz * cr_function(info), cr_function(info))
  )
  R12 <- R[c(1, 2), c(1, 2)]
  R135 <- R[c(1, 3, 5), c(1, 3, 5)]
  R1246 <- R[c(1, 2, 4, 6), c(1, 2, 4, 6)]
  if (check_positive_definite(R12) & check_positive_definite(R135) &
      check_positive_definite(R1246)) {
    temp <- 0
    x <- alpha
    while (temp < alpha + .Machine$double.eps) {
      x <- x + increment
      temp <- alpha_BAG(R12, R135, R1246, info, c_bound_1, c_bound_2, d_bound_type, alpha, x)$objective
      if (temp > alpha + .Machine$double.eps | x > 2 * alpha) {
        break
      }
    }
    x <- max(alpha, x - increment)
    alpha_BAG_opt <- alpha_BAG(R12, R135, R1246, info, c_bound_1, c_bound_2, d_bound_type, alpha, x)
    out <- c(rho_xy, rho_xz, rho_yz, t, c_bound_1, c_bound_2, alpha,
             x, unlist(alpha_BAG_opt))
  } else {
    out <- c(rho_xy, rho_xz, rho_yz, t, c_bound_1, c_bound_2, alpha,
             NA, NA, NA)
  }
  names(out) <- c("rho_xy", "rho_xz", "rho_yz", "t", "c_bound_1", "c_bound_2",
                  "alpha", "alpha_2", "delta", "alpha_BAG")
  return(out)
}


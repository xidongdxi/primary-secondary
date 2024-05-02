
source("/primary-secondary/util.R")
source("/primary-secondary/bonferroni.R")

## Simulation scenarios
rho_xy <- seq(0, 1, 0.01)
rho_xz <- seq(0, 1, 0.01)
rho_yz <- seq(0, 1, 0.01)
alpha <- 0.025
t <- 0.5
info <- c(t, 1)
solver_boundary_esf_function(info, esf_OBF_function(alpha, info))
c_bound_1 <- solver_boundary_esf_function(info, esf_OBF_function(alpha, info))[1]
c_bound_2 <- solver_boundary_esf_function(info, esf_OBF_function(alpha, info))[2]
d_bound_type <- "POC"
increment <- 1e-4
scen <- expand.grid(rho_xy, rho_xz, rho_yz, t, c_bound_1, c_bound_2, d_bound_type, alpha, increment)
colnames(scen) <- c("rho_xy", "rho_xz", "rho_yz", "t", "c_bound_1", "c_bound_2",
                    "d_bound_type", "alpha", "increment")
n_scen <- nrow(scen)

# Parallel search
library(future.apply)
# plan(multisession)
plan(cluster, workers = 124)
seed <- 10000
start_time <- Sys.time()
result <- future_lapply(1:n_scen, FUN = search_alpha_BAG, future.seed = seed,
                        future.packages = c("mvtnorm"), scen = scen)
end_time <- Sys.time()
end_time - start_time

results <- as.data.frame(do.call(rbind, result))
results <- results[order(results$alpha_2), ]
head(results)
source("/primary-secondary/util.R")
source("/primary-secondary/simes.R")

## Simulation scenarios
rho_xy <- seq(0, 0.9, 0.01)
rho_xz <- seq(0, 0.9, 0.01)
rho_yz <- seq(0, 0.9, 0.01)
t <- 0.5
c_bound_1 <- 2.537988
c_bound_2 <- 1.662092
d_bound_type <- "POC"
alpha <- 0.025
increment <- 1e-4
scen <- expand.grid(rho_xy, rho_xz, rho_yz, t, c_bound_1, c_bound_2, d_bound_type, alpha, increment)
colnames(scen) <- c("rho_xy", "rho_xz", "rho_yz", "t", "c_bound_1", "c_bound_2",
                    "d_bound_type", "alpha", "increment")
n_scen <- nrow(scen)

# Parallel
library(future.apply)
# plan(multisession)
plan(cluster, workers = 124)
seed <- 10000
start_time <- Sys.time()
result <- future_lapply(1:n_scen, FUN = search_alpha_SAG, future.seed = seed,
                        future.packages = c("mvtnorm"), scen = scen)
end_time <- Sys.time()
end_time - start_time

results <- as.data.frame(do.call(rbind, result))
write.csv(results, file = "results_simes_1e4_0.01_1e-5_0.025.csv")

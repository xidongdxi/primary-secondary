R code for the paper "Testing One Primary and Two Secondary Endpoints in a Two-Stage Group Sequential Trial with Extensions"

# Files included
**util.R**: R code containing functions for group sequential designs. Requires the `mvtnorm` package.
**bonferroni.R**: R code containing functions for Bonferroni-based normal theory boundaries.
**bonferroni_grid_search.R**: R code containing functions for the grid search of optimal boundaries for Bonferroni-based normal theory boundaries. The current code is configured for the computing instances on AWS.
**simes.R**: R code containing functions for Simes-based normal theory boundaries.
**simes_grid_search.R**: R code containing functions for the grid search of optimal boundaries for Simes-based normal theory boundaries. The current code is configured for the computing instances on AWS.

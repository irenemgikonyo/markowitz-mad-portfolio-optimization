# PACKAGES
required_packages = c("quantmod", "xts", "quadprog", "lpSolve",
                      "ggplot2", "dplyr", "tidyr", "scales")

installed = rownames(installed.packages())
for (p in required_packages) {
  if (!(p %in% installed)) install.packages(p)
}

library(quantmod)
library(xts)
library(quadprog)
library(lpSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

set.seed(2026)

# SETTINGS
tol = 1e-7
maximum_weight = 0.30

assets = c("AAPL", "MSFT", "AMZN", "GOOGL", "META", "NVDA",
           "JPM", "BAC", "JNJ", "PFE", "KO", "PEP", "WMT",
           "XOM", "CVX", "CAT", "DIS", "NFLX", "V", "MA")

start_date = "2021-01-01"
training_end = "2024-12-31"
testing_start = "2025-01-01"
testing_end = "2025-12-31"

# DATA DOWNLOAD
prices_list = list()

for (ticker in assets) {
  x = getSymbols(ticker, src = "yahoo", from = start_date,
                 to = "2026-01-01", auto.assign = FALSE)
  prices_list[[ticker]] = Ad(x)
}

prices_daily = do.call(merge, prices_list)
colnames(prices_daily) = assets
prices_daily = na.omit(prices_daily)

# WEEKLY PRICES AND RETURNS
prices_weekly = na.omit(apply.weekly(prices_daily, last))
returns_weekly = na.omit(prices_weekly / lag(prices_weekly) - 1)

R = coredata(returns_weekly)
colnames(R) = assets

# TRAINING / TESTING SPLIT
training_dates =
  index(returns_weekly) <= as.Date(training_end)

testing_dates =
  index(returns_weekly) >= as.Date(testing_start) &
  index(returns_weekly) <= as.Date(testing_end)

R_train = coredata(returns_weekly[training_dates])
R_test = coredata(returns_weekly[testing_dates])

colnames(R_train) = assets
colnames(R_test) = assets

n_assets = ncol(R_train)
T_train = nrow(R_train)
T_test = nrow(R_test)

if (n_assets != length(assets)) stop("Unexpected number of assets.")
if (T_train < 2) stop("Too few training observations.")
if (T_test < 2) stop("Too few testing observations.")

cat("\n===== DATA SPLIT =====\n")
cat("Training:", start_date, "to", training_end, "\n")
cat("Testing:", testing_start, "to", testing_end, "\n")
cat("Assets:", n_assets, "\n")
cat("Training observations:", T_train, "\n")
cat("Testing observations:", T_test, "\n")



# TRAINING-SET ASSET STATISTICS
mu = colMeans(R_train)
Sigma = cov(R_train)
names(mu) = assets

weekly_sd = apply(R_train, 2, sd)
annualised_return = (1 + mu)^52 - 1
annualised_volatility = weekly_sd * sqrt(52)

asset_statistics = data.frame(
  Asset = assets,
  Weekly_Mean = mu,
  Weekly_SD = weekly_sd,
  Annualised_Return = annualised_return,
  Annualised_Volatility = annualised_volatility
) %>% arrange(desc(Weekly_Mean))

if (anyNA(R_train) || anyNA(mu) || anyNA(Sigma))
  stop("Missing observations remain in training data.")

stopifnot(
  ncol(R_train) == length(assets),
  nrow(Sigma) == length(assets),
  ncol(Sigma) == length(assets),
  length(mu) == length(assets)
)

# CORRELATION ANALYSIS

cor_matrix = cor(R_train)

cor_long = as.data.frame(as.table(cor_matrix))
colnames(cor_long) = c("Asset1", "Asset2", "Correlation")

cor_high = cor_long %>%
  mutate(Asset1 = as.character(Asset1),
         Asset2 = as.character(Asset2)) %>%
  filter(Asset1 < Asset2) %>%
  arrange(desc(Correlation))

# INDIVIDUAL ASSET RETURN AND VOLATILITY

p1 = ggplot(asset_statistics,
            aes(x = Weekly_SD, y = Weekly_Mean, label = Asset)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  labs(title = "Weekly Return and Volatility of Individual Assets",
       x = "Weekly Standard Deviation",
       y = "Mean Weekly Return") +
  theme_minimal()

print(p1)

# COVARIANCE MATRIX VALIDATION

validate_covariance = function(Sigma) {
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma))
    stop("Sigma must be a square matrix.")
  
  if (!isTRUE(all.equal(Sigma, t(Sigma), tolerance = 1e-10)))
    stop("Covariance matrix is not symmetric.")
  
  eigenvalues = eigen(Sigma, symmetric = TRUE,
                      only.values = TRUE)$values
  
  if (min(eigenvalues) <= 0)
    stop("Covariance matrix is not positive definite.")
  
  invisible(TRUE)
}

validate_covariance(Sigma)



# MARKOWITZ OPTIMISATION


solve_markowitz = function(target_return, mu, Sigma,
                           max_weight = 0.30, tol = 1e-7) {
  
  n = length(mu)
  
  if (!is.matrix(Sigma) || nrow(Sigma) != n || ncol(Sigma) != n)
    stop(paste0("mu and Sigma dimensions do not agree. length(mu) = ",
                n, ", Sigma = ", nrow(Sigma), " x ", ncol(Sigma)))
  
  if (max_weight <= 0 || max_weight > 1)
    stop("max_weight must be greater than 0 and no greater than 1.")
  
  validate_covariance(Sigma)
  
  Dmat = 2 * Sigma
  dvec = rep(0, n)
  
  Amat = cbind(rep(1, n), mu, diag(n), -diag(n))
  bvec = c(1, target_return, rep(0, n), rep(-max_weight, n))
  
  result = tryCatch(
    solve.QP(Dmat = Dmat, dvec = dvec,
             Amat = Amat, bvec = bvec, meq = 2),
    error = function(e) {
      warning(paste("quadprog failed:", e$message))
      NULL
    }
  )
  
  if (is.null(result)) return(NULL)
  
  weights = as.numeric(result$solution)
  names(weights) = names(mu)
  weights[abs(weights) < tol] = 0
  
  achieved_return = sum(weights * mu)
  
  feasible = all(
    abs(sum(weights) - 1) <= 1e-6,
    abs(achieved_return - target_return) <= 1e-6,
    min(weights) >= -tol,
    max(weights) <= max_weight + tol
  )
  
  if (!feasible) return(NULL)
  
  weights
}



# MAD OPTIMISATION


solve_mad = function(target_return, R, mu, max_weight = 0.30) {
  
  T_obs = nrow(R)
  n = ncol(R)
  
  if (length(mu) != n) stop("R and mu dimensions do not agree.")
  if (T_obs < 2) stop("At least two observations are required.")
  if (max_weight <= 0 || max_weight > 1)
    stop("max_weight must be in (0, 1].")
  
  total_variables = n + T_obs
  
  objective = c(rep(0, n), rep(1 / T_obs, T_obs))
  
  constraint_list = list()
  direction_list = character()
  rhs_list = numeric()
  
  # Budget constraint
  row_sum = rep(0, total_variables)
  row_sum[1:n] = 1
  constraint_list[[length(constraint_list) + 1]] = row_sum
  direction_list = c(direction_list, "=")
  rhs_list = c(rhs_list, 1)
  
  # Target return constraint
  row_return = rep(0, total_variables)
  row_return[1:n] = mu
  constraint_list[[length(constraint_list) + 1]] = row_return
  direction_list = c(direction_list, "=")
  rhs_list = c(rhs_list, target_return)
  
  # Absolute deviation constraints
  for (t in seq_len(T_obs)) {
    
    deviation = R[t, ] - mu
    
    row_positive = rep(0, total_variables)
    row_positive[1:n] = deviation
    row_positive[n + t] = -1
    
    constraint_list[[length(constraint_list) + 1]] = row_positive
    direction_list = c(direction_list, "<=")
    rhs_list = c(rhs_list, 0)
    
    row_negative = rep(0, total_variables)
    row_negative[1:n] = -deviation
    row_negative[n + t] = -1
    
    constraint_list[[length(constraint_list) + 1]] = row_negative
    direction_list = c(direction_list, "<=")
    rhs_list = c(rhs_list, 0)
  }
  
  # Maximum weight constraints
  for (i in seq_len(n)) {
    row_upper = rep(0, total_variables)
    row_upper[i] = 1
    
    constraint_list[[length(constraint_list) + 1]] = row_upper
    direction_list = c(direction_list, "<=")
    rhs_list = c(rhs_list, max_weight)
  }
  
  # Non-negativity of deviation variables
  for (t in seq_len(T_obs)) {
    row_nonnegative = rep(0, total_variables)
    row_nonnegative[n + t] = 1
    
    constraint_list[[length(constraint_list) + 1]] = row_nonnegative
    direction_list = c(direction_list, ">=")
    rhs_list = c(rhs_list, 0)
  }
  
  const_matrix = do.call(rbind, constraint_list)
  
  result = tryCatch(
    lp(direction = "min",
       objective.in = objective,
       const.mat = const_matrix,
       const.dir = direction_list,
       const.rhs = rhs_list,
       all.int = FALSE,
       all.bin = FALSE,
       compute.sens = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(result) || result$status != 0) return(NULL)
  
  weights = as.numeric(result$solution[1:n])
  names(weights) = colnames(R)
  weights[abs(weights) < tol] = 0
  
  achieved_return = sum(weights * mu)
  
  feasible = all(
    abs(sum(weights) - 1) <= 1e-6,
    abs(achieved_return - target_return) <= 1e-6,
    min(weights) >= -tol,
    max(weights) <= max_weight + tol
  )
  
  if (!feasible) return(NULL)
  
  weights
}



# FEASIBLE EXPECTED-RETURN RANGE


calculate_extreme_return = function(sorted_mu, max_weight) {
  
  n = length(sorted_mu)
  
  if (max_weight * n < 1 - tol)
    stop("Maximum-weight constraint makes budget infeasible.")
  
  remaining = 1
  portfolio_return = 0
  
  for (x in sorted_mu) {
    allocation = min(max_weight, remaining)
    portfolio_return = portfolio_return + allocation * x
    remaining = remaining - allocation
    
    if (remaining <= tol) break
  }
  
  if (remaining > tol)
    stop("Could not allocate the full portfolio.")
  
  portfolio_return
}

minimum_feasible_return =
  calculate_extreme_return(sort(mu), maximum_weight)

maximum_feasible_return =
  calculate_extreme_return(sort(mu, decreasing = TRUE),
                           maximum_weight)

if (minimum_feasible_return >= maximum_feasible_return)
  stop("Feasible return range is empty.")



# GLOBAL MINIMUM-VARIANCE PORTFOLIO


solve_gmv = function(mu, Sigma, max_weight = 0.30, tol = 1e-7) {
  
  n = length(mu)
  
  Dmat = 2 * Sigma
  dvec = rep(0, n)
  
  Amat = cbind(rep(1, n), diag(n), -diag(n))
  bvec = c(1, rep(0, n), rep(-max_weight, n))
  
  result = tryCatch(
    solve.QP(Dmat = Dmat, dvec = dvec,
             Amat = Amat, bvec = bvec, meq = 1),
    error = function(e) NULL
  )
  
  if (is.null(result)) return(NULL)
  
  weights = as.numeric(result$solution)
  names(weights) = names(mu)
  weights[abs(weights) < tol] = 0
  
  weights
}

gmv_weights = solve_gmv(mu, Sigma, maximum_weight, tol)

if (is.null(gmv_weights))
  stop("GMV optimisation failed.")

gmv_return = sum(gmv_weights * mu)
gmv_sd = sqrt(as.numeric(t(gmv_weights) %*% Sigma %*% gmv_weights))



# TARGET RETURN GRID


target_margin = 0.000001

frontier_min_return = gmv_return + target_margin
frontier_max_return = maximum_feasible_return - target_margin

if (frontier_min_return >= frontier_max_return)
  stop("Efficient frontier range is empty.")

target_grid = seq(frontier_min_return,
                  frontier_max_return,
                  length.out = 50)



# SOLVER VALIDATION


markowitz_test = solve_markowitz(target_grid[1], mu, Sigma,
                                 maximum_weight)

mad_test = solve_mad(target_grid[1], R_train, mu,
                     maximum_weight)

if (is.null(markowitz_test) || is.null(mad_test))
  stop("One or more optimisation solvers failed validation.")



# EFFICIENT FRONTIERS


frontier_list = list()
markowitz_weights_list = list()
mad_weights_list = list()
successful_targets = numeric()

for (target in target_grid) {
  
  mw = solve_markowitz(target, mu, Sigma, maximum_weight)
  md = solve_mad(target, R_train, mu, maximum_weight)
  
  if (!is.null(mw) && !is.null(md)) {
    
    mw_returns = as.numeric(R_train %*% mw)
    md_returns = as.numeric(R_train %*% md)
    
    mw_return = sum(mw * mu)
    md_return = sum(md * mu)
    
    mw_sd = sd(mw_returns)
    md_sd = sd(md_returns)
    
    frontier_list[[length(frontier_list) + 1]] =
      data.frame(Target = target,
                 Markowitz_Return = mw_return,
                 Markowitz_SD = mw_sd,
                 MAD_Return = md_return,
                 MAD_SD = md_sd)
    
    markowitz_weights_list[[length(markowitz_weights_list) + 1]] = mw
    mad_weights_list[[length(mad_weights_list) + 1]] = md
    successful_targets = c(successful_targets, target)
  }
}

frontier_results = bind_rows(frontier_list)

if (nrow(frontier_results) < 2)
  stop("Too few feasible frontier points.")

target_grid = successful_targets



# RISK COMPARISON


frontier_results = frontier_results %>%
  mutate(SD_Difference = MAD_SD - Markowitz_SD,
         SD_Ratio = MAD_SD / Markowitz_SD)

mean_sd_difference = mean(frontier_results$SD_Difference)
mean_sd_ratio = mean(frontier_results$SD_Ratio)

mean_markowitz_frontier_sd = mean(frontier_results$Markowitz_SD)
mean_mad_frontier_sd = mean(frontier_results$MAD_SD)

percentage_sd_difference =
  ((mean_mad_frontier_sd - mean_markowitz_frontier_sd) /
     mean_markowitz_frontier_sd) * 100



# FRONTIER COMPARISON


frontier_comparison = data.frame(
  Mean_SD_Difference = mean(frontier_results$SD_Difference),
  Median_SD_Difference = median(frontier_results$SD_Difference),
  Mean_SD_Ratio = mean(frontier_results$SD_Ratio),
  Mean_Markowitz_SD = mean_markowitz_frontier_sd,
  Mean_MAD_SD = mean_mad_frontier_sd,
  Percentage_SD_Difference = percentage_sd_difference
)



# OBJECTIVE VERIFICATION


mad_objective_check = data.frame(
  Target = target_grid,
  Portfolio_MAD = sapply(mad_weights_list, function(w) {
    portfolio_returns = as.numeric(R_train %*% w)
    portfolio_mean = sum(w * mu)
    mean(abs(portfolio_returns - portfolio_mean))
  })
)

markowitz_objective_check = data.frame(
  Target = target_grid,
  Portfolio_Variance = sapply(markowitz_weights_list, function(w) {
    as.numeric(t(w) %*% Sigma %*% w)
  })
)

print(mad_objective_check)
print(markowitz_objective_check)



# RISK-RETURN FRONTIERS


frontier_long = bind_rows(
  frontier_results %>%
    transmute(Target, Model = "Markowitz",
              Return = Markowitz_Return, SD = Markowitz_SD),
  frontier_results %>%
    transmute(Target, Model = "MAD",
              Return = MAD_Return, SD = MAD_SD)
) %>% arrange(Model, Target)

p2 = ggplot(frontier_long, aes(x = SD, y = Return, colour = Model)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.8) +
  labs(title = "Markowitz and MAD Risk-Return Frontiers",
       x = "Weekly Portfolio Risk: Standard Deviation",
       y = "Weekly Expected Return",
       colour = "Model") +
  theme_minimal()

print(p2)



# COMMON TARGET PORTFOLIO


target_index = ceiling(length(target_grid) / 2)
common_target = target_grid[target_index]

m_weights = markowitz_weights_list[[target_index]]
mad_weights = mad_weights_list[[target_index]]



# TRAINING-PERIOD PERFORMANCE


markowitz_expected_return = sum(m_weights * mu)
mad_expected_return = sum(mad_weights * mu)

markowitz_portfolio_returns = as.numeric(R_train %*% m_weights)
mad_portfolio_returns = as.numeric(R_train %*% mad_weights)

markowitz_sd = sd(markowitz_portfolio_returns)
mad_sd = sd(mad_portfolio_returns)

markowitz_mad = mean(abs(markowitz_portfolio_returns -
                           markowitz_expected_return))

mad_portfolio_mad = mean(abs(mad_portfolio_returns -
                               mad_expected_return))



# PORTFOLIO WEIGHTS


weight_comparison = data.frame(
  Asset = assets,
  Markowitz = as.numeric(m_weights),
  MAD = as.numeric(mad_weights)
) %>%
  mutate(Difference = MAD - Markowitz,
         Absolute_Difference = abs(Difference),
         Max_Weight = pmax(Markowitz, MAD))



# PORTFOLIO DISTANCE


L1_distance = sum(abs(m_weights - mad_weights))
minimum_weight_reallocation = L1_distance / 2



# PORTFOLIO WEIGHTS PLOT


weight_long = weight_comparison %>%
  arrange(desc(Max_Weight)) %>%
  mutate(Asset = factor(Asset, levels = rev(Asset))) %>%
  select(Asset, Markowitz, MAD) %>%
  pivot_longer(cols = c(Markowitz, MAD),
               names_to = "Model",
               values_to = "Weight")

p3 = ggplot(weight_long, aes(x = Asset, y = Weight, fill = Model)) +
  geom_col(position = "dodge", width = 0.75) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Portfolio Weights at the Common Target Return",
       x = "Asset", y = "Portfolio Weight", fill = "Model") +
  theme_minimal()

print(p3)



# PORTFOLIO CONCENTRATION


calculate_concentration = function(weights, model) {
  
  positive = weights[weights > tol]
  HHI = sum(weights^2)
  
  data.frame(
    Model = model,
    Number_of_Holdings = length(positive),
    Maximum_Weight = max(weights),
    HHI = HHI,
    Effective_Number_of_Assets = 1 / HHI
  )
}

concentration_results = bind_rows(
  calculate_concentration(m_weights, "Markowitz"),
  calculate_concentration(mad_weights, "MAD")
)

# ANNUALISED GEOMETRIC RETURN
annualised_geometric_return = function(x) {
  prod(1 + x)^(52 / length(x)) - 1
}

# TRAINING PERFORMANCE RESULTS

performance_results = data.frame(
  Model = c("Markowitz", "MAD"),
  Mean_Weekly_Return = c(mean(markowitz_portfolio_returns),
                         mean(mad_portfolio_returns)),
  Weekly_SD = c(markowitz_sd, mad_sd),
  Annualised_Return_Geometric = c(
    annualised_geometric_return(markowitz_portfolio_returns),
    annualised_geometric_return(mad_portfolio_returns)
  ),
  Annualised_Volatility = c(markowitz_sd * sqrt(52),
                            mad_sd * sqrt(52)),
  Weekly_MAD = c(markowitz_mad, mad_portfolio_mad)
)

# MAXIMUM DRAWDOWN

calculate_drawdown = function(returns) {
  
  if (any(1 + returns <= 0))
    stop("Returns contain values <= -100%.")
  
  wealth = cumprod(1 + returns)
  running_max = cummax(wealth)
  
  wealth / running_max - 1
}

markowitz_drawdown = calculate_drawdown(markowitz_portfolio_returns)
mad_drawdown = calculate_drawdown(mad_portfolio_returns)

drawdown_results = data.frame(
  Model = c("Markowitz", "MAD"),
  Maximum_Drawdown = c(min(markowitz_drawdown),
                       min(mad_drawdown))
)

# OPTIMISATION PROBLEM DIMENSIONS

markowitz_variables = n_assets
markowitz_constraints = 2 + n_assets + n_assets

mad_variables = n_assets + T_train
mad_constraints = 2 + (2 * T_train) + n_assets + T_train

problem_dimensions = data.frame(
  Model = c("Markowitz", "MAD"),
  Variables = c(markowitz_variables, mad_variables),
  Constraints = c(markowitz_constraints, mad_constraints)
)

# COMPUTATIONAL PERFORMANCE

timing_repetitions = 10
timing_list = list()

for (j in seq_along(target_grid)) {
  
  target = target_grid[j]
  markowitz_times = numeric(timing_repetitions)
  mad_times = numeric(timing_repetitions)
  
  for (k in seq_len(timing_repetitions)) {
    
    markowitz_times[k] =
      system.time(solve_markowitz(target, mu, Sigma,
                                  maximum_weight))[["elapsed"]]
    
    mad_times[k] =
      system.time(solve_mad(target, R_train, mu,
                            maximum_weight))[["elapsed"]]
  }
  
  timing_list[[j]] = data.frame(
    Target = target,
    Markowitz_Mean_Time = mean(markowitz_times),
    Markowitz_Median_Time = median(markowitz_times),
    MAD_Mean_Time = mean(mad_times),
    MAD_Median_Time = median(mad_times)
  )
}

timing_results = bind_rows(timing_list)

mean_markowitz_time = mean(timing_results$Markowitz_Mean_Time)
mean_mad_time = mean(timing_results$MAD_Mean_Time)
time_ratio = mean_mad_time / mean_markowitz_time


# MACHINE INFORMATION

get_powershell = function(command) {
  tryCatch(
    system2("powershell",
            c("-NoProfile", "-Command", command),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
}

machine_information = list(
  Processor = get_powershell(
    "(Get-CimInstance Win32_Processor | Select-Object -ExpandProperty Name)"
  ),
  Cores = get_powershell(
    "(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfCores -Sum).Sum"
  ),
  Logical_Processors = get_powershell(
    "(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfLogicalProcessors -Sum).Sum"
  ),
  RAM_GB = round(
    as.numeric(get_powershell(
      "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"
    )) / 1024^3, 2
  )
)

print(machine_information)

# FINAL CONSTRAINT CHECKS

markowitz_constraints_ok = all(
  abs(sum(m_weights) - 1) <= tol,
  abs(sum(m_weights * mu) - common_target) <= 1e-6,
  min(m_weights) >= -tol,
  max(m_weights) <= maximum_weight + tol
)

mad_constraints_ok = all(
  abs(sum(mad_weights) - 1) <= tol,
  abs(sum(mad_weights * mu) - common_target) <= 1e-6,
  min(mad_weights) >= -tol,
  max(mad_weights) <= maximum_weight + tol
)

# OUT-OF-SAMPLE TESTING: 2025

cat("\n============================================\n")
cat("OUT-OF-SAMPLE TESTING: 2025\n")
cat("============================================\n")

# Apply TRAINING weights to TESTING returns
markowitz_test_returns = as.numeric(R_test %*% m_weights)
mad_test_returns = as.numeric(R_test %*% mad_weights)

# Realised performance
markowitz_test_mean = mean(markowitz_test_returns)
mad_test_mean = mean(mad_test_returns)

markowitz_test_sd = sd(markowitz_test_returns)
mad_test_sd = sd(mad_test_returns)

markowitz_test_mad =
  mean(abs(markowitz_test_returns - markowitz_test_mean))

mad_test_mad =
  mean(abs(mad_test_returns - mad_test_mean))

markowitz_test_annual_return =
  prod(1 + markowitz_test_returns)^(52 / length(markowitz_test_returns)) - 1

mad_test_annual_return =
  prod(1 + mad_test_returns)^(52 / length(mad_test_returns)) - 1

markowitz_test_annual_volatility = markowitz_test_sd * sqrt(52)
mad_test_annual_volatility = mad_test_sd * sqrt(52)



# OUT-OF-SAMPLE DRAWDOWN

markowitz_test_drawdown = calculate_drawdown(markowitz_test_returns)
mad_test_drawdown = calculate_drawdown(mad_test_returns)



# OUT-OF-SAMPLE RESULTS


out_of_sample_results = data.frame(
  Model = c("Markowitz", "MAD"),
  Mean_Weekly_Return = c(markowitz_test_mean, mad_test_mean),
  Weekly_SD = c(markowitz_test_sd, mad_test_sd),
  Annualised_Return = c(markowitz_test_annual_return,
                        mad_test_annual_return),
  Annualised_Volatility = c(markowitz_test_annual_volatility,
                            mad_test_annual_volatility),
  Weekly_MAD = c(markowitz_test_mad, mad_test_mad),
  Maximum_Drawdown = c(min(markowitz_test_drawdown),
                       min(mad_test_drawdown))
)



# OUT-OF-SAMPLE WEALTH


markowitz_test_wealth = cumprod(1 + markowitz_test_returns)
mad_test_wealth = cumprod(1 + mad_test_returns)

wealth_comparison = data.frame(
  Week = seq_along(markowitz_test_returns),
  Markowitz = markowitz_test_wealth,
  MAD = mad_test_wealth
)

wealth_long = wealth_comparison %>%
  pivot_longer(cols = c(Markowitz, MAD),
               names_to = "Model",
               values_to = "Wealth")

p4 = ggplot(wealth_long, aes(x = Week, y = Wealth, colour = Model)) +
  geom_line(linewidth = 1.1) +
  labs(title = "Out-of-Sample Portfolio Performance: 2025",
       x = "Week", y = "Growth of €1", colour = "Model") +
  theme_minimal()

print(p4)



# OUT-OF-SAMPLE COMPARISON


test_sd_difference = mad_test_sd - markowitz_test_sd
test_sd_ratio = mad_test_sd / markowitz_test_sd

test_return_difference =
  mad_test_annual_return - markowitz_test_annual_return

test_mad_difference = mad_test_mad - markowitz_test_mad

test_drawdown_difference =
  min(mad_test_drawdown) - min(markowitz_test_drawdown)

out_of_sample_comparison = data.frame(
  Metric = c("Annualised Return",
             "Annualised Volatility",
             "Weekly MAD",
             "Maximum Drawdown"),
  Markowitz = c(markowitz_test_annual_return,
                markowitz_test_annual_volatility,
                markowitz_test_mad,
                min(markowitz_test_drawdown)),
  MAD = c(mad_test_annual_return,
          mad_test_annual_volatility,
          mad_test_mad,
          min(mad_test_drawdown))
)





# SESSION INFORMATION


sessionInfo()



# FINAL RESULTS


cat("\n============================================\n")
cat("SAMPLE INFORMATION\n")
cat("============================================\n")
cat("Number of assets:", n_assets, "\n")
cat("Training observations:", T_train, "\n")
cat("Testing observations:", T_test, "\n")
cat("Minimum feasible weekly return:", minimum_feasible_return, "\n")
cat("Maximum feasible weekly return:", maximum_feasible_return, "\n")
cat("GMV weekly return:", gmv_return, "\n")
cat("GMV weekly SD:", gmv_sd, "\n")

cat("\n============================================\n")
cat("TRAINING PERIOD ASSET STATISTICS\n")
cat("============================================\n")
print(asset_statistics)

cat("\n============================================\n")
cat("HIGHEST CORRELATIONS\n")
cat("============================================\n")
print(head(cor_high, 10))

cat("\n============================================\n")
cat("FRONTIER RESULTS\n")
cat("============================================\n")
print(frontier_results)

cat("\n============================================\n")
cat("FRONTIER COMPARISON\n")
cat("============================================\n")
print(frontier_comparison)

cat("\n============================================\n")
cat("MAD OBJECTIVE CHECK\n")
cat("============================================\n")
print(mad_objective_check)

cat("\n============================================\n")
cat("MARKOWITZ OBJECTIVE CHECK\n")
cat("============================================\n")
print(markowitz_objective_check)

cat("\n============================================\n")
cat("COMMON TARGET\n")
cat("============================================\n")
print(common_target)

cat("\n============================================\n")
cat("PORTFOLIO WEIGHTS\n")
cat("============================================\n")
print(weight_comparison)

cat("\n============================================\n")
cat("CONCENTRATION\n")
cat("============================================\n")
print(concentration_results)

cat("\n============================================\n")
cat("TRAINING PERIOD PERFORMANCE\n")
cat("============================================\n")
print(performance_results)

cat("\n============================================\n")
cat("TRAINING PERIOD DRAWDOWN\n")
cat("============================================\n")
print(drawdown_results)

cat("\n============================================\n")
cat("PROBLEM DIMENSIONS\n")
cat("============================================\n")
print(problem_dimensions)

cat("\n============================================\n")
cat("COMPUTATIONAL TIME\n")
cat("============================================\n")
print(timing_results)

cat("\n============================================\n")
cat("OUT-OF-SAMPLE 2025 RESULTS\n")
cat("============================================\n")
print(out_of_sample_results)

cat("\n============================================\n")
cat("OUT-OF-SAMPLE COMPARISON\n")
cat("============================================\n")
print(out_of_sample_comparison)

cat("\n============================================\n")
cat("FINAL SUMMARY\n")
cat("============================================\n")

cat("Mean Markowitz training SD:", mean_markowitz_frontier_sd, "\n")
cat("Mean MAD training SD:", mean_mad_frontier_sd, "\n")
cat("Mean training SD difference:", mean_sd_difference, "\n")
cat("Mean training SD ratio:", mean_sd_ratio, "\n")
cat("Training percentage SD difference:",
    percentage_sd_difference, "%\n")

cat("Mean Markowitz optimisation time:",
    mean_markowitz_time, "\n")
cat("Mean MAD optimisation time:",
    mean_mad_time, "\n")
cat("Time ratio:", time_ratio, "\n")

cat("L1 distance:", L1_distance, "\n")
cat("Minimum weight reallocation:",
    minimum_weight_reallocation, "\n")

cat("Markowitz constraints satisfied:",
    markowitz_constraints_ok, "\n")
cat("MAD constraints satisfied:",
    mad_constraints_ok, "\n")

cat("\n============================================\n")
cat("2025 OUT-OF-SAMPLE SUMMARY\n")
cat("============================================\n")

cat("Markowitz 2025 annualised return:",
    markowitz_test_annual_return, "\n")
cat("MAD 2025 annualised return:",
    mad_test_annual_return, "\n")

cat("Markowitz 2025 annualised volatility:",
    markowitz_test_annual_volatility, "\n")
cat("MAD 2025 annualised volatility:",
    mad_test_annual_volatility, "\n")

cat("Markowitz 2025 MAD:",
    markowitz_test_mad, "\n")
cat("MAD 2025 MAD:",
    mad_test_mad, "\n")

cat("Markowitz 2025 maximum drawdown:",
    min(markowitz_test_drawdown), "\n")
cat("MAD 2025 maximum drawdown:",
    min(mad_test_drawdown), "\n")

cat("2025 return difference (MAD - Markowitz):",
    test_return_difference, "\n")

cat("2025 SD difference (MAD - Markowitz):",
    test_sd_difference, "\n")

cat("2025 MAD difference (MAD - Markowitz):",
    test_mad_difference, "\n")

cat("2025 drawdown difference (MAD - Markowitz):",
    test_drawdown_difference, "\n")

# =============================================================================
# SUPPLEMENTARY COMPUTATIONAL SCALING EXPERIMENT
# FIXED 2021--2024 SAMPLE, INCREASING NUMBER OF ASSETS
# =============================================================================



# 1. Fixed 100-stock universe


# The first 20 assets are exactly the same as the main empirical analysis.
# Each larger universe contains the preceding smaller universe.

assets_100 <- c(
  
  # Original 20-asset universe
  "AAPL","MSFT","AMZN","GOOGL","META","NVDA",
  "JPM","BAC","JNJ","PFE","KO","PEP","WMT",
  "XOM","CVX","CAT","DIS","NFLX","V","MA",
  
  # Additional 80 assets
  "AVGO","BRK-B","LLY","COST","PG","HD","ORCL","ABBV","CRM","MRK",
  "AMD","TMO","CSCO","ACN","MCD","LIN","IBM","WFC","ABT","GE",
  "INTU","QCOM","VZ","CMCSA","AMGN","TXN","PM","DHR","ISRG","NOW",
  "SPGI","GS","UNH","AXP","AMAT","LOW","RTX","BKNG","HON","COP",
  "BLK","SYK","MDT","TJX","LMT","ADP","GILD","DE","SCHW","C",
  "UPS","CB","BA","MO","USB","TGT","MMM","F","GM","DUK",
  "SO","NEE","CL","BDX","CI","CVS","MDLZ","ADI","MU","LRCX",
  "PLD","CME","ICE","ETN","AON","APD","EOG","SLB","GD","NOC"
)

stopifnot(
  length(assets_100) == 100,
  length(unique(assets_100)) == 100
)

# Confirm that the first 20 assets match the main study
stopifnot(
  identical(
    assets_100[1:20],
    assets
  )
)



# 2. Scaling design


n_values <- c(20, 40, 60, 80, 100)

# Five target-return levels for each asset-universe size
n_targets_scaling <- 5

# Batch repetitions for stable timing measurement
timing_repetitions_scaling <- 50

# Download a sufficiently wide daily-price window so that weekly
# aggregation is performed consistently with the main analysis.
# The scaling SAMPLE itself is subsequently restricted to the
# same 207 weekly training dates from 2021--2024.
scaling_start <- as.Date("2021-01-01")
scaling_end   <- as.Date("2025-01-10")



# 3. Download daily adjusted prices


prices_scaling_list <- lapply(
  assets_100,
  function(ticker) {
    
    x <- getSymbols(
      ticker,
      src = "yahoo",
      from = scaling_start,
      to = scaling_end,
      auto.assign = FALSE
    )
    
    Ad(x)
  }
)

names(prices_scaling_list) <- assets_100



# 4. Check downloads


failed_logical <- vapply(
  prices_scaling_list,
  is.null,
  logical(1)
)

failed_assets <- names(
  prices_scaling_list
)[failed_logical]

cat("\nDownload summary\n")
cat("----------------\n")
cat("Requested assets:", length(assets_100), "\n")
cat("Successful downloads:", sum(!failed_logical), "\n")
cat("Failed downloads:", sum(failed_logical), "\n")

if (length(failed_assets) > 0) {
  
  cat("\nFailed tickers:\n")
  print(failed_assets)
  
  stop(
    "At least one asset failed to download. ",
    "Replace the failed ticker(s) and rerun the scaling experiment."
  )
}


# 5. Merge daily adjusted prices
prices_scaling_daily <- do.call(
  merge,
  prices_scaling_list
)

colnames(prices_scaling_daily) <- assets_100



# 6. Retain common daily observations


# This is important: all 100 assets must have prices on every retained date.
prices_scaling_daily <- na.omit(
  prices_scaling_daily
)

cat("\nCommon daily observations:",
    nrow(prices_scaling_daily), "\n")



# 7. Construct weekly prices and weekly returns


prices_scaling_weekly <- na.omit(
  apply.weekly(
    prices_scaling_daily,
    last
  )
)

returns_scaling_weekly <- na.omit(
  prices_scaling_weekly /
    lag(prices_scaling_weekly) - 1
)

# EXACT weekly dates used in the main training analysis
main_training_dates <- index(
  returns_weekly[training_dates]
)

# Restrict scaling sample to those same weekly dates
common_scaling_dates <- intersect(
  index(returns_scaling_weekly),
  main_training_dates
)

returns_scaling_train <- returns_scaling_weekly[
  common_scaling_dates
]

R_scaling_full <- coredata(
  returns_scaling_train
)

colnames(R_scaling_full) <- assets_100




# 8. Define fixed T for the scaling experiment


T_scaling <- nrow(
  R_scaling_full
)

cat("Main training observations:", T_train, "\n")
cat("Scaling observations:", T_scaling, "\n")

stopifnot(
  T_scaling == T_train
)



# 9. Run computational scaling experiment


scaling_results_list <- list()

for (n_current in n_values) {
  
  cat(
    "\nRunning n =",
    n_current,
    "with T =",
    T_scaling,
    "\n"
  )
  
  
  # ---------------------------------------------------------------------------
  # 9.1 Nested asset universe
  # ---------------------------------------------------------------------------
  
  selected_assets <- assets_100[
    seq_len(n_current)
  ]
  
  R_n <- R_scaling_full[
    ,
    selected_assets,
    drop = FALSE
  ]
  
  # T must remain constant across every value of n
  stopifnot(
    nrow(R_n) == T_scaling,
    ncol(R_n) == n_current
  )
  
  
  # ---------------------------------------------------------------------------
  # 9.2 Estimate model inputs
  # ---------------------------------------------------------------------------
  
  mu_n <- colMeans(
    R_n
  )
  
  Sigma_n <- cov(
    R_n
  )
  
  names(mu_n) <- selected_assets
  
  validate_covariance(
    Sigma_n
  )
  
  
  # ---------------------------------------------------------------------------
  # 9.3 Global minimum-variance portfolio
  # ---------------------------------------------------------------------------
  
  gmv_weights_n <- solve_gmv(
    mu = mu_n,
    Sigma = Sigma_n,
    max_weight = maximum_weight,
    tol = tol
  )
  
  if (is.null(gmv_weights_n)) {
    
    warning(
      "GMV optimization failed for n = ",
      n_current
    )
    
    next
  }
  
  gmv_return_n <- sum(
    gmv_weights_n * mu_n
  )
  
  
  # ---------------------------------------------------------------------------
  # 9.4 Maximum feasible return
  # ---------------------------------------------------------------------------
  
  maximum_return_n <- calculate_extreme_return(
    sort(
      mu_n,
      decreasing = TRUE
    ),
    maximum_weight
  )
  
  lower_target_n <- gmv_return_n + 1e-6
  upper_target_n <- maximum_return_n - 1e-6
  
  if (lower_target_n >= upper_target_n) {
    
    warning(
      "Invalid target-return range for n = ",
      n_current
    )
    
    next
  }
  
  
  # ---------------------------------------------------------------------------
  # 9.5 Target-return grid
  # ---------------------------------------------------------------------------
  
  scaling_target_grid <- seq(
    lower_target_n,
    upper_target_n,
    length.out = n_targets_scaling
  )
  
  
  # ---------------------------------------------------------------------------
  # 9.6 Time both optimization models
  # ---------------------------------------------------------------------------
  
  target_timing_list <- list()
  
  for (j in seq_along(scaling_target_grid)) {
    
    target <- scaling_target_grid[j]
    
    solution_m <- NULL
    solution_mad <- NULL
    
    
    # -------------------------------------------------------------------------
    # Markowitz batch timing
    # -------------------------------------------------------------------------
    
    markowitz_total <- system.time({
      
      for (k in seq_len(timing_repetitions_scaling)) {
        
        solution_m <- solve_markowitz(
          target,
          mu_n,
          Sigma_n,
          maximum_weight
        )
      }
      
    })[["elapsed"]]
    
    
    # -------------------------------------------------------------------------
    # MAD batch timing
    # -------------------------------------------------------------------------
    
    mad_total <- system.time({
      
      for (k in seq_len(timing_repetitions_scaling)) {
        
        solution_mad <- solve_mad(
          target,
          R_n,
          mu_n,
          maximum_weight
        )
      }
      
    })[["elapsed"]]
    
    
    # Mean execution time per optimization
    markowitz_mean_time <-
      markowitz_total /
      timing_repetitions_scaling
    
    mad_mean_time <-
      mad_total /
      timing_repetitions_scaling
    
    
    # -------------------------------------------------------------------------
    # Check solver success
    # -------------------------------------------------------------------------
    
    if (
      is.null(solution_m) ||
      is.null(solution_mad)
    ) {
      
      warning(
        "Solver failure for n = ",
        n_current,
        ", target = ",
        target
      )
      
      next
    }
    
    
    # -------------------------------------------------------------------------
    # Store target-level timing result
    # -------------------------------------------------------------------------
    
    target_timing_list[[length(target_timing_list) + 1]] <-
      data.frame(
        n = n_current,
        T = T_scaling,
        Target = target,
        Markowitz_Mean_Time = markowitz_mean_time,
        MAD_Mean_Time = mad_mean_time
      )
  }
  
  
  scaling_results_list[[length(scaling_results_list) + 1]] <-
    bind_rows(
      target_timing_list
    )
}



# 10. Combine timing results


scaling_timing_results <- bind_rows(
  scaling_results_list
)

if (nrow(scaling_timing_results) == 0) {
  
  stop(
    "Scaling experiment produced no valid timing results."
  )
}



# 11. Summarize results by asset-universe size


scaling_summary <- scaling_timing_results |>
  group_by(
    n,
    T
  ) |>
  summarise(
    
    Markowitz_Mean_Time =
      mean(Markowitz_Mean_Time),
    
    MAD_Mean_Time =
      mean(MAD_Mean_Time),
    
    MAD_to_Markowitz_Ratio =
      MAD_Mean_Time /
      Markowitz_Mean_Time,
    
    .groups = "drop"
  )



# 12. Add optimization dimensions


scaling_summary <- scaling_summary |>
  mutate(
    
    Markowitz_Variables =
      n,
    
    MAD_Variables =
      n + T,
    
    Markowitz_Constraint_Conditions =
      2 + 2 * n,
    
    MAD_Constraint_Conditions =
      2 + 2 * T + n + T
  )



# 13. Display scaling results


cat("\n============================================\n")
cat("COMPUTATIONAL SCALING SUMMARY\n")
cat("============================================\n")

print(
  scaling_summary
)

# Reported dissertation scaling results


reported_scaling <- data.frame(
  n = c(20, 40, 60, 80, 100),
  Markowitz = c(
    0.000440,
    0.001600,
    0.003600,
    0.002680,
    0.010300
  ),
  MAD = c(
    0.0481,
    0.1640,
    0.1710,
    0.1180,
    0.3740
  )
)

reported_scaling$MAD_to_Markowitz_Ratio <-
  reported_scaling$MAD / reported_scaling$Markowitz



# 14. Plot reported dissertation optimization times


reported_scaling_plot_data <- reported_scaling |>
  select(
    n,
    Markowitz,
    MAD
  ) |>
  pivot_longer(
    cols = c(
      Markowitz,
      MAD
    ),
    names_to = "Model",
    values_to = "Mean_Time"
  )

p_scaling <- ggplot(
  reported_scaling_plot_data,
  aes(
    x = n,
    y = Mean_Time,
    colour = Model
  )
) +
  geom_line(
    linewidth = 1.1
  ) +
  geom_point(
    size = 2.5
  ) +
  scale_x_continuous(
    breaks = reported_scaling$n
  ) +
  labs(
    title = "Computational Scaling with Asset-Universe Size",
    x = "Number of Assets",
    y = "Mean Optimization Time (seconds)",
    colour = "Model"
  ) +
  theme_minimal()

print(
  p_scaling
)



# 15. Plot relative execution time


p_scaling_ratio <- ggplot(
  reported_scaling,
  aes(
    x = n,
    y = MAD_to_Markowitz_Ratio
  )
) +
  geom_line(
    linewidth = 1.1
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    breaks = reported_scaling$n
  ) +
  labs(
    title = "Relative Computational Time of MAD and Markowitz",
    x = "Number of Assets",
    y = "MAD / Markowitz Mean Time"
  ) +
  theme_minimal()

print(
  p_scaling_ratio
)



# 16. Save figures


dir.create(
  "figures",
  showWarnings = FALSE
)

ggsave(
  "figures/computational_scaling.pdf",
  plot = p_scaling,
  width = 8,
  height = 5.5
)

ggsave(
  "figures/computational_scaling_ratio.pdf",
  plot = p_scaling_ratio,
  width = 8,
  height = 5.5
)



# 17. Fit empirical power functions


fit_markowitz <- lm(
  log(Markowitz) ~ log(n),
  data = reported_scaling
)

fit_mad <- lm(
  log(MAD) ~ log(n),
  data = reported_scaling
)

summary(
  fit_markowitz
)

summary(
  fit_mad
)



# 18. Recover power-function coefficients


a_markowitz <- exp(
  coef(fit_markowitz)[1]
)

b_markowitz <- coef(
  fit_markowitz
)[2]

a_mad <- exp(
  coef(fit_mad)[1]
)

b_mad <- coef(
  fit_mad
)[2]

a_markowitz
b_markowitz
a_mad
b_mad

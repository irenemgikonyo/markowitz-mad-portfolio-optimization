# Markowitz and MAD Portfolio Optimization

This repository contains the R implementation used for an empirical comparison of the Markowitz mean-variance and Mean Absolute Deviation (MAD) portfolio optimization models.

The analysis compares the two approaches in terms of portfolio risk, portfolio composition, in-sample performance, out-of-sample performance, and computational performance.

## Data

The principal analysis uses 20 large-cap U.S. equities with adjusted price data obtained from Yahoo Finance through the `quantmod` package.

The sample period is divided into:

- Training period: 2021–2024
- Testing period: 2025

Daily adjusted prices are converted to weekly simple returns. The training sample contains 207 weekly observations and the out-of-sample testing period contains 53 weekly observations.

## Portfolio Models

Two portfolio optimization approaches are implemented:

### Markowitz Mean-Variance

The Markowitz model minimizes portfolio variance subject to:

- Full investment
- No short selling
- A maximum individual asset weight of 30%
- A specified target expected return

The quadratic programming problem is solved using `quadprog`.

### Mean Absolute Deviation

The MAD model minimizes mean absolute deviation from expected portfolio return under the same portfolio constraints.

The linear programming formulation introduces auxiliary deviation variables and is solved using `lpSolve`.

## Empirical Comparison

The models are evaluated over a common grid of 50 feasible target returns.

A representative portfolio is selected from the 25th point of this grid for detailed comparison of:

- Portfolio weights
- Portfolio concentration
- In-sample risk and return
- Out-of-sample performance during 2025

The 2025 analysis applies the portfolio weights estimated from the 2021–2024 training period without re-optimization.

## Computational Scaling

A supplementary experiment evaluates computational performance for asset universes of:

`n = 20, 40, 60, 80, 100`

The number of weekly observations is held fixed at `T = 207`.

For each asset-universe size, five feasible target returns are evaluated and each optimization model is solved 50 times per target.

The empirical timing observations reported in the dissertation are stored in:

`results/reported_scaling_timings.csv`

Execution times are machine- and run-dependent. Re-running the scaling experiment may therefore produce different absolute timings. The stored results preserve the observations used in the dissertation, while the R script provides the complete experimental procedure.

The fitted power functions in the analysis are descriptive empirical approximations and should not be interpreted as theoretical computational-complexity functions.

## Repository Structure

```text
markowitz-mad-portfolio-optimization/
├── README.md
├── .gitignore
├── portfolio_optimization_analysis.R
├── figures/
│   ├── risk_return_frontier.pdf
│   ├── computational_scaling.pdf
│   ├── computational_scaling_ratio.pdf
│   └── out_of_sample_wealth.pdf
└── results/
    └── reported_scaling_timings.csv

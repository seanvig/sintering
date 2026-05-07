# Sintering

An interactive R Shiny application for modeling the kinetics of rapid alumina sintering, based on the 1974 doctoral thesis of Dennis H. Taylor.

**[View the published app](https://seanvig.shinyapps.io/sintering/)**

## Background

Sintering is the process of forming a solid mass from ceramic powder by heating it below the melting point — atoms diffuse and fuse particles together without liquefaction. This app modernizes Dennis Taylor's original research on rapid sintering of alumina (Al₂O₃), examining how trace dopants (Y₂O₃, SnO₂) affect sintering kinetics across three alumina particle morphologies (Linde A/B, Cera, and Fused).

The project was rebuilt by Sean Taylor (2017–2021) to translate his father's Fortran-based dynamic programming algorithm into interactive R code and to make the original research accessible as a living, exploratory tool.

## Features

- **Qualitative Comparisons** — side-by-side plots of sintering curves across dopant types and alumina morphologies
- **Kinetics Modeling** — interactive curve fitting using a two-phase composite model:
  - Non-isothermal phase: exponential growth `Y = A·eᴷᵗ`
  - Isothermal phase: Johnson's equation `Y = 1 − B·e^(−Jt^n)`
- **Dynamic Programming Optimization** — automatically finds the inflection point that minimizes combined sum of squared residuals across both fitted curves
- **Manual Parameter Tuning** — adjust A, K, B, J coefficients interactively and observe the effect in real time
- **Vignettes** — embedded mathematical derivations, historical context, and the original 1974 Fortran code

## Data

Experimental measurements were digitized manually from figures in Dennis Taylor's published thesis. Measurements are normalized to dimensionless ratios:

- `y = ΔL / ΔL_m` (normalized shrinkage)
- `t = (t − t₀) / t₀` (normalized time)

The dataset covers multiple alumina samples across dopant concentrations and particle morphologies.

## Installation

```r
# Install dependencies
install.packages(c("shiny", "dplyr", "ggplot2", "knitr", "rmarkdown"))

# Clone the repo and open in RStudio, then run:
shiny::runApp("app.R")
```

## Usage

Launch the app and navigate the four tabs:

| Tab | Contents |
|-----|----------|
| Introduction | Background on alumina, sintering physics, and the original thesis |
| Qualitative Comparisons | Plots comparing sintering rates across dopants and alumina types |
| Modeling | Interactive curve fitting with parameter controls and R² display |
| More | Appendices, Fortran source code, R equivalents, gallery, and background story |

## Tech Stack

- **R / Shiny** — application framework
- **ggplot2** — scientific visualization
- **dplyr** — data manipulation
- **knitr / R Markdown** — embedded vignettes and mathematical notation
- **nls / lm** — nonlinear and linear least-squares fitting
- **Fortran** (historical) — original 1974 algorithm, preserved in the appendix

## Authors

- **Sean Taylor** — R implementation and application design ([seanvig@gmail.com](mailto:seanvig@gmail.com))
- **Dennis H. Taylor** — original research and 1974 doctoral thesis

## License

MIT

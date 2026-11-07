# IEAP RStudio Series 04: Linear and Non-linear Regressions

This repository contains the code, data pipelines, and final report for the IEAP RStudio Series 04 assignment.

## Project Overview

The project is divided into two main analytical parts:

1. **How does random noise affect linear regression?**
   We simulate three synthetic datasets with identical underlying linear relationships but different levels of random noise. We fit linear models to each and investigate the effect of noise on parameter estimation, predictive performance, and explained variance. We also propose a Monte Carlo simulation protocol to better match theoretical expectations.

2. **Finding a mathematical model for heart rate recovery after exercise**
   Using two publicly available physiological recordings (a 7.5-minute recording and a 3-hour recording), we investigate the mathematical models that describe heart rate recovery (HRR). We test linear, mono-exponential, and bi-exponential decay models, supported by a literature review and residual diagnostic analysis.

## Repository Structure

- `IEAP-Series04-Rstudio.qmd`: The master Quarto document that compiles the full report.
- `sections/`: Contains the sub-documents for each part of the assignment.
  - `01-noise.qmd`: Part A (Noise & Linear Regression).
  - `02-hrr-literature.qmd`: Part B (HRR Literature Review).
  - `03-hrr-short.qmd`: Part B (7.5-minute HR recording analysis).
  - `04-hrr-long.qmd`: Part B (3-hour HR recording analysis).
  - `05-git-workflow.qmd`: Documentation of the team's Git strategy.
- `data/`: Contains the download scripts and raw data (ignored in Git).
- `references.bib`: The bibliography file.
- `apa.csl`: The citation style formatting file.

## Team Workflow

This project was developed collaboratively using Git. We used feature branches to work in parallel on different sections. Work was merged back into the `main` branch via Pull Requests once complete. For a detailed breakdown of this workflow, please refer to `sections/05-git-workflow.qmd`.

## How to Render the Report

To render the full HTML report, open `IEAP-Series04-Rstudio.qmd` in RStudio and click **Render**, or run:
```bash
quarto render IEAP-Series04-Rstudio.qmd
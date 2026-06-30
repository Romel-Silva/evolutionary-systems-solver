# Evolutionary Solver for Systems of Equations

This repository contains the GNU Octave source code accompanying the paper:

**"A Novel Transformation Technique for Solving Highly Linear Systems of Equations via Evolutionary Algorithms."**

## Overview

The repository provides the implementations used in the numerical experiments reported in the paper. The codes illustrate how the proposed objective function ($F_{ob}$) combines deterministic linear-system techniques with a genetic algorithm to solve systems containing linear and nonlinear equations.

The examples include test problems presented in the manuscript, such as combustion, neurophysiology, and arithmetic systems.

## Repository Structure

* `/Combustion` – Combustion equilibrium problem.
* `/Neurophysiology` – Neurophysiology benchmark.
* `/Arithmetic` – Arithmetic benchmark.

*Note:* Additional scripts comparing the proposed objective function ($F_{ob}$) with the traditional formulation ($F_{ma}$) are included.

## Reproducibility

The repository contains the source code used to generate the numerical results reported in the paper. Algorithmic parameters, including population size, crossover, mutation, elitism, and stopping criteria, are provided directly in the scripts.

## Requirements

- [GNU Octave](https://octave.org/)
- `statistics` package (for the Wilcoxon rank-sum test).
- `ga` package (Octave-Forge), only in some scripts.

## Running the Code

1. Open the desired example directory.
2. Run the main script.

*Note:* The comparison scripts require the `ga` package to be installed and loaded.

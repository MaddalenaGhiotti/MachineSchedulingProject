# Machine Scheduling Project

This repository contains the implementation and analysis of heuristic algorithms for job scheduling problems, developed for two academic assignments in Operations Research (Politecnico di Torino).

## Homework 2 – Single-machine Scheduling

### Problem Statement
Minimize the maximum lateness (Lmax) in a single-machine scheduling problem with sequence-dependent setup times.

### Implemented Heuristics and models
- **Constructive Heuristics**:
  - EDD (Earliest Due Date)
  - SetUp (Minimize setup times)
  - Prod (Shortest Processing Time)

- **Iterative Heuristics**:
  - latenessMax
  - setupMax
  - Randk

- **Hybrid Heuristics**: Combination of constructive and iterative methods, e.g., `EDD + latenessMax`, `EDD + (setupMax + latenessMax) + Randk`.

- **MILP Model**: Used to benchmark heuristic solutions for small problem sizes.

### Tools
Implemented in MATLAB.

### Key Insights
EDD is the most balanced constructive heuristic in terms of effectiveness and execution time. Hybrid approaches starting from EDD perform significantly better than iterative methods alone.

---

## Homework 3 – Parallel Machine Scheduling

### Problem Statement
Extend the problem to multiple parallel machines, maintaining sequence-dependent setup times and minimizing Lmax.

### Implemented Heuristics
- **Constructive Heuristic**:
  - EDD (extended to parallel machines)

- **Iterative Heuristics**:
  - latenessMax
  - setupMax
  - Rand(h, k)

- **Hybrid Heuristics**: Analogous to HW2 but adapted to handle parallel machine configurations.

### Tools
Implemented in MATLAB.

### Key Insights
Hybrid methods show superior performance when starting from a smart initial configuration (e.g., EDD).

---

## Performance Evaluation
All algorithms are evaluated on:
- **Effectiveness**: Lmax value compared to MILP or best-found solution.
- **Efficiency**: Execution time.

### Summary Tables
Average results are provided for each heuristic, including % improvement and runtime. See the report PDFs for detailed plots and data.

---

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.
This project is part of an academic coursework and is intended for educational purposes.

## Authors
- Maddalena Ghiotti
- Aldo Sambo

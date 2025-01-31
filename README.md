# Sampling from Density power divergence-based Generalized posterior distribution via Stochastic optimization

This repository contains the simulation code used in the paper:

<a href="https://arxiv.org/abs/2501.07790">Sonobe, N., Momozaki, T., Nakagawa, T.(2025). Sampling from Density power divergence-based Generalized posterior distribution via Stochastic optimization.</a>

Explanations for each component of the code are provided below.

- `Simulation_1_1.R` : Simulation study demonstrating that our stochastic gradient approach performs equivalently to the standard LLB algorithm (Section 4.1)
- `Simulation_1_2.R` : Simulation study assessing the accuracy of our proposed method in sampling from the Density Power Divergence (DPD)-based posterior (Section 4.1)
- `Simulation_Compare_MultiNormal.R` : Computational time comparison between the proposed method and conventional method (Section 4.2)
- `Simulation_Compare_MultiNormal_MSE.R` : Mean Squared Error (MSE) comparison between the proposed method and conventional method (Section 4.2)
- `Simulation_Compare_PoisR.R` : Comparative simulation study of the proposed method with the conventional method (Section 4.3)

# Sampling from Density power divergence-based Generalized posterior distribution via Stochastic optimization

This repository contains the simulation code used in the paper:

<a href="https://arxiv.org/abs/2501.07790">Sonobe, N., Momozaki, T., Nakagawa, T.(2025). Sampling from Density power divergence-based Generalized posterior distribution via Stochastic optimization.</a>

Explanations for each component of the code are provided below.

- `Simulation_1_1.R` : Simulation study demonstrating that our stochastic gradient approach performs equivalently to the standard LLB algorithm (Section 4.1)
- `Simulation_1_2.R` : Simulation study assessing the accuracy of our proposed method in sampling from the Density Power Divergence (DPD)-based posterior (Section 4.1)
- `DPD_GB_Normal.stan` : Stan code used in `Simulation_1_2.R` (Section 4.1)
- `Simulation_Compare_MultiNormal.R` : Computational time comparison between the proposed method and conventional method (Section 4.2)
- `Simulation_Compare_MultiNormal_MSE.R` : Mean Squared Error (MSE) and posterior variance comparison between the proposed method and conventional method (Section 4.2)
- `Simulation_Compare_PoisR.R` : Comparative simulation study of the proposed method with the conventional method (Section 4.3.1)
- `PoisR_medpar` : Real data analysis for the data about hospital length of stay (Section 4.3.2)
- `DPD_LLB_IG.R` : Parameter estimation for an Inverse Gaussian distribution (Supplementary Materials, Section B)
- `DPD_LLB_Gompertz.R` : Parameter estimation for an Gompertz distribution (Supplementary Materials, Section B)
- `Comparison_StandardBayes_DPDLLB.R` : Univariate contamination and comparison with Standard Bayes (Supplementary Materials, Section D)
- `Comparison_LLBSGD_MCMC.R` : Two-dimensional credible regions and DPD-based posterior sampling (Supplementary Materials, Section D)

# LDPC Codes in 5G NR Implementation in MATLAB

## Description

This repository provides a comprehensive implementation of Low Density Parity Check (LDPC) codes for 5G New Radio (NR) in MATLAB. It encompasses the entire transmission chain, from base graph expansion and signal generation through modulation, channel modeling, and soft-decision decoding. The project demonstrates the performance of various decoding algorithms under different channel conditions.

## Key Features

* **Base Graph & Parity-Check Matrix Generation:** Supports two different 5G NR base graphs (**NR_2_6_52** and **NR_1_5_352**) with a configurable expansion factor (lifting size).
* **Encoding**: Systematic encoding using the double-diagonal structure and back-substitution based on the expanded parity-check matrix.
* **Modulation and Channel**: BPSK modulation over an AWGN channel model with adjustable Eb/N0.
* **Soft-Decision Decoding**: Implementation of two soft  (Belief Propagation) decoding algorithms:
  * Min-Sum Algorithm
  * Sum-Product Algorithm
* **Performance Visualization**: Scripts to plot error-probability curves, BER vs. Eb/N0 in both linear and logarithmic scales, and convergence behavior across iterations.

## Results

* **Error Probability vs. Eb/N0**: Visualization of decoding error probability for code rates 1/4, 1/3, 1/2, and 3/5 under hard and soft decoding.
* **BER Performance in Log Scale**: Comparison of BER curves against theoretical Shannon limits and normalization approximations.
* **Success Probability vs. Iterations**: Convergence plots showing iteration count required for successful decoding at various Eb/N0 levels for both min-sum and sum-product algorithms.
* **Algorithm Comparison**: Comparative analysis highlighting improved BER performance and smoother decoding convergence with the sum-product algorithm over min-sum.

---

> **Note**: Detailed implementation notes and theoretical background are provided in the project report (`Group_24_CT216_Project_Report.pdf`). Please refer to it for in-depth explanations and derivations.

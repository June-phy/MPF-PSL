# Quantitative Macromolecular Proton Fraction Imaging using Pulsed Spin-Lock - Simulation Code

This repository contains MATLAB simulation code for the paper:

**"Quantitative Macromolecular Proton Fraction Imaging using Pulsed Spin-Lock"**  
Authors: Qianxue Shan, Ziqiang Yu, Baiyan Jiang, Jian Hou, Qiuyi Shen, Winnie Chiu-Wing Chu, Vincent Wai-Sun Wong, Weitian Chen

---

## Description

This repository provides the simulation framework for quantitative macromolecular proton fraction (MPF) imaging using pulsed spin-lock techniques (MPF-PSL).  
**The main scripts are:**
- `Simulation_study_1.m`: Demo script for Simulation Study 1 in the MPF-PSL paper.  
  This script investigates the sensitivity of the quantitative measurements \( R_{mpfsl} \) and \( R_{mpfsl,pul} \) to various tissue parameters, including \( R_{1a} \), \( R_{2a} \), \( T_{2b} \), \( K_{ca} \), and \( f_b \).
- `Simulation_study_2.m`: Demo script for Simulation Study 2 in the MPF-PSL paper.  
  This script investigates how the Relative Measurement Precision (RMP) changes with increasing spin-lock duration.

Other functions and scripts provide:
- Definition of tissue and sequence parameters for liver and other biological tissues.
- Numerical simulation of magnetization evolution based on the Bloch-McConnell equations under various spin-lock and adiabatic pulse conditions.
- Generation of adiabatic half-passage (AHP) and reverse AHP pulses.

The code is intended for academic and research use, facilitating the quantitative analysis of MPF imaging using advanced spin-lock techniques.

---

## Citation

If you use this code or results in your research, please cite our paper:

> Qianxue Shan, Ziqiang Yu, Weitian Chen. Quantitative Macromolecular Proton Fraction Imaging using Pulsed Spin-Lock. *[MRM]*

---

## License

This code is for academic and research purposes only. Commercial use or redistribution is strictly prohibited without explicit written permission from the authors.

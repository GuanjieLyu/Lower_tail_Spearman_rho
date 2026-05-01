# Bernstein Smoothing and Bootstrap Inference for Lower-Tail Spearman's Rho

This repository contains the R code used to reproduce the simulation study for the paper:

> **Bernstein smoothing and bootstrap inference for lower-tail Spearman's rho**  
> Guanjie Lyu, Frédéric Ouimet, and Selim Orhun Susam

The paper studies a Bernstein-smoothed estimator of the lower-tail Spearman's rho curve and compares it with the empirical copula-based estimator. The code in this repository reproduces the Monte Carlo experiments, figures, and tables used to assess finite-sample performance, bootstrap confidence intervals, threshold selection, and sensitivity to the Bernstein polynomial degree.

## Overview

For a bivariate copula \(C\), the lower-tail Spearman's rho curve is defined as

\[
\rho_L(p)
=
\frac{
\int_0^p\int_0^p C(u,v)\,du\,dv - p^4/4
}{
p^3/3-p^4/4
},
\qquad 0<p\leq 1.
\]

This functional summarizes concordance in the lower-left region \([0,p]^2\) of the copula. The repository implements and compares two estimators:

1. the empirical copula-based estimator;
2. the Bernstein-smoothed estimator.

The simulations investigate three main questions:

1. whether the smooth Bernstein bootstrap provides reliable pointwise confidence intervals and supports stable threshold selection;
2. how the Bernstein estimator compares with the empirical copula-based estimator in terms of ISE, MISE, integrated squared bias, and integrated variance;
3. whether the rule-of-thumb Bernstein degree \(m=\lfloor n^{2/3}\rfloor\) is reasonable in finite samples.


# BalExMed

BalExMed implements Bayesian mediation analysis for a **compositional exposure**,
such as the relative abundances of gut microbial taxa, a continuous mediator and a
continuous outcome. Compositional mediation methods have mostly treated the
microbiome as the mediator. BalExMed treats it as the exposure.

The composition enters through a latent balance: the normalized log-ratio between
the geometric means of two groups of taxa. The two groups, the direct effect and the
indirect effect are estimated jointly:

$$
m_i = \alpha B(z, x_i) + \psi_m^\top c_i + e_i, \qquad
y_i = \gamma B(z, x_i) + \beta m_i + \psi_y^\top c_i + \varepsilon_i ,
$$

where $z_j \in \{-1, 0, 1\}$ places taxon $j$ in the denominator, outside the
balance, or in the numerator. For a unit increase in the balance the direct effect
is $\gamma$ and the indirect effect is $\alpha\beta$. The sampler returns
posterior draws of the effects and a posterior inclusion probability for each taxon.

## Installation

```r
# from a local copy of this folder
install.packages("path/to/BalExMed", repos = NULL, type = "source")

# or from GitHub
remotes::install_github("drjingma/BalExMed")
```

BalExMed needs R 4.0 or later and the `coda` package.

## Example

```r
library(BalExMed)
set.seed(1)
sim <- simulate_balexmed(n = 150, d = 30)

fit <- balexmed(sim$X, sim$m, sim$y, covariates = sim$covariates,
                eta = 1/5, n_iter = 4000, burn_in = 2000,
                n_chains = 4, cores = 4, z_init = "pca", seed = 1)

summary(fit)            # direct, indirect and total effects; Rhat and ESS
taxa_table(fit)         # posterior inclusion probability of each taxon
plot(fit)               # trace plots of the effects
selected_z(fit)         # balance configuration at PIP > 0.5
```

See `vignette("BalExMed")` for a walk-through.

## Practical notes

- **Zeros.** The balance uses logarithms, so zeros must be replaced before fitting,
  for example with a pseudocount or a model-based imputation.
- **Prior sparsity.** `eta` is the prior probability that a taxon enters each side of
  the balance. `eta = 1/3` is non-informative; smaller values, such as 1/5 or 1/20,
  give sparser balances and usually better taxon selection when many taxa are
  measured.
- **Samplers.** The default Gibbs sampler updates every taxon in each iteration and
  mixes well. The Metropolis-Hastings sampler (`sampler = "mh"`) is cheaper per
  iteration and can help when the number of taxa is very large.
- **Covariates.** An intercept is always included. Data frames are expanded with
  `model.matrix()`, so factors are allowed.
- **Interpretation.** The effects are causal only under the standard assumptions of
  mediation analysis, including no unmeasured confounding of the
  mediator-outcome relationship.

## Citation

Please cite the accompanying paper, *Mediation Analysis with Compositional
Exposures*; `citation("BalExMed")` gives the current reference.

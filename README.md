# Uncomputable Fractal Dimensions and the Chaitin Barrier

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23194419.svg)](https://doi.org/10.5281/zenodo.23194419)

Code, data, figures and manuscript for

> **R. Chen**, *Uncomputable Fractal Dimensions, Kolmogorov Complexity Cascades, and the Chaitin Barrier in
> Physical Systems* (2026). DOI: [10.5281/zenodo.23194419](https://doi.org/10.5281/zenodo.23194419)

## Summary

- **A computable set with uncomputable dimension.** A homogeneous Cantor set with stage ratios
  r_k = 2^(−1/d_k), d_k = 1/4 + Ω_k/2 (computable lower approximations of Chaitin's Ω), is a computable compact
  set whose box-counting and Hausdorff dimensions both equal 1/4 + Ω/2, a left-c.e. but uncomputable real.
  Self-similar sets with finitely many computable ratios always have computable dimension.
- **Three tiers of Kolmogorov complexity under scaling** (L = ln 1/ε): Θ(log L) for every computable set (smooth
  manifolds, self-similar fractals, and the set above); Θ(L) with one incompressible symbol per scale;
  e^Θ(L) with independent randomness per cell. Questions about finite discretisations are decidable;
  undecidability enters only through limits.
- **No "Planck firewall" against Gödel.** ZFC-independence already occurs at a description of order 10⁵ bits
  (Yedidia–Aaronson, 7910 states), whereas a 1 cm³ crystal at atomic resolution carries up to 10²⁴ bits and a
  one-dimensional random fractal crosses the barrier at L* ≈ 34.5, far above the Planck scale (L_P ≈ 80 for 1 m).
- **Numerics.** Binary lambda calculus gives a rigorous lower bound Ω_BLC ≥ 0.12018 (closed-term counts match
  OEIS A114852); LZMA compression and box counting for four classes of Cantor-type sets.

This is the third paper of a series; see also
[Holographic Entanglement of Fractal Regions](https://doi.org/10.5281/zenodo.23189575)
([code](https://github.com/Ruqing1963/holographic-fractal-entanglement)) and
[Multifractal Large Deviations in Quantum Entanglement](https://doi.org/10.5281/zenodo.23193512)
([code](https://github.com/Ruqing1963/multifractal-quantum-entanglement)).

## Repository structure

| Path | Contents |
|---|---|
| `paper/` | LaTeX source and compiled PDF of the manuscript |
| `code/chaitin_fractal.py` | Single script producing every number and the figure in the paper |
| `figures/` | Figure as vector PDF (used by the paper) and PNG |
| `data/` | Numerical data as CSV (metadata in `#` header lines) |
| `results/` | Console output of the script (the numbers quoted in the paper) |

## Reproducing the results

Requirements: Python ≥ 3.10 and the packages in `requirements.txt`
(tested with Python 3.12.4, NumPy 1.26.4, Matplotlib 3.8.4). Runtime is about 10 seconds.

```bash
pip install -r requirements.txt
python code/chaitin_fractal.py            # options: --nmax 24 (max BLC program length), --show
```

Run from the repository root; the figure goes to `figures/` and data to `data/` (override with `CFC_FIG_DIR`,
`CFC_DATA_DIR`). The random constructions use a fixed seed, so the output is reproducible.

| Data file | Content | Paper |
|---|---|---|
| `blc_closed_terms.csv` | closed-term counts, normalised/undecided terms, lower bound on Ω_BLC per length | §5.1 |
| `cantor_selfsimilar.csv`, `cantor_omega_dim.csv`, `cantor_random_level.csv`, `cantor_random_node.csv` | per stage: L, number of intervals, LZMA bits above baseline, box dimension, theoretical dimension | §5.2, Fig. 1 |
| `chaitin_barrier.csv` | illustrative c_F, L*, ε*, Planck L_P | §4 |

To rebuild the paper (pdfLaTeX, two passes):

```bash
cd paper
pdflatex Chen_2026_Uncomputable_Fractal_Dimensions.tex
pdflatex Chen_2026_Uncomputable_Fractal_Dimensions.tex
```

## Citation

```bibtex
@misc{Chen2026UncomputableFractal,
  author = {Chen, Ruqing},
  title  = {Uncomputable Fractal Dimensions, Kolmogorov Complexity Cascades,
            and the Chaitin Barrier in Physical Systems},
  year   = {2026},
  doi    = {10.5281/zenodo.23194419},
  url    = {https://doi.org/10.5281/zenodo.23194419}
}
```

## License

- **Code** (`code/`): [MIT License](LICENSE)
- **Manuscript, figures, data and results** (`paper/`, `figures/`, `data/`, `results/`):
  [CC BY 4.0](LICENSE-CC-BY-4.0.md)

## Contact

Ruqing Chen — GUT Geoservice Inc., Montreal — ruqing@hotmail.com

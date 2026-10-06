"""
Uncomputable fractal dimensions, Kolmogorov-complexity growth laws and the Chaitin barrier.

Part A  Toy halting probability of binary lambda calculus (BLC, Tromp):
        programs = closed lambda terms in de Bruijn binary encoding (00 = lambda, 01 = application,
        1^{n}0 = variable n); the encoding is self-delimiting, so the program set is prefix-free.
        Omega_BLC = sum over closed terms with a beta-normal form of 2^{-|p|}  (a left-c.e. real).
        We compute a rigorous interval [lower, upper] from all terms of length <= n_max:
          lower = sum of 2^{-|p|} over terms shown to normalise within the step budget,
          upper = lower + (undecided terms) + (Kraft mass of all longer closed terms).
        Only the bits common to both ends of the interval are certified.
        (We do not claim Omega_BLC is Martin-Loef random; that needs an optimal universal machine.)

Part B  Four Cantor-type sets on [0,1] (each interval keeps its two end pieces of ratio r):
        (1) self-similar, r = 1/4                      : computable,      K(S_eps) = O(log L)
        (2) r_k = 2^{-1/d_k}, d_k = 1/4 + Omega_k/2     : computable,      K = O(log L),
            but dimension D = 1/4 + Omega/2 is not computable when Omega is Chaitin's Omega
        (3) one random bit per level selects r_k         : K = Theta(L)       (stand-in for Omega bits)
        (4) one random bit per interval selects r        : K = Theta(N_eps) = e^{Theta(L)}
        For each level we rasterise on 2^22 pixels and report the LZMA-compressed size
        (an upper bound on K, up to an additive constant) and box-counting dimension.

Usage:  python chaitin_fractal.py [--nmax 24] [--show]

Figures are written to ./figures (PNG and PDF), data to ./data (CSV); override with the environment
variables CFC_FIG_DIR and CFC_DATA_DIR.

Companion code for: R. Chen, "Uncomputable Fractal Dimensions, Kolmogorov Complexity Cascades, and the
Chaitin Barrier in Physical Systems" (2026), DOI: 10.5281/zenodo.23194419
"""
import os
import sys
import lzma
import argparse
import numpy as np
import matplotlib.pyplot as plt
from functools import lru_cache

sys.setrecursionlimit(20000)
FIG_DIR = os.environ.get("CFC_FIG_DIR", "figures")
DATA_DIR = os.environ.get("CFC_DATA_DIR", "data")


def savefig(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=150)


def savedata(name, columns, cols, meta=()):
    os.makedirs(DATA_DIR, exist_ok=True)
    arr = np.column_stack([np.asarray(c, dtype=float) for c in cols])
    with open(os.path.join(DATA_DIR, f"{name}.csv"), "w", encoding="utf-8", newline="\n") as fh:
        for line in meta:
            fh.write(f"# {line}\n")
        fh.write(",".join(columns) + "\n")
        np.savetxt(fh, arr, delimiter=",", fmt="%.10g")


# ================================================================ Part A: binary lambda calculus
# terms: ('v', i) with 0-based de Bruijn index i;  ('l', body);  ('a', f, x)
def size(t):
    if t[0] == 'v':
        return t[1] + 2
    if t[0] == 'l':
        return 2 + size(t[1])
    return 2 + size(t[1]) + size(t[2])


def nodes(t):
    if t[0] == 'v':
        return 1
    if t[0] == 'l':
        return 1 + nodes(t[1])
    return 1 + nodes(t[1]) + nodes(t[2])


@lru_cache(maxsize=None)
def gen(s, k):
    """All terms of encoded length exactly s whose free indices are < k."""
    out = []
    if s >= 2 and s - 2 < k:                 # variable with index s-2
        out.append(('v', s - 2))
    if s >= 4:
        out += [('l', b) for b in gen(s - 2, k + 1)]
        for s1 in range(2, s - 3):
            fs = gen(s1, k)
            if not fs:
                continue
            xs = gen(s - 2 - s1, k)
            out += [('a', f, x) for f in fs for x in xs]
    return tuple(out)


def shift(t, d, c=0):
    if t[0] == 'v':
        return ('v', t[1] + d) if t[1] >= c else t
    if t[0] == 'l':
        return ('l', shift(t[1], d, c + 1))
    return ('a', shift(t[1], d, c), shift(t[2], d, c))


def subst(t, j, s):
    if t[0] == 'v':
        return s if t[1] == j else t
    if t[0] == 'l':
        return ('l', subst(t[1], j + 1, shift(s, 1)))
    return ('a', subst(t[1], j, s), subst(t[2], j, s))


def step(t):
    """One leftmost-outermost beta step; returns (term, reduced?)."""
    if t[0] == 'a':
        f, x = t[1], t[2]
        if f[0] == 'l':
            return shift(subst(f[1], 0, shift(x, 1)), -1), True
        f2, r = step(f)
        if r:
            return ('a', f2, x), True
        x2, r = step(x)
        return ('a', f, x2), r
    if t[0] == 'l':
        b, r = step(t[1])
        return ('l', b), r
    return t, False


def normalises(t, max_steps=400, max_nodes=3000):
    """True if a normal form is reached, None if undecided within the budget."""
    for _ in range(max_steps):
        t, r = step(t)
        if not r:
            return True
        if nodes(t) > max_nodes:
            return None
    return None


def omega_blc(nmax):
    lower, undecided, kraft = 0.0, 0.0, 0.0
    rows = []
    for s in range(4, nmax + 1):
        terms = gen(s, 0)
        h = u = 0
        for t in terms:
            r = normalises(t)
            if r is True:
                h += 1
            else:
                u += 1
        lower += h * 2.0 ** -s
        undecided += u * 2.0 ** -s
        kraft += len(terms) * 2.0 ** -s
        rows.append((s, len(terms), h, u, lower))
    upper = lower + undecided + (1.0 - kraft)          # Kraft: total mass of closed terms <= 1
    return lower, upper, rows


def common_bits(a, b, nbits=30):
    """Number of leading binary digits shared by a and b (both in [0,1))."""
    k = 0
    for _ in range(nbits):
        a, b = 2 * a, 2 * b
        if int(a) != int(b):
            break
        a, b = a - int(a), b - int(b)
        k += 1
    return k


def bits(x, n):
    out = ""
    for _ in range(n):
        x *= 2
        out += str(int(x))
        x -= int(x)
    return out


# ================================================================ Part B: Cantor-type sets
G = 2 ** 22                                   # pixels on [0,1]


def build(kind, levels, omega_seq=None, rng=None, ra=0.25, rb=0.10):
    """Return list over levels of (intervals array [N,2], eps_k, D_k_theory)."""
    iv = np.array([[0.0, 1.0]])
    out = []
    logsum, ln2k = 0.0, 0.0
    for k in range(1, levels + 1):
        n = len(iv)
        if kind == "selfsimilar":
            r = np.full(n, 0.25)
        elif kind == "omega_dim":
            d = 0.25 + omega_seq[min(k - 1, len(omega_seq) - 1)] / 2
            r = np.full(n, 2.0 ** (-1.0 / d))
        elif kind == "random_level":
            r = np.full(n, ra if rng.integers(2) else rb)
        elif kind == "random_node":
            r = np.where(rng.integers(2, size=n) == 1, ra, rb)
        lengths = iv[:, 1] - iv[:, 0]
        left = np.column_stack([iv[:, 0], iv[:, 0] + r * lengths])
        right = np.column_stack([iv[:, 1] - r * lengths, iv[:, 1]])
        iv = np.empty((2 * n, 2)); iv[0::2] = left; iv[1::2] = right
        lens = iv[:, 1] - iv[:, 0]
        if lens.min() * G < 1:
            break
        logsum += np.mean(np.log(1 / r)); ln2k += np.log(2)
        out.append((iv.copy(), float(np.exp(np.mean(np.log(lens)))), ln2k / logsum))
    return out


def raster(iv):
    b = np.zeros(G, dtype=np.uint8)
    lo = np.floor(iv[:, 0] * G).astype(np.int64)
    hi = np.ceil(iv[:, 1] * G).astype(np.int64)
    d = np.zeros(G + 1, dtype=np.int64)
    np.add.at(d, lo, 1); np.add.at(d, hi, -1)
    b[np.cumsum(d)[:-1] > 0] = 1
    return b


def lzma_bits(b):
    return 8 * len(lzma.compress(np.packbits(b).tobytes(), preset=9 | lzma.PRESET_EXTREME))


def box_dim(b, eps):
    """Box-counting slope over scales between the pixel size and ~eps^(1/2)."""
    js, ns = [], []
    for j in range(2, 23):
        m = 2 ** j
        if m > G:
            break
        occ = b.reshape(m, -1).max(axis=1).sum()
        js.append(j); ns.append(occ)
    js, ns = np.array(js), np.array(ns)
    sel = (2.0 ** -js >= eps) & (js >= 3)
    if sel.sum() < 3:
        return np.nan
    return np.polyfit(js[sel] * np.log(2), np.log(ns[sel]), 1)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=24, help="max BLC program length")
    ap.add_argument("--show", action="store_true")
    a = ap.parse_args()

    # ---------------- Part A
    lower, upper, rows = omega_blc(a.nmax)
    print("Part A: toy halting probability of binary lambda calculus")
    print("   len   #closed  #normalised  #undecided   lower bound")
    for s, n, h, u, lo in rows:
        if n:
            print(f"   {s:3d}  {n:8d}  {h:10d}  {u:10d}   {lo:.10f}")
    nc = common_bits(lower, upper)
    print(f"   Omega_BLC in [{lower:.10f}, {upper:.10f}]")
    print(f"   lower = 0.{bits(lower, 24)}...  upper = 0.{bits(upper, 24)}...")
    print(f"   certified leading bits: {nc}")
    R = np.array(rows, dtype=float)
    savedata("blc_closed_terms", ["length", "n_closed", "n_normalised", "n_undecided", "omega_lower_bound"],
             [R[:, i] for i in range(5)],
             ["binary lambda calculus, closed terms; normal-order reduction, budget 400 steps / 3000 nodes",
              f"Omega_BLC in [{lower:.12f}, {upper:.12f}] (upper bound trivial, see paper); certified bits: {nc}",
              "n_closed agrees with OEIS A114852"])
    omega_seq = [row[4] for row in rows]                  # increasing computable approximations

    # ---------------- Part B
    rng = np.random.default_rng(1963)
    kinds = [("selfsimilar", "(1) self-similar r=1/4 (computable)"),
             ("omega_dim", "(2) r_k from Omega_k: D = 1/4 + Omega/2 (computable approx.)"),
             ("random_level", "(3) random bit per level (stand-in for Omega bits)"),
             ("random_node", "(4) random bit per interval")]
    results = {}
    base = lzma_bits(np.zeros(G, dtype=np.uint8))
    print(f"\nPart B: compressed size minus empty-bitmap baseline ({base} bits), and box dimension")
    for kind, label in kinds:
        levs = build(kind, 40, omega_seq=omega_seq, rng=rng)
        Ls, Ks, Dbox, Dth, N = [], [], [], [], []
        for iv, eps, dth in levs:
            b = raster(iv)
            Ls.append(np.log(1 / eps)); Ks.append(lzma_bits(b) - base); Dbox.append(box_dim(b, eps))
            Dth.append(dth); N.append(len(iv))
        results[kind] = (label, np.array(Ls), np.array(Ks), np.array(Dbox), np.array(Dth), np.array(N))
        savedata(f"cantor_{kind}", ["level", "L", "n_intervals", "lzma_bits_above_baseline", "D_box", "D_theory_level"],
                 [np.arange(1, len(Ls) + 1), Ls, N, Ks, Dbox, Dth],
                 [label, f"raster 2^22 pixels; LZMA preset 9|extreme; empty-bitmap baseline {base} bits subtracted"])
        print(f"  {label}")
        for i in range(len(Ls)):
            print(f"     level {i + 1:2d}: L={Ls[i]:6.2f}  N={N[i]:7d}  LZMA={Ks[i]:8d} bits  "
                  f"D_box={Dbox[i]:.3f}  D_theory(level)={Dth[i]:.3f}")

    # Chaitin barrier (illustrative): description length of the 7910-state ZFC-independent TM
    cF = 7910 * 2 * (1 + 1 + 13)
    D4 = np.nanmean(results["random_node"][4][-3:])
    L_star = np.log(cF) / D4
    print(f"\nIllustrative Chaitin constant c_F ~ {cF} bits (7910-state TM, ~15 bits per transition).")
    print(f"For the random-per-interval set (D~{D4:.2f}), K ~ e^(D L) exceeds c_F at L* ~ {L_star:.1f},"
          f" i.e. eps* ~ e^-L* = {np.exp(-L_star):.1e} of the system size.")
    print(f"Planck cutoff for a 1 m system: L_P = ln(1 m / 1.6e-35 m) = {np.log(1 / 1.616e-35):.1f}")

    # ---------------- figures
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.4), constrained_layout=True)
    a0 = ax[0]
    ss = [r[0] for r in rows]
    a0.plot(ss, [r[4] for r in rows], "o-", label="lower bound (normalised within budget)")
    a0.axhline(upper, color="tab:red", ls="--", label="upper bound (undecided + Kraft remainder)")
    a0.set_xlabel("max program length n (bits)"); a0.set_ylabel("Omega_BLC")
    a0.set_title(f"Toy halting probability: {nc} bit(s) certified\n"
                 "(certifying bits needs proofs of non-halting, not more computation)")
    a0.legend(fontsize=8)

    a1 = ax[1]
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
    for (kind, _), c in zip(kinds, colors):
        label, Ls, Ks, Dbox, Dth, N = results[kind]
        a1.semilogy(Ls, Ks, "o-", color=c, label=label + "  [LZMA]")
        k = np.arange(1, len(Ls) + 1)
        theory = {"selfsimilar": np.log2(k) + 1, "omega_dim": np.log2(k) + 1,
                  "random_level": k.astype(float), "random_node": N.astype(float)}[kind]
        a1.semilogy(Ls, theory, ":", color=c, lw=1.5)
    a1.axhline(cF, color="purple", lw=1.5)
    a1.text(1.2, cF * 1.3, "Chaitin barrier c_F (illustrative, ~2e5 bits)", color="purple", fontsize=8)
    a1.set_xlabel("L = ln(1/eps)"); a1.set_ylabel("bits above empty-bitmap baseline")
    a1.set_title("Solid: LZMA (upper bound on K).  Dotted: K up to O(1):\n"
                 "log2 k (computable), k (bit per level), N_eps (bit per interval)")
    a1.legend(fontsize=7, loc="center right")

    a2 = ax[2]
    for (kind, _), c in zip(kinds, colors):
        label, Ls, Ks, Dbox, Dth, N = results[kind]
        a2.plot(Ls, Dth, "-", color=c, label=label.split(":")[0])
        a2.plot(Ls, Dbox, "x", color=c)
    D_lo, D_hi = 0.25 + lower / 2, 0.25 + upper / 2
    a2.axhspan(D_lo, D_hi, color="tab:orange", alpha=0.2, label="D = 1/4 + Omega_BLC/2 (certified interval)")
    a2.set_xlabel("L = ln(1/eps)"); a2.set_ylabel("dimension")
    a2.set_title("Running dimension (lines: theory per level; x: box counting)")
    a2.legend(fontsize=7)
    savefig(fig, "chaitin_fractal")
    savedata("chaitin_barrier", ["c_F_bits", "D_random_node", "L_star", "eps_star", "L_Planck_1m"],
             [[cF], [D4], [L_star], [np.exp(-L_star)], [np.log(1 / 1.616e-35)]],
             ["illustrative c_F = 7910 states x 2 transitions x ~15 bits (Yedidia-Aaronson 2016)"])
    print(f"figures written to {FIG_DIR}/, data to {DATA_DIR}/")
    if a.show:
        plt.show()


if __name__ == "__main__":
    main()

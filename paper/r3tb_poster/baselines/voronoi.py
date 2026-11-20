"""WF weight in atomic domains of a MODEL density (proxy for a periodic LOBA).

PROXY, not the real method: neither real-space Wannier functions nor a charge
density are available to a plugin. Both are modelled from the wout/EOSCAR data:

  WF density   isotropic Gaussian at the WF centre, per-axis variance spread/3,
               summed over periodic images, on a uniform periodic grid.
  domains      CFG["domain"]:
    nearest    no grid: the two nuclei nearest to the WF centre, weight of the
               Gaussian on either side of their bisector plane (analytic erf).
    voronoi    Voronoi cells of the nuclei (grid labels, nearest nucleus).
    bader      steepest-ascent (on-grid, 26 neighbours) basins of the model
               valence density sum_n occ_n G_n(r) (occ_n = #WF occupation when
               the block is complete, else n_pair), optionally plus a pseudo-
               core at each nucleus: one isotropic Gaussian of per-axis width
               CFG["core_sigma"] (A) holding CFG["core_q"] electrons
               ("zval": ZVAL of the element, "unit": 1, "none": no core). The
               core only makes every nucleus an attractor of the density, as
               the real core/cusp does; it carries no WF and no electrons in
               the count.
               Every attractor is polished by mean-shift on the analytic
               density and its basin given to the nucleus nearest to it
               (equidistant nuclei share the basin).
    potential  steepest-ascent basins of the bare ion-core electrostatic
               potential of point charges ZVAL (uniform background; exact Ewald
               sum: smooth part by FFT, erfc real-space part added near every
               nucleus): the zero-flux partition of the ion-core field; one
               basin per nucleus.
  ownership    a WF belongs to the atom with the largest weight; margin =
               first - second weight; atoms of the top atom's element within
               TIE_GAP of the top weight share the WF equally (only when every
               atom inside that window is of that element).
  NELECT cap   WFs are full (n_pair electrons) or empty. With a complete #WF
               block the least-occupied WFs are emptied first. Without it,
               CFG["surplus"] = "shape": the surplus is taken from atom-centred
               WFs (offset < PARITY * rms radius), donor charges equalised
               (minimum sum q^2), the rest from the most diffuse WFs;
               "diffuse": most diffuse WFs first.
  integers     largest remainders, electron total conserved, ties by atom order.

No electronegativity, radius or oxidation-state table is read; ZVAL comes from
the EOSCAR.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.spatial import cKDTree
from scipy.special import erfc

CFG = {'domain': 'voronoi', 'core_q': 'none', 'core_sigma': 0.2, 'tie_gap': 0.2, 'surplus': 'shape', 'h0': 0.1, 'nmax': 3000000}

PARITY = 1.0 / 8.0      # centred WF: offset below PARITY * rms radius (value of frozen r3-T-B)
FLOAT_EPS = 1e-6        # numerical equality
CUT_SIGMA = 6.0         # Gaussian evaluated out to this many per-axis sigma
RHO_VOID = 1e-6         # e/A^3: grid points below this are not labelled (vacuum)
DIST_TIE = 1e-3         # A: nuclei this close to equidistant from an attractor share its basin
SMEAR_H = 2.0           # potential variant: Ewald splitting width in grid steps (result independent of it)
TWO_PI_32 = (2.0 * math.pi) ** 1.5

META = {
    "name": "model-density domain ownership (proxy LOBA): " + str(CFG["domain"]),
    "config": CFG,
    "params": {
        "tie_gap": [CFG["tie_gap"], "weight gap below which same-element atoms share a WF "
                    "(0.2 = the 60/40 population threshold of LOBA); fixed a priori"],
        "core_q/core_sigma": [[CFG["core_q"], CFG["core_sigma"]],
                              "pseudo-core charge rule and Gaussian width (A); a priori, scanned"],
        "PARITY": [PARITY, "centred-WF test of the surplus rule; inherited from frozen r3-T-B"],
        "h0/nmax": [[CFG["h0"], CFG["nmax"]], "grid spacing target (A) and point cap; numerical"],
        "CUT_SIGMA/RHO_VOID/DIST_TIE/SMEAR_H/FLOAT_EPS":
            [[CUT_SIGMA, RHO_VOID, DIST_TIE, SMEAR_H, FLOAT_EPS], "numerical"],
    },
    "table_free": True,
    "proxy": True,
}

SHIFTS = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float)


# ---------------------------------------------------------------- setup
def _nice(n):
    """Smallest integer >= n with prime factors 2, 3, 5 only (fast FFT)."""
    n = max(int(n), 4)
    while True:
        m = n
        for p in (2, 3, 5):
            while m % p == 0:
                m //= p
        if m == 1:
            return n
        n += 1


def _setup(case):
    st = case.structure
    L = np.array(st.lattice, float)
    invL = np.linalg.inv(L)
    fa = np.array(st.frac, float)
    fa = fa - np.floor(fa)
    els = list(st.elements)
    zval = np.array([float(case.eoscar.zval[e]) for e in els])
    n_pair = 2 if case.ispin == 1 else 1
    occ_map = {(w.spin, w.index): w.occ for w in (case.eoscar.wf_detail or [])}
    cen, spr, chn, keys = [], [], [], []
    for ci, ch in enumerate(case.channels):
        for wi, (c, s) in enumerate(zip(ch.centres, ch.spreads)):
            cen.append(c)
            spr.append(s)
            chn.append(ci)
            keys.append((ci + 1, wi + 1))
    has_occ = len(keys) > 0 and all(k in occ_map for k in keys)
    if has_occ:
        occ = np.array([min(max(occ_map[k], 0.0), float(n_pair)) for k in keys])
    else:
        occ = np.full(len(keys), float(n_pair))
    fc = np.array(cen, float).reshape(-1, 3) @ invL
    fc = fc - np.floor(fc)
    vol = abs(float(np.linalg.det(L)))
    h = max(CFG["h0"], (vol / CFG["nmax"]) ** (1.0 / 3.0))
    N = [_nice(math.ceil(np.linalg.norm(L[k]) / h)) for k in range(3)]
    hmax = max(np.linalg.norm(L[k]) / N[k] for k in range(3))
    spread = np.maximum(np.array(spr, float), 3.0 * hmax * hmax)   # sigma >= one grid step
    return {"L": L, "invL": invL, "fa": fa, "els": els, "na": len(els), "zval": zval,
            "n_pair": n_pair, "occ": occ, "has_occ": has_occ, "fc": fc, "chn": chn,
            "keys": keys, "spread_raw": np.array(spr, float), "sigma": np.sqrt(spread / 3.0),
            "vol": vol, "ngrid": N, "hmax": hmax,
            "nelect": case.eoscar.nelect_wann}


def _min_image(df, L):
    """Minimum-image Cartesian vectors and squared lengths for fractional
    differences df (..., 3)."""
    df = df - np.round(df)
    best_v = df @ L
    best_d2 = (best_v * best_v).sum(-1)
    for s in SHIFTS:
        if not s.any():
            continue
        v = (df + s) @ L
        d2 = (v * v).sum(-1)
        m = d2 < best_d2
        best_v[m] = v[m]
        best_d2[m] = d2[m]
    return best_v, best_d2


# ---------------------------------------------------------------- grid
class _Grid:
    def __init__(self, L, N):
        self.L = L
        self.N = [int(n) for n in N]
        self.invL = np.linalg.inv(L)
        self.g = L @ L.T
        self.colnorm = np.linalg.norm(self.invL, axis=0)
        self.size = self.N[0] * self.N[1] * self.N[2]
        self.vol = abs(float(np.linalg.det(L)))
        self._g2 = None
        self._m = None

    def box(self, fc, radius):
        """Index box covering the sphere of the given radius around fc.
        Returns None when the box is so large that the FFT route is cheaper."""
        half = [int(math.ceil(radius * self.colnorm[k] * self.N[k])) for k in range(3)]
        if (2 * half[0] + 1) * (2 * half[1] + 1) * (2 * half[2] + 1) > 2 * self.size:
            return None
        idx, d, dup = [], [], False
        for k in range(3):
            i0 = int(math.floor(fc[k] * self.N[k] + 0.5))
            ii = np.arange(i0 - half[k], i0 + half[k] + 1)
            dup = dup or ii.size > self.N[k]
            idx.append(ii % self.N[k])
            d.append(ii / self.N[k] - fc[k])
        x, y, z = d[0][:, None, None], d[1][None, :, None], d[2][None, None, :]
        g = self.g
        r2 = (g[0, 0] * x * x + g[1, 1] * y * y + g[2, 2] * z * z
              + 2.0 * (g[0, 1] * x * y + g[0, 2] * x * z + g[1, 2] * y * z))
        return idx, r2, dup

    def recip(self):
        if self._g2 is None:
            m = [np.fft.fftfreq(n) * n for n in self.N]
            iL = self.invL
            g2 = 0.0
            for j in range(3):
                comp = (iL[j, 0] * m[0][:, None, None] + iL[j, 1] * m[1][None, :, None]
                        + iL[j, 2] * m[2][None, None, :])
                g2 = g2 + comp * comp
            self._g2 = (2.0 * math.pi) ** 2 * g2
            self._m = m
        return self._g2, self._m

    def periodic_gauss(self, fc, sigma):
        """Unit-charge periodic Gaussian on the whole grid (Fourier route)."""
        g2, m = self.recip()
        ph = (np.exp(-2j * math.pi * m[0] * fc[0])[:, None, None]
              * np.exp(-2j * math.pi * m[1] * fc[1])[None, :, None]
              * np.exp(-2j * math.pi * m[2] * fc[2])[None, None, :])
        c = np.exp(-0.5 * sigma * sigma * g2) * ph
        return np.fft.ifftn(c).real * (self.size / self.vol)

    def add_gauss(self, rho, fc, sigma, q):
        b = self.box(fc, CUT_SIGMA * sigma)
        if b is None:
            rho += q * self.periodic_gauss(fc, sigma)
            return
        idx, r2, dup = b
        vals = (q / (TWO_PI_32 * sigma ** 3)) * np.exp(-r2 / (2.0 * sigma * sigma))
        if not dup:
            rho[np.ix_(idx[0], idx[1], idx[2])] += vals
        else:
            flat = ((idx[0][:, None, None] * self.N[1] + idx[1][None, :, None]) * self.N[2]
                    + idx[2][None, None, :])
            rho += np.bincount(flat.ravel(), weights=vals.ravel(),
                               minlength=self.size).reshape(self.N)

    def add_point_tail(self, field, fc, sigma, q):
        """Add q * erfc(r / (sqrt(2) sigma)) / r: point-charge potential minus
        the potential of the Gaussian-smeared charge (exact Ewald real-space part)."""
        idx, r2, dup = self.box(fc, CUT_SIGMA * sigma)
        r = np.sqrt(np.maximum(r2, (1e-3 * sigma) ** 2))
        vals = q * erfc(r / (math.sqrt(2.0) * sigma)) / r
        if not dup:
            field[np.ix_(idx[0], idx[1], idx[2])] += vals
        else:
            flat = ((idx[0][:, None, None] * self.N[1] + idx[1][None, :, None]) * self.N[2]
                    + idx[2][None, None, :])
            field += np.bincount(flat.ravel(), weights=vals.ravel(),
                                 minlength=self.size).reshape(self.N)

    def domain_weights(self, labels, ndom, fc, sigma):
        """Integral of the unit Gaussian over each domain (bin ndom = unlabelled)."""
        b = self.box(fc, CUT_SIGMA * sigma)
        if b is None:
            vals = self.periodic_gauss(fc, sigma)
            w = np.bincount(labels.ravel(), weights=vals.ravel(), minlength=ndom + 1)
            return w * (self.vol / self.size)
        idx, r2, _dup = b
        vals = np.exp(-r2 / (2.0 * sigma * sigma))
        lab = labels[np.ix_(idx[0], idx[1], idx[2])]
        w = np.bincount(lab.ravel(), weights=vals.ravel(), minlength=ndom + 1)
        return w * (self.vol / self.size) / (TWO_PI_32 * sigma ** 3)

    def cart_points(self, mask):
        f = [np.arange(n) / n for n in self.N]
        ii = np.nonzero(mask.ravel())[0]
        i0, rem = np.divmod(ii, self.N[1] * self.N[2])
        i1, i2 = np.divmod(rem, self.N[2])
        fr = np.stack([f[0][i0], f[1][i1], f[2][i2]], axis=1)
        return ii, fr


def _ascent(field, grid, active):
    """On-grid steepest ascent: every active grid point follows the neighbour
    (26 directions) of largest positive slope until a local maximum. Returns
    the flat index of the maximum reached (own index for inactive points)."""
    N = grid.N
    flat = np.arange(grid.size).reshape(N)
    best = np.zeros(N)
    ptr = flat.copy()
    step = np.array([[1.0 / N[0], 0, 0], [0, 1.0 / N[1], 0], [0, 0, 1.0 / N[2]]]) @ grid.L
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                if di == 0 and dj == 0 and dk == 0:
                    continue
                dist = float(np.linalg.norm(di * step[0] + dj * step[1] + dk * step[2]))
                sh = (-di, -dj, -dk)
                s = (np.roll(field, sh, axis=(0, 1, 2)) - field) / dist
                m = s > best
                best[m] = s[m]
                ptr[m] = np.roll(flat, sh, axis=(0, 1, 2))[m]
    ptr = ptr.ravel()
    act = active.ravel()
    ptr[~act] = np.nonzero(~act)[0]
    while True:
        nxt = ptr[ptr]
        if np.array_equal(nxt, ptr):
            break
        ptr = nxt
    return ptr


def _polish(x0, gc, gs, ga, L, invL):
    """Mean-shift ascent of the analytic Gaussian-mixture density from the
    Cartesian points x0 (m, 3). gc: fractional centres, gs: sigmas, ga: charges."""
    if len(x0) == 0 or len(gc) == 0:
        return x0
    out = np.array(x0, float)
    amp = ga / (TWO_PI_32 * gs ** 3) / (gs * gs)
    chunk = max(1, int(1e6 // max(1, len(gc))))
    for a in range(0, len(out), chunk):
        x = out[a:a + chunk]
        fx = x @ invL
        v, _ = _min_image(gc[None, :, :] - fx[:, None, :], L)      # (m, n, 3) from x to centre
        c = x[:, None, :] + v                                      # image-resolved centres
        for _ in range(2000):
            d2 = ((x[:, None, :] - c) ** 2).sum(-1)
            w = amp[None, :] * np.exp(-d2 / (2.0 * gs[None, :] ** 2))
            tot = w.sum(1)
            ok = tot > 0.0
            new = x.copy()
            new[ok] = (w[ok, :, None] * c[ok]).sum(1) / tot[ok, None]
            moved = float(np.abs(new - x).max()) if len(x) else 0.0
            x = new
            if moved < 1e-8:
                break
        out[a:a + chunk] = x
    return out


def _nearest_nuclei(x, g):
    """Share matrix (points x atoms): the nucleus nearest to each Cartesian
    point; nuclei within DIST_TIE of the nearest share equally."""
    fx = x @ g["invL"]
    _, d2 = _min_image(g["fa"][None, :, :] - fx[:, None, :], g["L"])
    d = np.sqrt(d2)
    share = (d <= d.min(1, keepdims=True) + DIST_TIE).astype(float)
    return share / share.sum(1, keepdims=True)


# ---------------------------------------------------------------- domains
def _model_density(g, grid, with_core):
    rho = np.zeros(grid.N)
    for n in range(len(g["fc"])):
        if g["occ"][n] > 0.0:
            grid.add_gauss(rho, g["fc"][n], g["sigma"][n], g["occ"][n])
    if with_core:
        q = _core_charges(g)
        for a in range(g["na"]):
            if q[a] > 0.0:
                grid.add_gauss(rho, g["fa"][a], CFG["core_sigma"], q[a])
    return rho


def _core_charges(g):
    if CFG["core_q"] == "zval":
        return g["zval"].copy()
    if CFG["core_q"] == "unit":
        return np.ones(g["na"])
    return np.zeros(g["na"])


def _weights_nearest(g):
    """Analytic two-nucleus Voronoi weights (no grid)."""
    L = g["L"]
    nw, na = len(g["fc"]), g["na"]
    W = np.zeros((nw, na))
    if nw == 0:
        return W
    v, d2 = _min_image(g["fa"][None, :, :] - g["fc"][:, None, :], L)   # (nw, na, 3)
    for n in range(nw):
        o = np.argsort(d2[n], kind="stable")
        a = o[0]
        if na == 1:
            W[n, a] = 1.0
            continue
        b = o[1]
        ab = v[n, b] - v[n, a]
        lab = float(np.linalg.norm(ab))
        # signed distance of the WF centre from the a-b bisector plane (a side > 0)
        delta = (d2[n, b] - d2[n, a]) / (2.0 * lab) if lab > 0 else 0.0
        wa = 0.5 * (1.0 + math.erf(delta / (math.sqrt(2.0) * g["sigma"][n])))
        W[n, a], W[n, b] = wa, 1.0 - wa
    return W


def _weights_grid(g, info):
    grid = _Grid(g["L"], g["ngrid"])
    domain = CFG["domain"]
    na = g["na"]
    with_core = domain == "bader" and CFG["core_q"] != "none"
    rho = _model_density(g, grid, with_core)
    active = rho > RHO_VOID
    info["grid"] = list(grid.N)
    info["active_fraction"] = float(active.mean())

    if domain == "voronoi":
        ii, fr = grid.cart_points(active)
        img = ((g["fa"][None, :, :] + SHIFTS[:, None, :]).reshape(-1, 3)) @ g["L"]
        owner = np.tile(np.arange(na), len(SHIFTS))
        _, k = cKDTree(img).query(fr @ g["L"])
        labels = np.full(grid.size, na, dtype=np.int64)
        labels[ii] = owner[k]
        ndom, M = na, np.eye(na)
    else:
        if domain == "potential":
            chg = np.zeros(grid.N)
            sg = SMEAR_H * g["hmax"]
            for a in range(na):
                grid.add_gauss(chg, g["fa"][a], sg, g["zval"][a])
            g2, _m = grid.recip()
            ck = np.fft.fftn(chg)
            ck[0, 0, 0] = 0.0
            g2s = g2.copy()
            g2s[0, 0, 0] = 1.0
            field = np.fft.ifftn(4.0 * math.pi * ck / g2s).real
            for a in range(na):
                grid.add_point_tail(field, g["fa"][a], sg, g["zval"][a])
            act = np.ones(grid.N, dtype=bool)
        else:
            field, act = rho, active
        ptr = _ascent(field, grid, act)
        roots, inv = np.unique(ptr[act.ravel()], return_inverse=True)
        labels = np.full(grid.size, len(roots), dtype=np.int64)
        labels[act.ravel()] = inv
        ndom = len(roots)
        i0, rem = np.divmod(roots, grid.N[1] * grid.N[2])
        i1, i2 = np.divmod(rem, grid.N[2])
        x = np.stack([i0 / grid.N[0], i1 / grid.N[1], i2 / grid.N[2]], axis=1) @ g["L"]
        if domain == "bader":
            gc, gs, ga = g["fc"], g["sigma"], g["occ"]
            keep = ga > 0.0
            gc, gs, ga = gc[keep], gs[keep], ga[keep]
            if with_core:
                q = _core_charges(g)
                gc = np.vstack([gc, g["fa"]])
                gs = np.concatenate([gs, np.full(na, CFG["core_sigma"])])
                ga = np.concatenate([ga, q])
            x = _polish(x, gc, gs, ga, g["L"], g["invL"])
        M = _nearest_nuclei(x, g) if ndom else np.zeros((0, na))
        info["n_attractors"] = int(ndom)
        info["attractor_xyz"] = x
        info["attractor_share"] = M
    labels = labels.reshape(grid.N)
    nw = len(g["fc"])
    W = np.zeros((nw, na))
    captured = np.zeros(nw)
    centre_dom = []
    for n in range(nw):
        w = grid.domain_weights(labels, ndom, g["fc"][n], g["sigma"][n])
        captured[n] = w[:ndom].sum()
        W[n] = w[:ndom] @ M
        i = [int(math.floor(g["fc"][n][k] * grid.N[k] + 0.5)) % grid.N[k] for k in range(3)]
        centre_dom.append(int(labels[i[0], i[1], i[2]]))
    info["captured"] = captured
    # diagnostic only: atoms owning the domain that contains the WF centre
    info["centre_hosts"] = [sorted(int(a) for a in np.nonzero(M[d] > 0)[0]) if d < ndom else []
                            for d in centre_dom]
    return W


# ---------------------------------------------------------------- ownership
def _ownership(g, W):
    els, na = g["els"], g["na"]
    Wn = W.copy()
    tot = Wn.sum(1)
    fallback = _weights_nearest(g)
    hosts, margin, top2 = [], [], []
    for n in range(len(Wn)):
        row = Wn[n] / tot[n] if tot[n] > 1e-12 else fallback[n]
        o = np.argsort(-row, kind="stable")
        a = int(o[0])
        w1 = float(row[a])
        w2 = float(row[o[1]]) if na > 1 else 0.0
        window = [int(b) for b in o if w1 - row[b] < CFG["tie_gap"]]
        if len(window) > 1 and all(els[b] == els[a] for b in window):
            hosts.append({b: 1.0 / len(window) for b in window})
        else:
            hosts.append({a: 1.0})
        margin.append(w1 - w2)
        top2.append((a, w1, int(o[1]) if na > 1 else -1, w2))
        Wn[n] = row
    return hosts, margin, top2, Wn


def _equalise(q, cap, members, amount):
    """Remove `amount` electrons from `members`, raising their charges to a
    common level (minimum sum of q^2); atom a gives at most cap[a]."""
    if amount <= 0.0 or not members:
        return {}

    def removed(level):
        return sum(min(max(level - q[a], 0.0), cap[a]) for a in members)

    knots = sorted(set([q[a] for a in members] + [q[a] + cap[a] for a in members]))
    lev = knots[-1]
    for lo, hi in zip(knots, knots[1:]):
        r_hi = removed(hi)
        if r_hi >= amount:
            r_lo = removed(lo)
            lev = lo + (hi - lo) * (amount - r_lo) / (r_hi - r_lo) if r_hi > r_lo else hi
            break
    return {a: min(max(lev - q[a], 0.0), cap[a]) for a in members}


def _drain(order_keys, left, hosts, amount, lost):
    """Empty WFs tier by tier (equal keys within FLOAT_EPS form a tier, emptied
    equally per electron held). order_keys: list of (key, wf index), ascending."""
    i = 0
    while amount > FLOAT_EPS and i < len(order_keys):
        j = i
        while j < len(order_keys) and abs(order_keys[j][0] - order_keys[i][0]) <= FLOAT_EPS:
            j += 1
        tier = [n for _k, n in order_keys[i:j]]
        have = sum(left[n] for n in tier)
        if have > FLOAT_EPS:
            take = min(amount, have)
            for n in tier:
                f = take * left[n] / have
                for a, s in hosts[n].items():
                    lost[a] += f * s
                left[n] -= f
            amount -= take
        i = j
    return amount


def _populations(g, hosts, info):
    na, n_pair, nw = g["na"], g["n_pair"], len(hosts)
    ne = np.zeros(na)
    for h in hosts:
        for a, s in h.items():
            ne[a] += n_pair * s
    lost = np.zeros(na)
    info["mode"] = "closed"
    info["excess"] = 0
    if g["nelect"] is None:
        return ne, lost
    excess = float(n_pair * nw - int(round(g["nelect"])))
    info["excess"] = excess
    if excess <= 0:
        return ne, lost
    left = [float(n_pair)] * nw
    if g["has_occ"]:
        info["mode"] = "occupations"
        keys = sorted((float(g["occ"][n]) / n_pair, n) for n in range(nw))
        excess = _drain(keys, left, hosts, excess, lost)
    else:
        info["mode"] = "no occupations / " + CFG["surplus"]
        centred = []
        if CFG["surplus"] == "shape":
            _, d2 = _min_image(g["fa"][None, :, :] - g["fc"][:, None, :], g["L"])
            for n in range(nw):
                if len(hosts[n]) == 1:
                    a = next(iter(hosts[n]))
                    if math.sqrt(d2[n, a]) < PARITY * math.sqrt(max(g["spread_raw"][n], 0.0)):
                        centred.append((n, a))
            donors = sorted({a for _n, a in centred})
            cap = np.zeros(na)
            for n, a in centred:
                cap[a] += left[n]
            q = g["zval"] - ne
            x = _equalise(q, cap, donors, excess)
            for a, v in x.items():
                lost[a] += v
                excess -= v
            for n, _a in centred:
                left[n] = 0.0        # exhausted whenever a remainder is left
            info["donors"] = donors
        if excess > FLOAT_EPS:
            keys = sorted((-float(g["spread_raw"][n]), n) for n in range(nw) if left[n] > FLOAT_EPS)
            excess = _drain(keys, left, hosts, excess, lost)
    return ne - lost, lost


def _apportion(ne):
    """Largest remainders with the total conserved; ties by atom order."""
    total = int(round(float(ne.sum())))
    base = [int(math.floor(x + FLOAT_EPS)) for x in ne]
    rem = [x - b for x, b in zip(ne, base)]
    rest = total - sum(base)
    order = sorted(range(len(ne)), key=lambda a: (-round(rem[a], 6), a))
    out = list(base)
    k = 0
    while rest > 0 and order:
        out[order[k % len(order)]] += 1
        rest -= 1
        k += 1
    order_low = sorted(range(len(ne)), key=lambda a: (round(rem[a], 6), -a))
    k = 0
    while rest < 0 and order_low:
        out[order_low[k % len(order_low)]] -= 1
        rest += 1
        k += 1
    return out


def analyse(case):
    """Everything: weights, hosts, margins, populations, integers."""
    g = _setup(case)
    info = {}
    # numpy/Accelerate raises spurious floating-point flags in matmul on this
    # platform; results are checked for finiteness instead.
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        if CFG["domain"] == "nearest":
            W = _weights_nearest(g)
        else:
            W = _weights_grid(g, info)
    if not np.isfinite(W).all():
        raise FloatingPointError("non-finite domain weights")
    hosts, margin, top2, Wn = _ownership(g, W)
    ne, lost = _populations(g, hosts, info)
    n_int = _apportion(ne)
    os_ = [int(round(z)) - k for z, k in zip(g["zval"], n_int)]
    info.update({"g": g, "weights": Wn, "hosts": hosts, "margin": margin, "top2": top2,
                 "populations": ne, "lost": lost, "os": os_})
    return info


def assign(case):
    return analyse(case)["os"]

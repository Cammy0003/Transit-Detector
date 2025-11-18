Got it. Here’s a **lean, resume-friendly** pipeline that detects transits in TESS PDCSAP light curves **without** the heavy spline detrender. It’s accurate enough to find obvious candidates, easy to reason about end-to-end, and small enough to implement cleanly in Rust or Python.

---

# A simple transit-detection recipe (no fancy detrending)

## 0) Inputs & hygiene (one screen of code)

* Use columns: **TIME**, **PDCSAP_FLUX**, **PDCSAP_FLUX_ERR**, filter **QUALITY == 0**.
* Sort by time, drop NaNs/Infs.
* Let ($\Delta t$) = median cadence (robust).
* **Segment on gaps**: start a new segment whenever ($t_{i}-t_{i-1} > 5,\Delta t$). This prevents gap edges from screwing with the search.

## 1) Minimal normalization (exactly one line per segment)

For each segment (s):
$$
F^{(s)}_{\text{norm}}(t) = \frac{F^{(s)}_{\text{raw}}(t)}{\operatorname{median}\big(F^{(s)}_{\text{raw}}\big)} - 1
$$
That’s it. Now your baseline is ~0 and “depth” is a **negative fraction** (e.g., ($-500\ \text{ppm} \approx -5\times10^{-4})$).

Optional but helpful (still simple): **clip wild outliers** with a one-sided rule to protect dips:

* Estimate noise via MAD: ($\sigma \approx 1.4826 \cdot \operatorname{MAD}(F^{(s)}_{\text{norm}})$).
* Replace points where ($F^{(s)}_{\text{norm}} > +6\sigma$) with the segment median (ignore negative side to keep transits intact).

> That’s your whole “detrend”: **none**. PDCSAP already removes common trends; segment + normalize + clip is enough for a resume project.

## 2) Box Least Squares (BLS) on a fixed grid (the core)

You’re fitting a simple **box model** to a **folded** light curve for each trial period. Very small math footprint, easy to implement.

**Grids**

* Period: ($P \in [0.5, 20]$) days (or whatever you want)
* Duration: ($d \in {0.5, 1, 2, 3, 4}$) hours (convert to days)
* Phase: continuous in ($[0, P)$) but we’ll find phase efficiently (below)

**Fold**
For a trial ($P$), compute phases ($\phi_i = \big(t_i \bmod P\big)/P \in [0,1)$).

**Fast in-transit search with prefix sums**

1. Bin phases into ($B$) bins (e.g., $(B=200)$), keep:

    * mean flux per bin ($y_b$) and weights ($w_b = 1/\sigma_b^2$) (or equal weights if you want it simpler).
2. For a trial duration ($d$) (in phase units ($q = d/P)$), you want the **contiguous arc** ($[\phi_0, \phi_0+q)$) whose average flux is most negative.

    * On the circular array ($y_b$), use a **sliding window** over ($L=\lfloor qB \rfloor$) bins with prefix sums to find the minimum weighted average in ($O(B)$).
3. The “depth” estimate for that arc is simply the mean inside window minus the mean outside:
   $$
   \hat{\delta} ;=; \frac{\sum_{\text{in}} w_b y_b}{\sum_{\text{in}} w_b} ;-; \frac{\sum_{\text{out}} w_b y_b}{\sum_{\text{out}} w_b}
   $$
   It should be **negative** for transits.

**Power / SNR**
Estimate per-point noise (\sigma) from the **out-of-window** bins (MAD → ($\sigma$)). Then a simple score:
$$
\text{SNR}(P,d) = \frac{-\hat{\delta}}{\sigma}\sqrt{N_{\text{in}}}
$$
where ($N_{\text{in}}$) is the number of points (or total weight) inside the best window. Keep the **best SNR** over phases (found by the sliding window) for each ($(P,d)$), and then keep the best over ($d$) for each ($P$). Your **periodogram** is ($P \mapsto \max_d \text{SNR}(P,d)$).

This is the essence of BLS, implemented with basic arrays and prefix sums. No external libs required, no spline anything.

## 3) Pick the candidate & refine

* Select the **top peak** in the BLS periodogram: ($(P^*, d^*, \phi^*)$) with max SNR.
* **Fine-tune** around ($P^*$) with a smaller step (e.g., ±2% range, 5–10× finer resolution). Re-estimate ($(P,d,\phi)$).
* Report:

    * ($P$) (days), duration ($d$) (hours), depth ($\delta$) (ppm), SNR
    * Number of observed transits (how many windows land on data)
    * A quick **phase-folded plot** (bin and show in-transit vs out-of-transit)

## 4) Sanity checks (still simple, no new methods)

* **Odd/even**: Fold on ($P$) and ($2P$). If the two dips at (2P) differ in depth by ($\gtrsim 3\sigma$), it smells like an eclipsing binary → flag.
* **Secondary search** at phase ~0.5 on ($P$). If there’s a shallower “secondary” that’s still strong → flag as EB.
* **Depth stability** across segments: compute depth per segment; if wildly inconsistent → flag (likely systematics).

## 5) Minimal outputs (good for a resume)

* One CSV row per target with: TIC (if you have it), best ($P$), ($d$), ($\delta$), SNR, #transits, flags (odd/even, secondary, depth-stability).
* Plots:

    * Periodogram ( $\text{SNR}(P)$ )
    * Phase-folded light curve at ($P^*$) with the best-fit box overlaid
* A tiny README: “what I did, what each figure means, caveats.”

---

## Why this works (and why it’s simple)

* **PDCSAP** already knocked out most shared trends. You don’t have to fight the star’s full systematics budget.
* The only “trend handling” you do is **segment + normalize + outlier clip**—all of which you can explain in two lines of math.
* BLS is just “find the darkest repeating box after folding,” which you can implement with **prefix sums and a sliding window** (1 page of Rust).

---

## Pseudocode you can drop into Rust

```rust
// 0) Load & hygiene
let (t, f, fe) = load_tess_pdcsap(file)?;
let mask = quality == 0 && f.is_finite() && fe.is_finite();
let (t, f, fe) = filter_and_sort(t[mask], f[mask], fe[mask]);

// segment on gaps
let dt_med = median_diff(&t)?;
let segments = segment_on_gaps(&t, 5.0 * dt_med);

// 1) Normalize per segment & clip positive outliers
let mut t_all = Vec::new();
let mut y_all = Vec::new();
for (s, e) in segments {
    let fs = &f[s..e];
    let ts = &t[s..e];
    let med = median(fs)?;
    let mut ys: Vec<f64> = fs.iter().map(|&v| v/med - 1.0).collect();
    let sigma = 1.4826 * mad(&ys)?;
    for y in ys.iter_mut() {
        if *y > 6.0 * sigma { *y = 0.0; } // recenter outliers to baseline
    }
    t_all.extend_from_slice(ts);
    y_all.extend_from_slice(&ys);
}

// 2) BLS period scan (prefix sums on a phase histogram)
let periods = linspace(0.5, 20.0, 4000); // days
let durations_hr = [0.5, 1.0, 2.0, 3.0, 4.0];
let nbins = 200;

let mut best = None; // (snr, P, d_hr, phi0, depth)
for &P in &periods {
    let phase = fold_phases(&t_all, P); // [0,1)
    let (yb, wb) = bin_phase(&phase, &y_all, nbins); // averages & weights
    let (prefix_yw, prefix_w) = prefix_sums(&yb, &wb); // for O(1) window sums

    for &d_hr in &durations_hr {
        let q = (d_hr / 24.0) / P; // phase duration
        let L = (q * nbins as f64).max(1.0) as usize;
        let (phi0_bin, depth, snr) = best_window_on_circle(&prefix_yw, &prefix_w, L);
        if best.map_or(true, |(bs, _,_,_,_)| snr > bs) {
            let phi0 = (phi0_bin as f64) / (nbins as f64);
            best = Some((snr, P, d_hr, phi0, depth));
        }
    }
}

if let Some((snr, p, d_hr, phi0, depth)) = best {
    println!("BEST: P={:.5} d, dur={:.2} h, depth={:.0} ppm, SNR={:.1}, phi0={:.3}",
             p, d_hr, depth*1e6, snr, phi0);
}
```

> Every function in there is short: `segment_on_gaps`, `median_diff`, `mad`, `fold_phases`, `bin_phase`, `prefix_sums`, `best_window_on_circle`. Nothing exotic.

---

## What to write on the resume

* **“Implemented a from-scratch Box Least Squares transit searcher for TESS PDCSAP light curves in Rust.”**
* Segmentation & normalization only (no heavyweight detrending); **O(N)** prefix-sum BLS with phase binning.
* Outputs: periodogram, best period/duration/depth/SNR, phase-folded plot; flags for odd/even & secondary eclipses.
* (Optional) Compared against a known planet to demonstrate it recovers the published period (include the figure).

---

## If you want *one* optional upgrade that stays simple

Add a **fixed, single-parameter high-pass**: a rolling median with a **fixed window** of, say, **1.5 days** per segment,
$$
B(t) = \operatorname{rolling_median}\big(F_{\text{raw}},, W=1.5\ \text{days}\big),\quad
F_{\text{hp}} = \frac{F_{\text{raw}}}{B}-1.
$$
Don’t grid-search it; just fix (W). If you do this, keep the rest identical. It’s still explainable in two lines and gains you robustness on slowly varying stars.

---

If you want, I can turn the pseudocode above into a clean Rust module (no external math crates), plus a tiny CLI that ingests a FITS file and spits out the CSV + plots.
